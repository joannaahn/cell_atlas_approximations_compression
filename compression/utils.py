'''
Utility functions for the compression
'''
import os
import gc
import pathlib
import yaml
import gzip
import numpy as np
import pandas as pd
import h5py

import scanpy as sc


root_repo_folder = pathlib.Path(__file__).parent.parent
# Try to put the output in the API repo if available
output_folder = root_repo_folder / '..' / 'cell_atlas_approximations_API' / 'web' / 'static' / 'atlas_data'
if not output_folder.is_dir():
    output_folder = root_repo_folder / 'data' / 'atlas_approximations'


def load_config(species):
    config_path = pathlib.Path(__file__).parent / "organism_configs" / (species + ".yml")
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Propagate cell supertypes and order for simplicity
    if 'cell_annotations' in config:
        for mt in config['measurement_types']:
            config_mt = config[mt]
            if 'cell_annotations' not in config_mt:
                config_mt['cell_annotations'] = config['cell_annotations']
            else:
                for st in config['cell_annotations']:
                    config_mt['cell_annotations'][st] = config['cell_annotations'][st]

    for mt in config["measurement_types"]:
        config_mt = config[mt]

        if ("path" not in config_mt) and ("path_global" not in config_mt):
            config_mt["path"] = species + '.h5ad'

        if ("path" in config_mt) and isinstance(config_mt["path"], str):
            config_mt["path"] = {t: config_mt["path"] for t in config_mt["tissues"]}

        # Use absolute paths
        root_fdn = root_repo_folder / 'data' / 'full_atlases' / mt / species
        for key in ["path_global", "path_metadata_global"]:
            if key in config_mt:
                config_mt[key] = root_fdn / config_mt[key]
        if "path" in config_mt:
            for tissue in config_mt["path"]:
                config_mt["path"][tissue] = root_fdn / config_mt["path"][tissue]

        if "filter_cells" not in config_mt:
            config_mt["filter_cells"] = {}

        if "require_subannotation" not in config_mt["cell_annotations"]:
            config_mt["cell_annotations"]["require_subannotation"] = []

        if "subannotation_kwargs" not in config_mt["cell_annotations"]:
            config_mt["cell_annotations"]["subannotation_kwargs"] = {}

        if "blacklist" not in config_mt["cell_annotations"]:
            config_mt["cell_annotations"]["blacklist"] = {}

        celltype_order = []
        for supertype in config_mt["cell_annotations"]["supertype_order"]:
            ct_order_supertype = (
                supertype, config_mt["cell_annotations"]["cell_supertypes"][supertype],
            )
            celltype_order.append(ct_order_supertype)
        config_mt["cell_annotations"]["celltype_order"] = celltype_order

        del config_mt["cell_annotations"]["supertype_order"]
        del config_mt["cell_annotations"]["cell_supertypes"]

        if "feature_annotation" not in config_mt:
            config_mt["feature_annotation"] = False
        else:
            # FIXME
            config_mt["feature_annotation"] = root_repo_folder / 'data' / 'gene_annotations' / config_mt["feature_annotation"]

    return config


def filter_cells(adata, config_mt):
    """Filter cells according to some parameter dictionary."""
    filter_dict = config_mt.get("filter_cells", {})

    if len(filter_dict) == 0:
        return adata

    ncells_orig = adata.shape[0]

    if "min_cells_per_type" in filter_dict:
        nmin = filter_dict["min_cells_per_type"]

        column = config_mt["cell_annotations"]["column"]
        ct_counts = adata.obs[column].value_counts()
        ct_abundant = ct_counts.index[ct_counts >= nmin]
        adata = adata[adata.obs[column].isin(ct_abundant)]

    if "unannotated" in filter_dict:
        unanno_values = filter_dict["unannotated"]
        if isinstance(unanno_values, str):
            unanno_values = [unanno_values]

        column = config_mt["cell_annotations"]["column"]
        # If a list is given, take the first one you find or fail
        if not isinstance(column, str):
            columns = column
            for column in columns:
                if column in adata.obs.columns:
                    break
            else:
                raise ValueError(
                    f"None of the cell annotation columns found: {columns}",
                )
        adata = adata[~adata.obs[column].isin(unanno_values)]

    ncells_new = adata.shape[0]

    if ncells_new < ncells_orig:
        delta = ncells_orig - ncells_new
        print(f'Filtered out {delta} cells, originally {ncells_orig} cells, {ncells_new} remaining')

    return adata


def subannotate(adata,
                species, annotation,
                markers,
                bad_prefixes=None,
                verbose=True,
                trash_unknown=True,
                skip_subannotation=False):
    '''This function subannotates a coarse annotation from an atlasi.

    This is ad-hoc, but that's ok for now. Examples are 'lymphocyte', which is
    a useless annotation unless you know what kind of lymphocytes these are, or
    if it's a mixed bag.
    '''
    # If skipping, return list of empty annotations - basically blacklisting
    if skip_subannotation:
        return [""] * adata.shape[0]

    if bad_prefixes is None:
        bad_prefixes = []

    markersi = markers.get(annotation, None)
    if markersi is None:
        raise ValueError(
            f'Cannot subannotate without markers for {species}, {annotation}')

    adata = adata.copy()
    sc.pp.log1p(adata)

    genes, celltypes = [], []
    for celltype, markers_ct in markersi.items():
        celltypes.append(celltype)
        for gene in markers_ct:
            if gene in adata.var_names:
                genes.append(gene)
            elif verbose:
                print('Missing gene:', gene)

    adatam = adata[:, genes].copy()

    # No need for PCA because the number of genes is small

    # Get neighbors
    sc.pp.neighbors(adatam)

    # Get communities
    sc.tl.leiden(adatam)

    adata.obs['subleiden'] = adatam.obs['leiden']
    sc.tl.rank_genes_groups(
        adata,
        'subleiden',
        method='t-test_overestim_var',
    )
    top_marker = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(2)

    subannos = {}
    for cluster, genestop in top_marker.items():
        found = False
        for gene in genestop:
            if found:
                break
            found_bad_prefix = False
            for bad_pfx in bad_prefixes:
                if gene.startswith(bad_pfx):
                    found_bad_prefix = True
                    break
            if found_bad_prefix:
                subannos[cluster] = ''
                continue
            for celltype, markers_ct in markersi.items():
                if gene in markers_ct:
                    subannos[cluster] = celltype
                    found = True
                    break
            else:
                # FIXME: trash clusters with unknown markers for now
                if not trash_unknown:
                    import ipdb; ipdb.set_trace()
                    raise ValueError('Marker not found:', gene)
                else:
                    subannos[cluster] = ''
        if not found:
            subannos[cluster] = ''

    new_annotations = adata.obs['subleiden'].map(subannos)

    return new_annotations


def correct_annotations(
    adata,
    column,
    species,
    tissue,
    rename_dict,
    require_subannotation,
    blacklist=None,
    subannotation_kwargs=None,
):
    '''Correct cell types in each tissue according to known dict'''
    # If a list is given, take the first one you find or fail
    if not isinstance(column, str):
        columns = column
        for column in columns:
            if column in adata.obs.columns:
                break
        else:
            raise ValueError(
                f"None of the cell annotation columns found: {columns}",
            )

    # Ignore cells with NaN in the cell.type column
    idx = adata.obs[column].isin(
            adata.obs[column].value_counts().index)
    adata = adata[idx].copy()

    gc.collect()

    adata.obs[column + '_lowercase'] = adata.obs[column].str.lower()

    if blacklist is None:
        blacklist = {}
    if subannotation_kwargs is None:
        subannotation_kwargs = {}

    celltypes_new = np.asarray(adata.obs[column + '_lowercase']).copy()

    # Exclude blacklisted
    if tissue in blacklist:
        for ctraw in blacklist[tissue]:
            celltypes_new[celltypes_new == ctraw] = ''

    # Rename according to standard dict
    if 'cell_types' in rename_dict:
        for ctraw, celltype in rename_dict['cell_types'].items():
            # one can use brain:neuron for renaming in specific tissues only
            if isinstance(ctraw, str) and (':' not in ctraw):
                celltypes_new[celltypes_new == ctraw] = celltype
            else:
                # Organ-specific renames
                if isinstance(ctraw, str) and ':' in ctraw:
                    organraw, ctraw = ctraw.split(':')
                else:
                    organraw, ctraw = ctraw
                if organraw == tissue:
                    celltypes_new[celltypes_new == ctraw] = celltype

    ct_found = np.unique(celltypes_new)

    # In some data sets, some unnotated clusters are denoted by a digit
    for ctraw in ct_found:
        if ctraw.isdigit():
            celltypes_new[celltypes_new == ctraw] = ''
    ct_found = np.unique(celltypes_new)

    # Look for coarse annotations
    ctnew_list = set(celltypes_new)
    for celltype in ctnew_list:
        if celltype in require_subannotation:
            idx = celltypes_new == celltype
            adata_coarse_type = adata[idx]
            print(f'Subannotating {celltype}')
            subannotations = subannotate(
                adata_coarse_type, species, celltype,
                **subannotation_kwargs,
            )

            # Ignore reclustering into already existing types, we have enough
            for subanno in subannotations:
                if subanno in ct_found:
                    subannotations[subannotations == subanno] = ''
            print('Subannotation done')

            celltypes_new[idx] = subannotations

    adata.obs['cellType'] = celltypes_new

    # Eliminate cell types with less than 3 cells
    ncells = adata.obs['cellType'].value_counts()
    rare_celltypes = ncells.index[ncells < 3]
    adata.obs.loc[adata.obs['cellType'].isin(rare_celltypes), 'cellType'] = ''

    # Correction might declare some cells as untyped/low quality
    # they have an empty string instead of an actual annotation
    if (adata.obs['cellType'] == '').sum() > 0:
        idx = adata.obs['cellType'] != ''
        adata= adata[idx]

    return adata


def get_celltype_order(celltypes_unordered, celltype_order):
    '''Use global order to reorder cell types for this tissue'''
    celltypes_ordered = []
    for broad_type, celltypes_broad_type in celltype_order:
        for celltype in celltypes_broad_type:
            if celltype in celltypes_unordered:
                celltypes_ordered.append(celltype)

    celltypes_found = []
    missing_celltypes = False
    for celltype in celltypes_unordered:
        if celltype not in celltypes_ordered:
            if not missing_celltypes:
                missing_celltypes = True
                print('Missing celltypes:')
            print(celltype)
        else:
            celltypes_found.append(celltype)
    if missing_celltypes:
        print('Cell types found:')
        for celltype in celltypes_found:
            print(celltype)

    if missing_celltypes:
        raise IndexError("Missing cell types!")

    return celltypes_ordered


def collect_gene_annotations(anno_fn, genes):
    '''Collect gene annotations from GTF file'''
    if anno_fn is False:
        return None

    featype = 'gene'

    with gzip.open(anno_fn, 'rt') as gtf:
        gene_annos = []
        for line in gtf:
            if f'\t{featype}\t' not in line:
                continue
            fields = line.split('\t')
            if fields[2] != featype:
                continue
            attrs = fields[-1].split(';')

            gene_name = None
            transcript_id = None
            for attr in attrs:
                if 'gene_name' in attr:
                    gene_name = attr.split(' ')[-1][1:-1]
                elif " gene " in attr:
                    gene_name = attr.split(' ')[-1][1:-1]
                elif "ID=gene:" in attr:
                    gene_name = attr.split(';')[0].split(':')[-1]
                # for drosophila, flybase uses "gene_symbol"
                elif 'gene_symbol' in attr:
                    gene_name = attr.split(" ")[-1][1:-1]
                elif 'transcript_id' in attr:
                    transcript_id = attr.split(' ')[-1][1:-1]
                elif 'Name=' in attr:
                    gene_name = attr.split('=')[1]
                    transcript_id = gene_name

            if (gene_name is not None) and (transcript_id is None):
                transcript_id = gene_name

            if (gene_name is None) or (transcript_id is None):
                continue
            gene_annos.append({
                'transcript_id': transcript_id,
                'gene_name': gene_name,
                'chromosome_name': fields[0],
                'start_position': int(fields[3]),
                'end_position': int(fields[4]),
                'strand': 1 if fields[6] == '+' else -1,
                'transcription_start_site': int(fields[3]) if fields[6] == '+' else int(fields[4]),
                })
    gene_annos = pd.DataFrame(gene_annos)

    # NOTE: some species like zebrafish don't really have a transcript id (yet?)
    #assert gene_annos['transcript_id'].value_counts()[0] == 1

    # FIXME: choose the largest transcript or something. For this particular
    # repo it's not that important
    gene_annos = (gene_annos.drop_duplicates('gene_name')
                            .set_index('gene_name', drop=False))

    genes_missing = list(set(genes) - set(gene_annos['gene_name'].values))
    gene_annos_miss = pd.DataFrame([], index=genes_missing)
    gene_annos_miss['transcript_id'] = gene_annos_miss.index
    gene_annos_miss['start_position'] = -1
    gene_annos_miss['end_position'] = -1
    gene_annos_miss['strand'] = 0
    gene_annos_miss['chromosome_name'] = ''
    gene_annos_miss['transcription_start_site'] = -1
    gene_annos = pd.concat([gene_annos, gene_annos_miss])
    gene_annos = gene_annos.loc[genes]
    gene_annos['strand'] = gene_annos['strand'].astype('i2')

    return gene_annos


def collect_feature_annotations(anno_fn, features, measurement_type="gene_expression"):
    if measurement_type == "gene_expression":
        return collect_gene_annotations(anno_fn, features)
    return None


def store_compressed_feature_sequences(
    fn_out,
    feature_sequences,
    measurement_type,
    compression=19,
    ):
    '''Store compressed features into h5 file.'''
    # Optional zstd compression using hdf5plugin
    if compression:
        import hdf5plugin
        # NOTE: decompressing zstd is equally fast no matter how much compression.
        # As for compression speed, levels 1-19 are normal, 20-22 "ultra".
        # A quick runtime test shows *faster* access for clevel=22 than clevel=3,
        # while the file size is around 10% smaller. Compression speed is significantly
        # slower, but (i) still somewhat faster than actually averaging the data and
        # (ii) compresses whole human RNA+ATAC is less than 1 minute. That's nothing
        # considering these approximations do not change that often.
        comp_kwargs = hdf5plugin.Zstd(clevel=compression)
    else:
        comp_kwargs = {}


    with h5py.File(fn_out, 'a') as h5_data:
        me = h5_data[measurement_type]
        group = me.create_group('feature_sequences')
        group.attrs["type"] = feature_sequences["type"]
        # Compression here is really important, even though slow. Decompression should be
        # fast no matter what. Doing more than level 19 requires a LOT of memory
        group.create_dataset(
            "sequences", data=feature_sequences["sequences"].values.astype('S'),
            **comp_kwargs,
        )


def store_compressed_feature_annotations(
    fn_out,
    feature_annos,
    measurement_type,
    ):
    '''Store compressed feature annotations.'''
    with h5py.File(fn_out, 'a') as h5_data:
        me = h5_data[measurement_type]

        group = me.create_group('feature_annotations')
        group.create_dataset(
                'feature_name',
                data=feature_annos.index.values.astype('S'))
        if 'transcription_start_site' in feature_annos.columns:
            group.create_dataset(
                    'transcription_start_site',
                    data=feature_annos['transcription_start_site'].values, dtype='i8')
        group.create_dataset(
                'chromosome_name',
                data=feature_annos['chromosome_name'].astype('S'))
        group.create_dataset(
                'start_position',
                data=feature_annos['start_position'].values, dtype='i8')
        group.create_dataset(
                'end_position',
                data=feature_annos['end_position'].values, dtype='i8')
        group.create_dataset(
                'strand', data=feature_annos['strand'].values, dtype='i2')


def store_compressed_atlas(
        fn_out,
        compressed_atlas,
        tissues,
        celltype_order,
        measurement_type='gene_expression',
        compression=22,
        quantisation="chromatin_accessibility",
        #chunked=True,
        ):
    '''Store compressed atlas into h5 file.

    Args:
        fn_out: The h5 file with the compressed atlas.
        compressed_atlas: The dict with the result.
        tissues: A list of tissues covered.
        celltype_order: The order of cell types.
        measurement_type: What type of data this is (gene expression, chromatin accessibility, etc.).
        quantisation: If not None, average measurement is quantised with these bins.
        compression: Use zstd compression of the data arrays (avg and frac). Levels are 1-22,
            whereas 0 or False means no compression. No performace decrease is observed.
    '''
    add_kwargs = {}

    # Optional zstd compression using hdf5plugin
    if compression:
        import hdf5plugin
        # NOTE: decompressing zstd is equally fast no matter how much compression.
        # As for compression speed, levels 1-19 are normal, 20-22 "ultra".
        # A quick runtime test shows *faster* access for clevel=22 than clevel=3,
        # while the file size is around 10% smaller. Compression speed is significantly
        # slower, but (i) still somewhat faster than actually averaging the data and
        # (ii) compresses whole human RNA+ATAC is less than 1 minute. That's nothing
        # considering these approximations do not change that often.
        comp_kwargs = hdf5plugin.Zstd(clevel=compression)
    else:
        comp_kwargs = {}

    # Data can be quantised for further compression (typically ATAC-Seq)
    if (quantisation == True) or (quantisation == measurement_type):
        if measurement_type == "chromatin_accessibility":
            # NOTE: tried using quantiles for this, but they are really messy
            # and subject to aliasing effects. 8 bits are more than enough for
            # most biological questions given the noise in the data
            qbits = 8
            bins = np.array([-0.001, 1e-8] + np.logspace(-4, 0, 2**qbits - 1).tolist()[:-1] + [1.1])
            # bin "centers", aka data quantisation
        elif measurement_type == "gene_expression":
            # Counts per ten thousand quantisation
            qbits = 16
            bins = np.array([-0.001, 1e-8] + np.logspace(-2, 4, 2**qbits - 1).tolist()[:-1] + [1.1e4])
        else:
            raise ValueError(f"Quantisation for {measurement_type} not set.")

        quantisation_array = [0] + np.sqrt(bins[1:-2] * bins[2:-1]).tolist() + [1]

        qbytes = qbits // 8
        # Add a byte if the quantisation is not optimal
        if qbits not in (8, 16, 32, 64):
            qbytes += 1
        avg_dtype = f"u{qbytes}"
        quantisation = True
    else:
        avg_dtype = "f4"
        quantisation = False

    for tissue, group in compressed_atlas.items():
        features = group['features'].tolist()
        break

    with h5py.File(fn_out, 'a') as h5_data:
        me = h5_data.create_group(measurement_type)
        me.create_dataset('features', data=np.array(features).astype('S'))
        if quantisation:
            me.create_dataset('quantisation', data=np.array(quantisation_array).astype('f4'))

        me.create_dataset('tissues', data=np.array(tissues).astype('S'))
        supergroup = me.create_group('by_tissue')
        for tissue in tissues:
            tgroup = supergroup.create_group(tissue)
            #for label in ['celltype', 'celltype_dataset_timepoint']:
            for label in ['celltype']:
                group = tgroup.create_group(label)

                # Number of cells
                ncells = compressed_atlas[tissue][label]['ncells']
                group.create_dataset(
                    'cell_count', data=ncells.values, dtype='i8')

                # Average in a cell type
                avg = compressed_atlas[tissue][label]['avg']
                if quantisation:
                    # pd.cut wants one dimensional arrays so we ravel -> cut -> reshape
                    avg_vals = (pd.cut(avg.values.ravel(), bins=bins, labels=False)
                                .reshape(avg.shape)
                                .astype(avg_dtype))
                    avg = pd.DataFrame(
                        avg_vals, columns=avg.columns, index=avg.index,
                    )

                # TODO: manual chunking might increase performance a bit, the data is
                # typically accessed only vertically (each feature its own island)
                #if chunked:
                #    # Chunk each feature on its own: this is perfect for ATAC-Seq 
                #    add_kwargs['chunks'] = (1, len(features))

                # Cell types
                group.create_dataset(
                    'index', data=avg.columns.values.astype('S'))
                group.create_dataset(
                    'average', data=avg.T.values, dtype=avg_dtype,
                    **add_kwargs,
                    **comp_kwargs,
                )
                if measurement_type == 'gene_expression':
                    # Fraction detected in a cell type
                    frac = compressed_atlas[tissue][label]['frac']
                    group.create_dataset(
                        'fraction', data=frac.T.values, dtype='f4',
                        **add_kwargs,
                        **comp_kwargs,
                    )

                # Local neighborhoods
                neid = compressed_atlas[tissue][label]['neighborhood']
                neigroup = group.create_group('neighborhood')
                ncells = neid['ncells']
                neigroup.create_dataset(
                    'cell_count', data=ncells.values, dtype='i8')
                avg = neid['avg']
                if quantisation:
                    # pd.cut wants one dimensional arrays so we ravel -> cut -> reshape
                    avg_vals = (pd.cut(avg.values.ravel(), bins=bins, labels=False)
                                .reshape(avg.shape)
                                .astype(avg_dtype))
                    avg = pd.DataFrame(
                        avg_vals, columns=avg.columns, index=avg.index,
                    )
                neigroup.create_dataset(
                    'index', data=avg.columns.values.astype('S'))
                neigroup.create_dataset(
                    'average', data=avg.T.values, dtype=avg_dtype,
                    **add_kwargs,
                    **comp_kwargs,
                )
                if measurement_type == 'gene_expression':
                    # Fraction detected in a cell type
                    frac = neid['frac']
                    neigroup.create_dataset(
                        'fraction', data=frac.T.values, dtype='f4',
                        **add_kwargs,
                        **comp_kwargs,
                    )

                # TODO: centroid coordinates and convex hulls


        ct_group = me.create_group('celltypes')
        supertypes = np.array([x[0] for x in celltype_order])
        ct_group.create_dataset(
                'supertypes',
                data=supertypes.astype('S'),
                )
        for supertype, subtypes in celltype_order:
            ct_group.create_dataset(
                supertype,
                data=np.array(subtypes).astype('S'),
            )


def sanitise_gene_names(genes):
    genes_new = []
    for gene in genes:
        gene_new = gene.replace(',', ';')
        gene_new = gene_new.split(' ')[0]
        genes_new.append(gene_new)

    if len(set(genes_new)) != len(set(genes)):
        raise ValueError("Gene names are not unique after sanitisation.")

    return genes_new


def normalise_counts(adata_tissue, input_normalisation, measurement_type="gene_expression"):
    """Normalise counts no matter what the input normalisation is."""
    if measurement_type == "gene_expression":
        if input_normalisation not in (
                "cptt", "raw", "cpm", "cpm+log", "cptt+log", "to-raw", "to-raw+cptt+log"):
            raise ValueError("Input normalisation not recognised: {input_normalisation}")

        if input_normalisation in ("to-raw", "to-raw+cptt+log"):
            adata_tissue = adata_tissue.raw.to_adata()

        if input_normalisation in ("cpm+log", "cptt+log", "to-raw+cptt+log"):
            adata_tissue.X = np.expm1(adata_tissue.X)

        if input_normalisation in ("raw", "cpm", "cpm+log", "to-raw"):
            sc.pp.normalize_total(
                adata_tissue,
                target_sum=1e4,
                key_added='coverage',
            )

        return adata_tissue

    elif measurement_type == "chromatin_accessibility":
        if input_normalisation not in ("binary", "to-binary"):
            raise ValueError("Input normalisation not recognised: {input_normalisation}")

        if input_normalisation == "to-binary":
            adata_tissue.X.data[:] = 1

        return adata_tissue

    raise ValueError("measurement type not recognised")


def compress_tissue(
    adata_tissue,
    celltype_order,
    measurement_type="gene_expression",
    max_neighborhoods=5,
):
    """Compress atlas for one tissue after data is clean, normalised, and reannotated."""
    celltypes = get_celltype_order(
        adata_tissue.obs['cellType'].value_counts().index,
        celltype_order,
    )

    features = adata_tissue.var_names
    avg = pd.DataFrame(
            np.zeros((len(features), len(celltypes)), np.float32),
            index=features,
            columns=celltypes,
            )
    if measurement_type == "gene_expression":
        frac = pd.DataFrame(
                np.zeros((len(features), len(celltypes)), np.float32),
                index=features,
                columns=celltypes,
                )
    ncells = pd.Series(
            np.zeros(len(celltypes), np.int64), index=celltypes,
            )
    for celltype in celltypes:
        idx = adata_tissue.obs['cellType'] == celltype
        
        # Number of cells
        ncells[celltype] = idx.sum()

        # Average across cell type
        Xidx = adata_tissue[idx].X
        avg[celltype] = np.asarray(Xidx.mean(axis=0))[0]
        if measurement_type == "gene_expression":
            frac[celltype] = np.asarray((Xidx > 0).mean(axis=0))[0]

    # Local neighborhoods
    neid = _compress_neighborhoods(
        ncells,
        adata_tissue,
        max_neighborhoods=max_neighborhoods,
        measurement_type=measurement_type,
    )


    result = {
        'features': features,
        'celltype': {
            'avg': avg,
            'neighborhood': {
                'avg': neid['avg'],
            },
            'ncells': ncells,
        },
    }
    if measurement_type == "gene_expression":
        result['celltype']['frac'] = frac
        result['celltype']['neighborhood']['frac'] = neid['frac']

    return result


def _compress_neighborhoods(
    ncells,
    adata,
    max_neighborhoods=5,
    max_cells_per_type=300,
    measurement_type='gene_expression',
):
    """Compress local neighborhood of a single cell type."""
    # Try something easy first, like k-means
    from sklearn.cluster import KMeans
    from scipy.spatial import ConvexHull

    features = adata.var_names

    celltypes = list(ncells.keys())
    nei_columns = []
    nei_ncells = pd.Series(
            np.zeros(len(celltypes) * max_neighborhoods, np.int64),
            )
    nei_avg = pd.DataFrame(
            np.zeros((len(features), len(celltypes) * max_neighborhoods), np.float32),
            index=features,
            )
    nei_coords = pd.DataFrame(
            np.zeros((2, len(celltypes) * max_neighborhoods), np.float32),
            index=['x', 'y'],
            )
    convex_hulls = []
    if measurement_type == "gene_expression":
        nei_frac = pd.DataFrame(
                np.zeros((len(features), len(celltypes) * max_neighborhoods), np.float32),
                index=features,
                )

    # Tune neighborhood number for rare cell types
    # NOTE: the following lines need to be in order, obviously
    n_neighborhoods = ncells.copy()
    n_neighborhoods[:] = max_neighborhoods
    #n_neighborhoods[ncells < 150] = 5
    n_neighborhoods[ncells < 75] = 4
    n_neighborhoods[ncells < 25] = 3

    # Subsample with some regard for cell typing
    cell_ids = []
    for celltype, ncell in ncells.items():
        cell_ids_ct = adata.obs_names[adata.obs['cellType'] == celltype]
        if ncell > max_cells_per_type:
            idx_rand = np.random.choice(range(ncell), size=max_cells_per_type, replace=False)
            cell_ids_ct = cell_ids_ct[idx_rand]
        cell_ids.extend(list(cell_ids_ct))
    adata = adata[cell_ids].copy()

    # Log
    sc.pp.log1p(adata)

    # Select features
    sc.pp.highly_variable_genes(adata)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]

    # Create embedding, a proxy for cell states broadly
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    points = adata.obsm['X_umap']

    # Back to all features for storage
    adata = adata.raw.to_adata()
    adata.obsm['X_umap'] = points

    # Back to cptt or equivalent for storage
    adata.X.data = np.expm1(adata.X.data)

    for celltype, n_nei in n_neighborhoods.items():
        adata_ct = adata[adata.obs['cellType'] == celltype]
        points_ct = adata_ct.obsm['X_umap']

        # Assign cells to mutually exclusive states (for now)
        # NOTE: reduce the K until all clusters have at least 3 cells
        for n_neii in range(n_nei, 0, -1):
            kmeans = KMeans(
                n_clusters=n_neii,
                random_state=0,
                n_init='auto',
            ).fit(points_ct) 

            ncells_ct = pd.Series(kmeans.labels_).value_counts()
            if ncells_ct.min() >= 3:
                break
        else:
            raise ValueError(f"Celltype with < 3 cells total: {celltype}")

        for i in range(kmeans.n_clusters):
            idx = kmeans.labels_ == i
            # Add the number of cells
            nei_ncells.iloc[:, len(nei_columns) + i] = idx.sum()

            # Add the average expression
            nei_avg.iloc[:, len(nei_columns) + i] = np.asarray(adata_ct.X[idx].mean(axis=0))[0]
            # Add the fraction expressing
            if measurement_type == "gene_expression":
                nei_frac.iloc[:, len(nei_columns) + i] = np.asarray((adata_ct.X[idx] > 0).mean(axis=0))[0]

            # Add the coordinates of the center
            points_i = points_ct[idx]
            nei_coords.iloc[:, len(nei_columns) + i] = points_i.mean(axis=0)

            # Add the convex hull
            hull = ConvexHull(points_i)
            convex_hulls.append(points_i[hull.vertices])

        # Housekeeping
        nei_columns.extend([celltype] * kmeans.n_clusters)

    # Clean up
    del adata
    gc.collect()

    nei_ncells = nei_ncells.iloc[:len(nei_columns)]
    nei_avg = nei_avg.iloc[:, :len(nei_columns)]
    nei_coords = nei_coords.iloc[:, :len(nei_columns)]
    nei_avg.columns = nei_columns
    nei_coords.columns = nei_columns
    if measurement_type == "gene_expression":
        nei_frac = nei_frac.iloc[:, :len(nei_columns)]
        nei_frac.columns = nei_columns

    neid = {
        'ncells': nei_ncells,
        'avg': nei_avg,
        'coords_centroid': nei_coords,
        'convex_hull': convex_hulls,
    }
    if measurement_type == "gene_expression":
        neid['frac'] = nei_frac

    return neid


def collect_feature_sequences(config_mt, features, measurement_type, species):
    """Collect sequences of features to enable cross-species comparisons."""
    # CREDIT NOTE: FROM BIOPYTHON
    def SimpleFastaParser(handle):
        """Iterate over Fasta records as string tuples.
    
        Arguments:
         - handle - input stream opened in text mode
    
        For each record a tuple of two strings is returned, the FASTA title
        line (without the leading '>' character), and the sequence (with any
        whitespace removed). The title line is not divided up into an
        identifier (the first word) and comment or description.
    
        >>> with open("Fasta/dups.fasta") as handle:
        ...     for values in SimpleFastaParser(handle):
        ...         print(values)
        ...
        ('alpha', 'ACGTA')
        ('beta', 'CGTC')
        ('gamma', 'CCGCC')
        ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
        ('delta', 'CGCGC')
    
        """
        # Skip any text before the first record (e.g. blank lines, comments)
        for line in handle:
            if line[0] == ">":
                title = line[1:].rstrip()
                break
        else:
            # no break encountered - probably an empty file
            return
    
        # Main logic
        # Note, remove trailing whitespace, and any internal spaces
        # (and any embedded \r which are possible in mangled files
        # when not opened in universal read lines mode)
        lines = []
        for line in handle:
            if line[0] == ">":
                yield title, "".join(lines).replace(" ", "").replace("\r", "")
                lines = []
                title = line[1:].rstrip()
                continue
            lines.append(line.rstrip())
    
        yield title, "".join(lines).replace(" ", "").replace("\r", "")

    path = config_mt['feature_sequences']['path']
    path = root_repo_folder / 'data' / 'full_atlases' / measurement_type / species / path

    seqs = {fea: "" for fea in features}
    with gzip.open(path, 'rt') as f:
        for gene, seq in SimpleFastaParser(f):
            # Sometimes they need a gene/id combo from biomart
            if '|' in gene:
                gene = gene.split('|')[0]
            if gene == '':
                continue

            if gene in seqs:
                seqs[gene] = seq

    seqs = pd.Series(seqs).loc[features]

    return {
        'sequences': seqs,
        'type': config_mt['feature_sequences']['type'],
    }
