# vim: fdm=indent
'''
author:     Fabio Zanini
date:       22/09/23
content:    Compress isodiametra pulchra, an acoel.
'''
import os
import sys
import gc
import pathlib
import gzip
import h5py
import numpy as np
import pandas as pd

import anndata
import scanpy as sc

from utils import (
    load_config,
    root_repo_folder,
    output_folder,
    get_tissue_data_dict,
    subannotate,
    normalise_counts,
    correct_annotations,
    get_celltype_order,
    compress_tissue,
    collect_feature_annotations,
    collect_feature_sequences,
    store_compressed_feature_sequences,
    store_compressed_atlas,
    filter_cells,
    )


if __name__ == '__main__':

    species_list = [
        # Multi-organ species
        'h_sapiens',
        'm_musculus',
        'm_myoxinus',
        'd_melanogaster',
        'x_laevis',

        # Single-organ species
        'l_minuta',
        'h_miamia',
        'a_queenslandica',
        'c_elegans',
        'i_pulchra',
        'a_queenslandica',
        'd_rerio',
        't_adhaerens',
        's_mediterranea',
        's_mansoni',
        's_lacustris',
        'm_leidyi',
        'n_vectensis',
    ]

    for species in species_list:
        print('--------------------------------')
        print(species)
        print('--------------------------------')

        config = load_config(species)
        atlas_data_folder = root_repo_folder / 'data' / 'full_atlases' / 'RNA' / species
        fn_out = output_folder / f'{species}.h5'

        # Remove existing compressed atlas file if present
        if os.path.isfile(fn_out):
            os.remove(fn_out)

        # Iterate over gene expression, chromatin accessibility, etc.
        for measurement_type in config["measurement_types"]:
            # FIXME
            if measurement_type == "gene_expression":
                continue

            compressed_atlas = {}
            config_mt = config[measurement_type]
            celltype_order = config_mt["cell_annotations"]["celltype_order"]

            load_params = {}
            if 'load_params' in config_mt:
                load_params.update(config_mt['load_params'])

            if "path_global" in config_mt:
                print(f"Read full atlas")
                adata = anndata.read(config_mt["path_global"], **load_params)

            if "path_metadata_global" in config_mt:
                print("Read global metadata separately")
                meta = pd.read_csv(config_mt["path_metadata_global"], sep='\t', index_col=0).loc[adata.obs_names]

                if ("filter_cells_global" in config_mt) and ("metadata" in config_mt["filter_cells_global"]):
                    column, value = config_mt["filter_cells_global"]["metadata"]
                    meta = meta.loc[meta[column] == value]

                if 'tissues' in config_mt['cell_annotations']['rename_dict']:
                    tissues_raw = meta['tissue'].value_counts().index.tolist()
                    tdict = config_mt['cell_annotations']['rename_dict']['tissues']
                    tmap = {t: tdict.get(t, t) for t in tissues_raw}
                    meta['tissue'] = meta['tissue'].map(tmap)
                    del tdict, tmap

                tissues = meta['tissue'].value_counts().index.tolist()
                tissues = sorted([t for t in tissues if t != ''])
            else:
                tissues = config_mt["tissues"]

            # Iterate over tissues
            for tissue in tissues:
                if tissue == '':
                    continue
                print(tissue)

                if "path_metadata_global" in config_mt:
                    meta_tissue = meta.loc[meta['tissue'] == tissue]

                if "path_global" not in config_mt:
                    print(f"Read full atlas for {tissue}")
                    adata_tissue = anndata.read(config_mt["path"][tissue], **load_params)
                else:
                    print(f'Slice cells for {tissue}')
                    if "path_metadata_global" in config_mt:
                        adata_tissue = adata[meta_tissue.index]
                    else:
                        adata_tissue = adata[adata.obs['tissue'] == tissue]

                if ("load_params" in config_mt) and ("backed" in config_mt["load_params"]):
                    adata_tissue = adata_tissue.to_memory()
                
                if "path_metadata_global" in config_mt:
                    adata_tissue.obs = meta_tissue.copy()

                print("Filter cells")
                adata_tissue = filter_cells(adata_tissue, config_mt["filter_cells"])

                print("Normalise")
                adata_tissue = normalise_counts(
                    adata_tissue,
                    config_mt['normalisation'],
                    measurement_type,
                )

                print("Correct cell annotations")
                adata_tissue = correct_annotations(
                    adata_tissue,
                    config_mt['cell_annotations']['column'],
                    species,
                    tissue,
                    config_mt['cell_annotations']['rename_dict'],
                    config_mt['cell_annotations']['require_subannotation'],
                    blacklist=config_mt['cell_annotations']['blacklist'],
                    subannotation_kwargs=config_mt['cell_annotations']['subannotation_kwargs'],
                )

                print("Compress atlas")
                compressed_atlas[tissue] = compress_tissue(
                    adata_tissue, celltype_order,
                )

            # TODO: harmonise across tissues, sometimes (e.g. fly) that is not a given
            features = adata_tissue.var_names

            del adata_tissue
            if "path_global" in config_mt:
                del adata
            if "path_metadata_global" in config_mt:
                del meta

            print('Garbage collection before storing compressed atlas')
            gc.collect()

            print('Store compressed atlas')
            store_compressed_atlas(
                    fn_out,
                    compressed_atlas,
                    tissues,
                    celltype_order,
                    measurement_type=measurement_type,
            )

            del compressed_atlas
            del tissues
            del celltype_order

            if "feature_sequences" in config_mt:
                print('Garbage collection before storing feature sequences')
                gc.collect()

                print('Collect feature sequences')
                feature_sequences = collect_feature_sequences(
                    config_mt,
                    features,
                    measurement_type, species,
                )

                print('Store feature sequences')
                store_compressed_feature_sequences(
                    fn_out,
                    feature_sequences,
                    measurement_type,
                )

            del feature_sequences

            # FIXME
            if False:
                print('Garbage collection before storing feature annotations')
                gc.collect()

                print('Collect feature annotations')
                feature_annos = collect_feature_annotations(
                        config_mt['feature_annotation'],
                        features,
                        measurement_type,
                )

                if feature_annos is not None:
                    print('Store feature annotations')
                    store_compressed_feature_annotations(
                        fn_out,
                        feature_annos,
                        measurement_type,
                    )

                del feature_annos

            print('Garbage collection at the end of a species and measurement type')
            gc.collect()

            # FIXME
            # OK, now onwards to ATAC-Seq
            sys.exit()
