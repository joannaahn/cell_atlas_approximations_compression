# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/05/22
content:    Compress C Elegans (Cao et al 2017)
'''
import os
import sys
import pathlib
import gzip
import h5py
import numpy as np
import pandas as pd

import anndata
import scanpy as sc

from utils import (
    root_repo_folder,
    output_folder,
    get_tissue_data_dict,
    subannotate,
    fix_annotations,
    get_celltype_order,
    collect_gene_annotations,
    store_compressed_atlas,
    sanitise_gene_names,
    )


species = 's_lacustris'
atlas_data_folder = root_repo_folder / 'data' / 'full_atlases' / 'RNA' / 's_lacustris'
# There cannot be an annotation yet since the transcriptome was assembled here
#anno_fn = root_repo_folder / 'data' / 'gene_annotations' / 'c_elegans.PRJNA13758.WS287.annotations.gff3.gz'
fn_out = output_folder / f'{species}.h5'

# TODO
rename_dict = {
    'cell_types': {
        'myopeptidocytes1': 'myopeptidocyte',
        'myopeptidocytes2': 'myopeptidocyte',
        'mesocytes 1': 'mesocyte',
        'mesocytes 2': 'mesocyte',
        'mesocytes 3': 'mesocyte',
        'metabolocytes1': 'metabolocyte',
        'metabolocytes2': 'metabolocyte',
        'choanoblasts1': 'choanoblast',
        'choanoblasts2': 'choanoblast',
        'incurrent pinacocytes1': 'incurrent pinacocyte',
        'incurrent pinacocytes2': 'incurrent pinacocyte',
        'apendopinacocytes1': 'apendopinacocyte',
        'apendopinacocytes2': 'apendopinacocyte',
        'archaeocytes': 'archaeocyte',
        'choanocytes': 'choanocyte',
        'apopylar cells': 'apopylar',
        'sclerocytes': 'sclerocyte',
        'granulocytes': 'granulocyte',
        'basopinacocytes': 'basopinacocyte',
        'lophocytes': 'lophocyte',
        'amoebocytes': 'amoebocyte',
        'sclerophorocytes': 'sclerophorocyte',
    },
}

coarse_cell_types = [
]


celltype_order = [
    ('immune', [
        'granulocyte',
        'amoebocyte',
    ]),
    ('epithelial', [
        'choanocyte',
        'apopylar',
        'choanoblast',
        'incurrent pinacocyte',
        'basopinacocyte',
        'apendopinacocyte',
    ]),
    ('endothelial', [
    ]),
    ('mesenchymal', [
        'archaeocyte',
        'myopeptidocyte',
        'metabolocyte',
        'mesocyte',
        'sclerocyte',
        'sclerophorocyte',
        'lophocyte',
    ]),
    ('other', [
        'neuroid',
    ]),
    ('unknown', [
        'unknown',
    ])
]


if __name__ == '__main__':

    # Remove existing compressed atlas file if present
    if os.path.isfile(fn_out):
        os.remove(fn_out)

    compressed_atlas = {}

    tissue_sources = get_tissue_data_dict(species, atlas_data_folder)
    tissues = list(tissue_sources.keys())
    for it, (tissue, full_atlas_fn) in enumerate(tissue_sources.items()):
        print(tissue)

        adata_tissue = anndata.read(full_atlas_fn)

        # Ignore cells with NaN in the cell.type column
        idx = adata_tissue.obs['cell_type'].isin(
                adata_tissue.obs['cell_type'].value_counts().index)
        adata_tissue = adata_tissue[idx]

        # cptt throughout
        sc.pp.normalize_total(
            adata_tissue,
            target_sum=1e4,
            key_added='coverage',
        )

        # Fix cell type annotations
        adata_tissue.obs['cell_type'] = adata_tissue.obs['cell_type'].str.lower()
        adata_tissue.obs['cellType'] = fix_annotations(
            adata_tissue, 'cell_type', 's_lacustris', tissue,
            rename_dict, coarse_cell_types,
        )

        # Age
        adata_tissue.obs['age'] = '8day'

        # Correction might declare some cells as untyped/low quality
        # they have an empty string instead of an actual annotation
        if (adata_tissue.obs['cellType'] == '').sum() > 0:
            idx = adata_tissue.obs['cellType'] != ''
            adata_tissue = adata_tissue[idx]

        celltypes = get_celltype_order(
            adata_tissue.obs['cellType'].value_counts().index,
            celltype_order,
        )

        print('Add data to celltype group')
        genes = adata_tissue.var_names
        avg_ge = pd.DataFrame(
                np.zeros((len(genes), len(celltypes)), np.float32),
                index=genes,
                columns=celltypes,
                )
        frac_ge = pd.DataFrame(
                np.zeros((len(genes), len(celltypes)), np.float32),
                index=genes,
                columns=celltypes,
                )
        ncells_ge = pd.Series(
                np.zeros(len(celltypes), np.int64), index=celltypes,
                )
        for celltype in celltypes:
            idx = adata_tissue.obs['cellType'] == celltype
            Xidx = adata_tissue[idx].X
            avg_ge[celltype] = np.asarray(Xidx.mean(axis=0))[0]
            frac_ge[celltype] = np.asarray((Xidx > 0).mean(axis=0))[0]
            ncells_ge[celltype] = idx.sum()

        print('Add data to celltype-timepoint group')
        ages = adata_tissue.obs['age'].value_counts().index.tolist()
        ages.sort()
        columns_age = []
        for ct in celltypes:
            for age in ages:
                columns_age.append('_'.join([ct, 'Slacustris', str(age)]))

        # Averages
        genes = adata_tissue.var_names
        avg_ge_tp = pd.DataFrame(
                np.zeros((len(genes), len(celltypes) * len(ages)), np.float32),
                index=genes,
                columns=columns_age,
                )
        frac_ge_tp = pd.DataFrame(
                np.zeros((len(genes), len(celltypes) * len(ages)), np.float32),
                index=genes,
                columns=columns_age,
                )
        ncells_ge_tp = pd.Series(
                np.zeros(len(columns_age), np.int64), index=columns_age,
                )
        for celltype in celltypes:
            adata_ct = adata_tissue[adata_tissue.obs['cellType'] == celltype]
            for age in ages:
                idx_age = (adata_ct.obs['age'] == age).values.nonzero()[0]
                if len(idx_age) == 0:
                    continue
                Xct_age = adata_ct.X[idx_age]
                label = '_'.join([celltype, 'Musser', str(age)])
                avg_ge_tp[label] = np.asarray(Xct_age.mean(axis=0))[0]
                frac_ge_tp[label] = np.asarray((Xct_age > 0).mean(axis=0))[0]
                ncells_ge_tp[label] = len(idx_age)

        genes = sanitise_gene_names(genes)

        compressed_atlas[tissue] = {
            'features': np.asarray(genes),
            'celltype': {
                'avg': avg_ge,
                'frac': frac_ge,
                'ncells': ncells_ge,
            },
            'celltype_dataset_timepoint': {
                'avg': avg_ge_tp,
                'frac': frac_ge_tp,
                'ncells': ncells_ge_tp,
            },
        }

    print('Consolidate gene list across tissues')
    needs_union = False
    genes = None
    for tissue, tdict in compressed_atlas.items():
        genest = list(tdict['features'])
        if genes is None:
            genes = genest
            continue
        if genest != genes:
            needs_union = True
            genes = set(genes) | set(genest)

    if needs_union:
        raise NotImplementedError('TODO: make union of features')

    #print('Get gene annotations')
    #gene_annos = collect_gene_annotations(anno_fn, genes)
    gene_annos = None  # Not available yet

    print('Store compressed atlas to file')
    store_compressed_atlas(
            fn_out,
            compressed_atlas,
            tissues,
            gene_annos,
            celltype_order,
    )
