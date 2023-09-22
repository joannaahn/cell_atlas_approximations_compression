# vim: fdm=indent
'''
author:     Fabio Zanini
date:       22/09/23
content:    Compress isodiametra pulchra, an acoel.
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
    )


species = 'i_pulchra'
atlas_data_folder = root_repo_folder / 'data' / 'full_atlases' / 'RNA' / species
#anno_fn = root_repo_folder / 'data' / 'gene_annotations' / 'l_minuta??z'
fn_out = output_folder / f'{species}.h5'


rename_dict = {
    'cell_types': {
        'neurons': 'neuron',
        'secretory cells': 'secretory',
        'cilia related': 'ciliated',
        'sensory neurons': 'neuron',
        'stem cells': 'stem',
        'epithelial i': 'epithelial',
        'epithelial ii': 'epithelial',
        'digestive i': 'digestive',
        'digestive ii': 'digestive',
        'trp+ neurons': 'neuron',
    }
}

celltype_tissue_blacklist = {
        'whole': [
            'uncharacterized',
            'contamination',
        ],
}


coarse_cell_types = []


celltype_order = [
    ('epithelial', [
        'epithelial',
        'ciliated',
    ]),
    ('mesenchymal', [
        'digestive',
        'muscle',
    ]),
    ('other', [
        'stem',
        'neuron',
        'chemosensory',
        'secretory',
    ]),
]


if __name__ == '__main__':

    # Remove existing compressed atlas file if present
    if os.path.isfile(fn_out):
        os.remove(fn_out)

    compressed_atlas = {}

    tissues = ['whole']
    for tissue in tissues:
        adata_tissue = anndata.read(atlas_data_folder / (species + '.h5ad'))

        # It's already in raw counts

        # cptt throughout
        sc.pp.normalize_total(
            adata_tissue,
            target_sum=1e4,
            key_added='coverage',
        )

        # Fix cell type annotations
        adata_tissue.obs['cell_type'] = adata_tissue.obs['manual_annotations'].str.lower()
        adata_tissue.obs['cellType'] = fix_annotations(
            adata_tissue, 'cell_type', species, tissue,
            rename_dict, coarse_cell_types,
            blacklist=celltype_tissue_blacklist,
        )


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

        compressed_atlas[tissue] = {
            'features': genes,
            'celltype': {
                'avg': avg_ge,
                'frac': frac_ge,
                'ncells': ncells_ge,
            },
        }

    print('No gene annotations available')
    gene_annos = None

    print('Store compressed atlas to file')
    store_compressed_atlas(
            fn_out,
            compressed_atlas,
            tissues,
            gene_annos,
            celltype_order,
    )
