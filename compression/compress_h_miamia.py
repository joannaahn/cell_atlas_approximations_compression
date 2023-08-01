# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/05/22
content:    Compress h miamia.
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


species = 'h_miamia'
dataset_name = 'BoWang'
atlas_data_folder = root_repo_folder / 'data' / 'full_atlases' / 'RNA' / species
#anno_fn = root_repo_folder / 'data' / 'gene_annotations' / 'Schistosoma_mansoni.Smansoni_v7.57.gff3.gz'
fn_out = output_folder / f'{species}.h5'

# TODO
rename_dict = {
    'tissues': {},
    'cell_types': {
        # Reannotated by Dania in Bo Wang's lab.
        # see original at https://www.nature.com/articles/s41467-023-38016-4
        'germ_cells_1': 'germline',
        'germ_cells_2': 'germline',
        'muscle': 'striated muscle',
        'endo_progen': 'endo progenitor',
        'neural': 'neuron',
        'neural-2': 'neuron',
        'epi_pro': 'epidermal progenitor',
        'neo-2': 'pluripotent',
        'endo_like-I': 'digestive progenitor',
        'endo_like-III': 'endo tail',
        'endo_like-II': 'diffuse endo',
    }
}

celltype_tissue_blacklist = {
    'whole': [
        # Unknowns, question marks
        'uk-1', 'uk-2', 'uk-3',
        'GSC?',
    ],
}

coarse_cell_types = [
]


celltype_order = [
    ('immune', [
    ]),
    ('epithelial', [
        'epidermal',
        'epidermal progenitor',
    ]),
    ('endothelial', [
        'endo progenitor',
        'diffuse endo',
        'endo tail',
        'digestive progenitor',
        'digestive',
    ]),
    ('mesenchymal', [
        #'muscle progenitors',
        'striated muscle',
        #'parenchymal progenitor',
        #'parenchymal',
        'secretory',
    ]),
    ('other', [
        'neuron',
        'pluripotent',
        'germline',
    ]),
]


if __name__ == '__main__':

    # It's juvenile, but the cell types are the same according to Bo
    adata = anndata.read(
       atlas_data_folder / '20230628_hm_adult_annotated.h5ad',
    )

    # FIXME: do something about the awkward double gene IDs/names

    # Take raw matrix
    adata = adata.raw.to_adata()

    # TODO: this normalisation looks really weird, it's indrop data that
    # was generated using customised perl scripts so... let's take it
    # easy for now
    
    # cptt throughout
    sc.pp.normalize_total(
        adata,
        target_sum=1e4,
        key_added='coverage',
    )

    # Set tissue
    adata.obs['tissue'] = 'whole'

    tissues = adata.obs['tissue'].value_counts().index
    compressed_atlas = {}
    for it, tissue in enumerate(tissues):
        print(tissue)
        adata_tissue = adata[adata.obs['tissue'] == tissue]

        # Ignore cells with NaN in the cell.type column
        idx = adata_tissue.obs['cell_types'].isin(
                adata_tissue.obs['cell_types'].value_counts().index)
        adata_tissue = adata_tissue[idx]

        # Fix cell type annotations
        adata_tissue.obs['cellType'] = fix_annotations(
            adata_tissue, 'cell_types', species, tissue,
            rename_dict, coarse_cell_types,
            blacklist=celltype_tissue_blacklist,
        )

        # Age
        adata_tissue.obs['age'] = 'adult'

        # Correction might declare some cells as untyped/low quality
        # they have an empty string instead of an actual annotation
        if (adata_tissue.obs['cellType'] == '').sum() > 0:
            idx = adata_tissue.obs['cellType'] != ''
            adata_tissue = adata_tissue[idx]

        celltypes = get_celltype_order(
            adata_tissue.obs['cellType'].value_counts().index,
            celltype_order,
        )

        print('Average')
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

    print('Consolidate gene list across tissues')
    genes = adata.var_names

    print('Get gene annotations')
    #gene_annos = collect_gene_annotations(anno_fn, genes)
    gene_annos = None

    # Remove existing compressed atlas file if present
    if os.path.isfile(fn_out):
        os.remove(fn_out)

    print('Store compressed atlas to file')
    store_compressed_atlas(
            fn_out,
            compressed_atlas,
            tissues,
            gene_annos,
            celltype_order,
    )
