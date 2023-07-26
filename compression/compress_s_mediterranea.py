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
    )


species = 's_mediterranea'
dataset_name = 'Plass'
atlas_data_folder = root_repo_folder / 'data' / 'full_atlases' / 'RNA' / 's_mediterranea'
#atlas_matrix_fn = atlas_data_folder / 'dge.txt.gz'
#atlas_obs_fn = atlas_data_folder / 'Planaria_Seurat_annot.csv'
# Annotations are already in the h5ad
#anno_fn = root_repo_folder / 'data' / 'gene_annotations' / ...
fn_out = output_folder / f'{species}.h5'

# TODO
rename_dict = {
    'tissues': {},
    'cell_types': {
        # see also https://www.science.org/doi/10.1126/science.aaq1723#sec-10
        'neoblast 1': 'pluripotent',
        'neoblast 4': 'pluripotent',
        'neoblast 6': 'pluripotent',
        'neoblast 9': 'pluripotent',
        'neoblast 5': 'pluripotent',
        'neoblast 8': 'neural progenitor',
        'neoblast 13': 'pluripotent',
        'neoblast 12': 'neural progenitor',
        'neural progenitors': 'neural progenitor',
        'pharynx cell type': 'pharyngeal',  # some kind of epi it seems
        'pharynx cell type progenitors': 'pharyngeal progenitor',  # some kind of epi it seems
        'epidermis': 'epidermal',
        'late epidermal progenitors 1': 'epidermal',
        'late epidermal progenitors 2': 'epidermal',
        'muscle body': 'striated muscle',
        'muscle pharynx': 'striated muscle',
        'ChAT neurons 1': 'neuron',
        'ChAT neurons 2': 'photoreceptor',
        'GABA neurons': 'neuron',
        'cav-1+ neurons': 'neuron',
        'npp-18+ neurons': 'neuron',
        'spp-11+ neurons': 'neuron',
        'phagocytes': 'macrophage',
        'glia': 'glial',
        'protonephridia': 'flame',
        'goblet cells': 'goblet',
        'early epidermal progenitors': 'epidermal progenitor',
        'late epidermal progenitors': 'epidermal progenitor',
        'epidermal neoblasts': 'epidermal progenitor',
        'parenchymal progenitors': 'parenchymal progenitor',
        'aqp+ parenchymal cells': 'aqp+ parenchymal',
        'otf+ cells 1': 'otf+ 1',
        'otf+ cells 2': 'otf+ 2',
        'ldlrr-1+ parenchymal cells': 'parenchymal progenitor',
        'pgrn+ parenchymal cells': 'parenchymal progenitor',
        'psap+ parenchymal cells': 'psap+ parenchymal',
        'gut progenitors': 'macrophage progenitor',
        'psd+ cells': 'psd+ support',  # human: neurons, bladder muscle, or endometrium fibros (?)
        'secretory 1': 'secretory',
        'secretory 2': 'secretory',
        'secretory 3': 'secretory',
        'secretory 4': 'secretory',
    },
}

celltype_tissue_blacklist = {
    'whole': [
        #f'unk_{x}' for x in range(1, 22),
        # A bunch of neoblast-like cells that are rare and poorly described
        'neoblast 2',
        'neoblast 3',
        'neoblast 7',
        'neoblast 10',  # what is this??
        'neoblast 11',  # not found in the wild-type only t-SNE (S20)
        'activated early epidermal progenitors', # from the PAGA this does not look clear
        # Ignore DV boundary for now, unclear if we can trust the classification
        'epidermis DVb',
        'epidermis DVb neoblast',
    ],
}

coarse_cell_types = [
]


celltype_order = [
    ('immune', [
        'macrophage',
        'macrophage progenitor',
        'glial',
    ]),
    ('epithelial', [
        #'comb',
        'epidermal',
        'epidermal progenitor',
        'pharyngeal',
        'pharyngeal progenitor',
        'flame',
        'goblet',
    ]),
    ('endothelial', [
    ]),
    ('mesenchymal', [
        #'digestive',
        'muscle progenitors',
        'striated muscle',
        #'smooth muscle',
        'pigment',
        'parenchymal progenitor',
        'aqp+ parenchymal',
        'psap+ parenchymal',
    ]),
    ('other', [
        #'venom',
        #'lens',
        'neuron',
        'neural progenitor',
        'photoreceptor',
        'pluripotent',
        'otf+ 1',
        'otf+ 2',
        'psd+ support',
        'secretory',
    ]),
]


if __name__ == '__main__':

    # Remove existing compressed atlas file if present
    if os.path.isfile(fn_out):
        os.remove(fn_out)

    compressed_atlas = {}

    adata = anndata.read(
       atlas_data_folder / 's_mediterranea_with_cell_and_transcript_annotations.h5ad',
    )
    # Rename for the sake of mental sanity
    adata.obs.rename(columns={'final_Id': 'cell_type'}, inplace=True)

    # Only keep fresh samples, not the regenerating ones
    # Most cell types are well represented anyway, the other ones appear questionable
    # in the first place
    adata = adata[adata.obs['run'].str.startswith('plan')]

    # Set tissue
    adata.obs['tissue'] = 'whole'

    tissues = adata.obs['tissue'].value_counts().index
    for it, tissue in enumerate(tissues):
        print(tissue)
        adata_tissue = adata[adata.obs['tissue'] == tissue]

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
        adata_tissue.obs['cellType'] = fix_annotations(
            adata_tissue, 'cell_type', species, tissue,
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

    print('Consolidate gene list across tissues')
    genes = adata.var_names

    print('Get gene annotations')
    # This was a TOTAL pain but could be done via
    # https://planmine.mpinat.mpg.de/planmine/bagDetails.do?bagName=All%20contigs%20in%20dd_Smed
    gene_annos = adata.var.rename(
        columns={
            'Transcript > Id': 'trascript_id',
            'Transcript > Chromosome Location > Strand': 'strand',
            'Transcript > Chromosome Location > Start': 'start_position',
            'Transcript > Chromosome Location > End': 'end_position',
        })
    gene_annos['chromosome_name'] = (
        gene_annos.index
                  .str
                  .split('_', expand=True)
                  .get_level_values(3))
    # FIXME: is this vaguely correct??
    gene_annos['gene_name'] = gene_annos.index

    print('Store compressed atlas to file')
    store_compressed_atlas(
            fn_out,
            compressed_atlas,
            tissues,
            gene_annos,
            celltype_order,
    )
