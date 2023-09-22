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


species = 'c_elegans'
atlas_data_folder = root_repo_folder / 'data' / 'full_atlases' / 'RNA' / 'c_elegans'
anno_fn = root_repo_folder / 'data' / 'gene_annotations' / 'c_elegans.PRJNA13758.WS287.annotations.gff3.gz'
fn_out = output_folder / f'{species}.h5'


rename_dict = {
    'cell_types': {
        'Unclassified glia': 'glia',
        'gABAergic neurons': 'GABAergic neuron',
        'am/ph sheath cells': 'sheath',
        'pharyngeal epithelia': 'pharyngeal epi',
        'distal tip cells': 'distal tip',
        'socket cells': 'socket',
        'excretory cells': 'excretory',
        'somatic gonad precursors': 'somatic gonad',
        'unclassified neurons': 'neuron',
        'dopaminergic neurons': 'dopaminergic neuron',
        'cholinergic neurons': 'cholinergic neuron',
        'ciliated sensory neurons': 'ciliated sensory neuron',
        'canal associated neurons': 'canal associated neuron',
        'touch receptor neurons': 'touch receptor neuron',
        'pharyngeal neurons': 'pharyngeal neuron',
        'oxygen sensory neurons': 'oxygen sensory neuron',
        'flp-1(+) interneurons': 'flp-1(+) interneuron',
        'other interneurons': 'neuron',
        'vulval precursors': 'vulval precursor',
        'coelomocytes': 'coelomocyte',
        'seam cells': 'seam',
        'sex myoblasts': 'sex myoblast',
        'gabaergic neurons': 'GABAergic neuron',
        'unclassified glia': 'glia',
    },
}

celltype_tissue_blacklist = {
    'whole': ['Failed QC', 'failed qc'],
}

coarse_cell_types = [
]


celltype_order = [
    ('immune', [
        'glia',
    ]),
    ('epithelial', [
        'seam',
        'non-seam hypodermis',
        'pharyngeal epi',
        'coelomocyte',
        'distal tip',
    ]),
    ('endothelial', [
    ]),
    ('mesenchymal', [
        'body wall muscle',
        'pharyngeal muscle',
        'intestinal/rectal muscle',
        'sex myoblast',
        'sheath',
        'socket',
        'pharyngeal gland',
        'excretory',
        'rectum',
    ]),
    ('other', [
        'dopaminergic neuron',
        'cholinergic neuron',
        'ciliated sensory neuron',
        'GABAergic neuron',
        'canal associated neuron',
        'touch receptor neuron',
        'pharyngeal neuron',
        'oxygen sensory neuron',
        'flp-1(+) interneuron',
        'neuron',
        'germline',
        'somatic gonad',
        'vulval precursor',
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
        idx = adata_tissue.obs['cell.type'].isin(
                adata_tissue.obs['cell.type'].value_counts().index)
        adata_tissue = adata_tissue[idx]

        # cptt throughout
        sc.pp.normalize_total(
            adata_tissue,
            target_sum=1e4,
            key_added='coverage',
        )

        # Fix cell type annotations
        adata_tissue.obs['cell.type'] = adata_tissue.obs['cell.type'].str.lower()
        adata_tissue.obs['cellType'] = fix_annotations(
            adata_tissue, 'cell.type', 'c_elegans', tissue,
            rename_dict, coarse_cell_types,
            blacklist=celltype_tissue_blacklist,
        )

        # Age
        adata_tissue.obs['age'] = 'L2'

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

    print('Get gene annotations')
    gene_annos = collect_gene_annotations(anno_fn, genes)

    print('Store compressed atlas to file')
    store_compressed_atlas(
            fn_out,
            compressed_atlas,
            tissues,
            gene_annos,
            celltype_order,
    )
