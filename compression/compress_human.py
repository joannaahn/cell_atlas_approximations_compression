# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/05/22
content:    Compress Tabula Sapiens.
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


species = 'h_sapiens'
atlas_data_folderd = {
    'RNA': root_repo_folder / 'data' / 'full_atlases' / 'RNA' / 'tabula_sapiens',
    'ATAC': root_repo_folder / 'data' / 'full_atlases' / 'ATAC' / 'h_sapiens',
}
atlas_fn_atac = atlas_data_folderd['ATAC'] / 'Ren_lab_cell_by_cCRE_matrix.h5ad'
atlas_fn_atac_meta = atlas_data_folderd['ATAC'] / 'Ren_lab_Cell_metadata.tsv.gz'
anno_fn = root_repo_folder / 'data' / 'gene_annotations' / 'Homo_sapiens.GRCh38.109.gtf.gz'
fn_out = output_folder / 'tabula_sapiens.h5'


rename_dict = {
    'tissues': {
        # RNA
        'Large_Intestine': 'Colon',

        # ATAC
        'pancreas': 'Pancreas',
        'esophagus_mucosa': 'Esophagus',
        'esophagus_muscularis': 'Esophagus',
        'artery_aorta': 'Heart',
        'muscle': 'Muscle',
        'thyroid': 'Thyroid',
        'colon_sigmoid': 'Colon',
        'stomach': 'Stomach',
        'heart_lv': 'Heart',
        'LungMap': '',
        'esophagus_ge_junction': 'Esophagus',
        'mammary_tissue': 'Mammary',
        'colon_transverse': 'Colon',
        'Human_brain': 'Brain',
        'adrenal_gland': 'Adrenal',
        'small_intestine': 'Colon',
        'ovary': 'Ovary',
        'uterus': 'Uterus',
        'heart_atrial_appendage': 'Heart',
        'liver': 'Liver',
        'nerve_tibial': 'Nerve',
        'artery_tibial': '',
        'lung': 'Lung',
        'skin': 'Skin',
        'vagina': 'Vagina',
        'islet': 'Pancreas',
        'skin_sun_exposed': 'Skin',
        'adipose_omentum': 'Fat', # White fat
        # NOTE: should we keep these??
        'heart_ra_CARE181213': 'Heart',
        'heart_rv_CARE181213': 'Heart',
        'heart_lv_CARE191122': 'Heart',
        'heart_lv_CARE190307': 'Heart',
        'heart_rv_CARE190307': 'Heart',
        'heart_la_CARE191122': 'Heart',
        'heart_rv_CARE181125': 'Heart',
        'heart_lv_CARE181125': 'Heart',
        'heart_rv_CARE190331': 'Heart',
        'heart_la_CARE190307': 'Heart',
        'heart_lv_CARE190331': 'Heart',
        'heart_la_CARE181125': 'Heart',
        'heart_ra_CARE190307': 'Heart',
    },
    'cell_types': {
        # RNA
        'cd24 neutrophil': 'neutrophil',
        'cd4-positive, alpha-beta t cell': 'T',
        'cd8-positive, alpha-beta t cell': 'T',
        'erythroid progenitor': 'erythroid',
        'nk cell': 'NK',
        'hematopoietic stem cell': 'HSC',
        'nampt neutrophil': 'neutrophil',
        'memory b cell': 'B',
        'naive b cell': 'B',
        'myeloid progenitor': 'myeloid',
        'plasmablast': 'plasma cell',
        'enterocyte of epithelium of large intestine': 'enterocyte',
        'immature enterocyte': 'enterocyte',
        'paneth cell of epithelium of large intestine': 'paneth',
        'mature enterocyte': 'enterocyte',
        'b cell': 'B',
        'large intestine goblet cell': 'goblet',
        'transit amplifying cell of large intestine': 'transit amp',
        'goblet cell': 'goblet',
        'intestinal crypt stem cell': 'crypt',
        'intestinal crypt stem cell of large intestine': 'crypt',
        'intestinal enteroendocrine cell': 'enteroendocrine',
        'gut endothelial cell': 'endothelial',
        'mast cell': 'mast',
        'intestinal tuft cell': 'brush',
        'cardiac muscle cell': 'cardiomyocyte',
        'cardiac endothelial cell': 'coronary',
        'fibroblast of cardiac tissue': 'fibroblast',
        'smooth muscle cell': 'smooth muscle',
        'cd4-positive helper t cell': 'T',
        'kidney epithelial cell': 'epithelial',
        'endothelial cell': 'endothelial',
        'type i pneumocyte': 'AT1 epi',
        'type ii pneumocyte': 'AT2 epi',
        'basal cell': 'basal',
        'classical monocyte': 'monocyte',
        'club cell': 'club',
        'non-classical monocyte': 'monocyte',
        'capillary endothelial cell': 'capillary',
        'respiratory goblet cell': 'goblet',
        'lung ciliated cell': 'ciliated',
        'capillary aerocyte': 'CAP2',
        'vein endothelial cell': 'venous',
        'lung microvascular endothelial cell': 'capillary',
        'adventitial cell': 'fibroblast',
        'dendritic cell': 'dendritic',
        'intermediate monocyte': 'monocyte',
        'pericyte cell': 'pericyte',
        'endothelial cell of artery': 'arterial',
        'cd4-positive alpha-beta t cell': 'T',
        'bronchial smooth muscle cell': 'smooth muscle',
        'vascular associated smooth muscle cell': 'vascular smooth muscle',
        'cd8-positive alpha-beta t cell': 'T',
        'endothelial cell of lymphatic vessel': 'lymphatic',
        'bronchial vessel endothelial cell': 'capillary',
        'pulmonary ionocyte': 'ionocyte',
        'plasmacytoid dendritic cell': 'plasmacytoid',
        'mesothelial cell': 'mesothelial',
        'serous cell of epithelium of bronchus': 'serous',
        'myofibroblast cell': 'smooth muscle',
        'respiratory mucous cell': 'mucous',
        'pancreatic acinar cell': 'acinar',
        'pancreatic ductal cell': 'ductal',
        'myeloid cell': 'myeloid',
        't cell': 'T',
        'pancreatic stellate cell': 'stellate',
        'pancreatic beta cell': 'beta',
        'pancreatic pp cell': 'PP',
        'pancreatic alpha cell': 'alpha',
        'pancreatic delta cell': 'delta',
        'epithelial cell': 'epithelial',
        'tongue muscle cell': 'striated muscle',
        'schwann cell': 'schwann',
        # ATAC
        'Transitional Zone Cortical Cell': 'cortical',
        'Zona Fasciculata Cortical Cell': 'cortical',
        'Zona Glomerulosa Cortical Cell': 'cortical',
        'Cortical Epithelial-like': 'cortical epi-like',
        'Fibroblast (Liver Adrenal)': 'fibroblast',
        'Macrophage (General)': 'macrophage',
        'Endothelial (Exocrine Tissues)': 'capillary',
        'Pericyte (General) 4': 'pericyte',
        'Fibroblast (General)': 'fibroblast',
        'Fibroblast (Peripheral Nerve)': 'fibroblast',
        'T Lymphocyte 1 (CD8+)': 'T',
        'Macrophage (General,Alveolar)': 'macrophage',
        'Endothelial Cell (General) 1': 'capillary',
        'Lymphatic Endothelial Cell': 'lymphatic',
        'Schwann Cell (General)': 'schwann',
        'Luteal Cell (Ovarian)': 'luteal',
        'Adipocyte': 'adipocyte',
        'Naive T cell': 'T',
        'Cardiac Pericyte 4': 'pericyte',
        'Natural Killer T Cell': 'NK',
        'Pericyte (General) 3': 'pericyte',
        'Memory B Cell': 'B',
        'Cardiac Pericyte 3': 'pericyte',
        'T lymphocyte 2 (CD4+)': 'T',
        'Endothelial Cell (General) 2': 'capillary',
        'Pericyte (General) 1': 'pericyte',
        'Mast Cell': 'mast',
        'Smooth Muscle (Vaginal)': 'smooth muscle',
        'Vascular Smooth Muscle 2': 'vascular smooth muscle',
        'Smooth Muscle (General)': 'smooth muscle',
        'CNS,Enteric Neuron': 'neuron',
        'Glutamatergic Neuron 1': 'neuron',
        'Glutamatergic Neuron 2': 'neuron',
        'Oligodendrocyte': 'oligodendrocyte',
        'GABAergic Neuron 1': 'neuron',
        'GABAergic Neuron 2': 'neuron',
        'Oligodendrocyte Precursor': 'OPC',
        'Microglia': 'glial',
        'Astrocyte 1': 'astrocyte',
        'Astrocyte 2': 'astrocyte',
        'Blood Brain Barrier Endothelial Cell': 'capillary',
        'Cardiac Fibroblasts': 'fibroblast',
        'Cardiac Pericyte 2': 'pericyte',
        'Endothelial Cell (Myocardial)': 'capillary',
        'Pericyte (General) 2': 'pericyte',
        'Colon Epithelial Cell 1': 'epithelial',
        'Small Intestinal Enterocyte': 'enterocyte',
        'Smooth Muscle (Colon) 1': 'smooth muscle',
        'Fibroblast (Gastrointestinal)': 'fibroblast',
        'Smooth Muscle (Colon) 2': 'smooth muscle',
        'Colonic Goblet Cell': 'goblet',
        'Plasma Cell': 'plasma cell',
        'Small Intestinal Goblet Cell': 'goblet',
        'Colon Epithelial Cell 2': 'epithelial',
        'Colon Epithelial Cell 3': 'epithelial',
        'Enterochromaffin Cell': 'enteroendocrine',
        'Smooth Muscle (General Gastrointestinal)': 'smooth muscle',
        'Tuft Cell': 'tuft',
        'Paneth Cell': 'paneth',
        'Smooth Muscle (Esophageal Muscularis) 3': 'smooth muscle',
        'Endothelial Cell (General) 3': 'capillary',
        'Pericyte (Esophageal Muscularis)': 'pericyte',
        'Smooth Muscle (Esophageal Mucosal)': 'smooth muscle',
        'Smooth Muscle (GE Junction)': 'smooth muscle',
        'Smooth Muscle (Esophageal Muscularis) 1': 'smooth muscle',
        'Smooth Muscle (Uterine)': 'smooth muscle',
        'Vascular Smooth Muscle 1': 'vascular smooth muscle',
        'Smooth Muscle (Esophageal Muscularis) 2': 'smooth muscle',
        'Fibroblast (Epithelial)': 'fibroblast',
        'Fibroblast (Sk Muscle Associated)': 'fibroblast',
        'Gastric Neuroendocrine Cell': 'neuroendocrine',
        'Alveolar Capillary Endothelial Cell': 'capillary',
        'Keratinocyte 1': 'keratinocyte',
        'Mesothelial Cell': 'mesothelial',
        'Melanocyte': 'melanocyte',
        'Esophageal Epithelial Cell': 'epithelial',
        'Foveolar Cell': 'foveolar',
        'Type I Skeletal Myocyte': 'striated muscle',
        'Myoepithelial (Skin)': 'myoepithelial',
        'Satellite Cell': 'satellite',
        'Granular Epidermal (Skin)': 'epidermal',
        'Basal Epidermal (Skin)': 'basal',
        'Club Cell': 'club',
        'Type II Skeletal Myocyte': 'striated muscle',
        'Keratinocyte 2': 'keratinocyte',
        'Chief Cell': 'chief',
        'Ventricular Cardiomyocyte': 'cardiomyocyte',
        'Atrial Cardiomyocyte': 'cardiomyocyte',
        'Cardiac Pericyte 1': 'pericyte',
        'Endocardial Cell': 'endocardial',
        'Hepatocyte': 'hepatocyte',
    },
}


celltype_tissue_blacklist = {
    'Adrenal': ['Luteal Cell (Ovarian)'],
    'Brain': ['Peripheral Nerve Stromal'],
    'Colon': [
        'Mammary Luminal Epithelial Cell 2', 'Pancreatic Delta,Gamma cell',
        'Ductal Cell (Pancreatic)',
        'Chief Cell', # stomach
    ],
    'Esophagus': [
        'Airway Goblet Cell',
        'Peripheral Nerve Stromal',
        'Mammary Luminal Epithelial Cell 2',
        'Basal Epithelial (Mammary)',
        'Thyroid Follicular Cell',
        'Pancreatic Acinar Cell',
    ],
    'Fat': [
        'Alverolar Type 2,Immune',
        'Peripheral Nerve Stromal',
        'Pancreatic Acinar Cell',
        'Ductal Cell (Pancreatic)',
    ],
    'Heart': [
        'Peripheral Nerve Stromal',
        'Ductal Cell (Pancreatic)',
        'Mammary Luminal Epithelial Cell 2',
        'Pancreatic Acinar Cell',
        'Alveolar Type 2 (AT2) Cell',
        'Thyroid Follicular Cell',
        'Alveolar Type 1 (AT1) Cell',
    ],
    'Liver': [
        'Ductal Cell (Pancreatic)',
        'Mammary Luminal Epithelial Cell 2',
    ],
}


coarse_cell_types = [
    'endothelial',
    'immune cell',
]



celltype_order = [
    ('immune', [
        'HSC',
        'neutrophil',
        'basophil',
        'granulocyte',
        'mast',
        'myeloid',
        'monocyte',
        'macrophage',
        'dendritic',
        'erythroid',
        'erythrocyte',
        'B',
        'plasma cell',
        'T',
        'NK',
        'plasmacytoid',
        'glial',
    ]),
    ('epithelial', [
        'epithelial',
        'goblet',
        'brush',
        'crypt',
        'transit amp',
        'enterocyte',
        'paneth',
        'AT1 epi',
        'AT2 epi',
        'club',
        'ciliated',
        'ductal',
        'acinar',
        'keratinocyte',
        'basal',
        'serous',
        'mucous',
        'cortical epi-like',
        'tuft',
        'melanocyte',
        'foveolar',
        'myoepithelial',
        'chief',
        'epidermal',
    ]),
    ('endothelial', [
        'arterial',
        'venous',
        'coronary',
        'capillary',
        'CAP2',
        'lymphatic',
        'endocardial',
    ]),
    ('mesenchymal', [
        'fibroblast',
        'alveolar fibroblast',
        'cardiomyocyte',
        'stellate',
        'striated muscle',
        'smooth muscle',
        'vascular smooth muscle',
        'pericyte',
        'mesothelial',
        'satellite',
    ]),
    ('other', [
        'enteroendocrine',
        'neuroendocrine',
        'hepatocyte',
        'ionocyte',
        'alpha',
        'beta',
        'PP',
        'delta',
        'schwann',
        'adipocyte',
        'cortical',
        'luteal',
        'neuron',
        'oligodendrocyte',
        'OPC',
        'astrocyte',
    ]),
    ('unknown', [
        'unknown',
    ])
]


if __name__ == '__main__':

    # Remove existing compressed atlas file if present
    if os.path.isfile(fn_out):
        os.remove(fn_out)

    if False:
        print('RNA')
        compressed_atlas = {}
        tissue_sources = get_tissue_data_dict(
                'human', atlas_data_folder, rename_dict)
        tissues = list(tissue_sources.keys())
        for it, (tissue, full_atlas_fn) in enumerate(tissue_sources.items()):
            print(tissue)

            adata_tissue = anndata.read(full_atlas_fn)

            # Restart from raw data and renormalize
            adata_tissue = adata_tissue.raw.to_adata()

            # cptt throughout
            sc.pp.normalize_total(
                adata_tissue,
                target_sum=1e4,
                key_added='coverage',
            )

            # Fix cell type annotations
            adata_tissue.obs['cellType'] = fix_annotations(
                adata_tissue, 'cell_ontology_class', 'human', tissue,
                rename_dict, coarse_cell_types,
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

            print('Add data to celltype-timepoint group')
            # NOTE: see supplementary table 1 of the Science paper
            adata_tissue.obs['age'] = adata_tissue.obs['Life stage']
            ages = adata_tissue.obs['age'].value_counts().index.tolist()
            ages.sort()
            columns_age = []
            for ct in celltypes:
                for age in ages:
                    columns_age.append('_'.join([ct, 'TS', str(age)]))

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
                    label = '_'.join([celltype, 'TS', str(age)])
                    avg_ge_tp[label] = np.asarray(Xct_age.mean(axis=0))[0]
                    frac_ge_tp[label] = np.asarray((Xct_age > 0).mean(axis=0))[0]
                    ncells_ge_tp[label] = len(idx_age)

            compressed_atlas[tissue] = {
                'features': genes,
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

        print('Add gene annotations')
        gene_annos = collect_gene_annotations(anno_fn, genes)

        print('Store compressed atlas to file')
        store_compressed_atlas(
                fn_out,
                compressed_atlas,
                tissues,
                gene_annos,
                celltype_order,
                measurement_type='gene_expression',
        )

    if True:
        print('ATAC')
        compressed_atlas = {}
        adata_all = anndata.read(atlas_fn_atac, backed='r')

        print('Load metadata')
        meta = pd.read_csv(atlas_fn_atac_meta, sep='\t', index_col=0).loc[adata_all.obs_names]

        print('Ignore fetal data for now')
        meta = meta.loc[meta['Life stage'] == 'Adult']

        print('Rename tissues and exclude tissues with no name')
        meta['tissue.orig'] = meta['tissue']
        meta['tissue'] = pd.Series(
                [x[:x.rfind('_')] for x in meta['tissue'].values],
                index=meta.index,
                ).map(rename_dict['tissues'])
        meta = meta.loc[meta['tissue'] != '']

        print('Get list of (renames) tissues')
        tissues = np.sort(meta['tissue'].unique())

        print('Iterate over tissues')
        # FIXME
        for it, tissue in enumerate(tissues[7:]):
            print(tissue)

            meta_tissue = meta.loc[meta['tissue'] == tissue]

            # Load into memory
            adata_tissue = adata_all[meta_tissue.index].to_memory()
            adata_tissue.obs = meta_tissue.copy()

            # Binarise ATAC counts
            adata_tissue.X.data[:] = 1

            # Fix cell type annotations
            adata_tissue.obs['cellType'] = fix_annotations(
                adata_tissue, 'cell type', 'human', tissue,
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
            features = adata_tissue.var_names
            avg_ca = pd.DataFrame(
                    np.zeros((len(features), len(celltypes)), np.float32),
                    index=features,
                    columns=celltypes,
                    )
            ncells_ca = pd.Series(
                    np.zeros(len(celltypes), np.int64), index=celltypes,
                    )
            for celltype in celltypes:
                idx = adata_tissue.obs['cellType'] == celltype
                Xidx = adata_tissue[idx].X
                avg_ca[celltype] = np.asarray(Xidx.mean(axis=0))[0]
                ncells_ca[celltype] = idx.sum()

            print('Add data to celltype-timepoint group')
            # NOTE: see supplementary table 1 of the Science paper
            adata_tissue.obs['age'] = adata_tissue.obs['Life stage'].copy()
            ages = adata_tissue.obs['age'].value_counts().index.tolist()
            ages.sort()
            columns_age = []
            for ct in celltypes:
                for age in ages:
                    columns_age.append('_'.join([ct, 'BR:ATAC', str(age)]))

            # Averages
            features = adata_tissue.var_names
            avg_ca_tp = pd.DataFrame(
                    np.zeros((len(features), len(celltypes) * len(ages)), np.float32),
                    index=features,
                    columns=columns_age,
                    )
            ncells_ca_tp = pd.Series(
                    np.zeros(len(columns_age), np.int64), index=columns_age,
                    )
            for celltype in celltypes:
                adata_ct = adata_tissue[adata_tissue.obs['cellType'] == celltype]
                for age in ages:
                    idx_age = (adata_ct.obs['age'] == age).values.nonzero()[0]
                    if len(idx_age) == 0:
                        continue
                    Xct_age = adata_ct.X[idx_age]
                    label = '_'.join([celltype, 'BR:ATAC', str(age)])
                    avg_ca_tp[label] = np.asarray(Xct_age.mean(axis=0))[0]
                    ncells_ca_tp[label] = len(idx_age)

            compressed_atlas[tissue] = {
                'features': features,
                'celltype': {
                    'avg': avg_ca,
                    'ncells': ncells_ca,
                },
                'celltype_dataset_timepoint': {
                    'avg': avg_ca_tp,
                    'ncells': ncells_ca_tp,
                },
            }

        #print('Add peak annotations')
        #feature_annos = collect_gene_annotations(anno_fn, features)

        print('Store compressed atlas to file')
        store_compressed_atlas(
                fn_out,
                compressed_atlas,
                tissues,
                None,
                celltype_order,
                measurement_type='chromatin_accessibility',
        )

