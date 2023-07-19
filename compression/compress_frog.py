# vim: fdm=indent
"""
author:     Fabio Zanini
date:       06/12/22
content:    Compress frog data.
"""
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


species = "x_laevis"
atlas_data_folder = root_repo_folder / "data" / "full_atlases" / "RNA" / species
atlas_data_fn = atlas_data_folder / "Xenopus_Figure1.h5ad"
anno_fn = (
    root_repo_folder / "data" / "gene_annotations" / "GCF_017654675.1_Xenopus_laevis_v10.1_genomic.gtf.gz"
)
fn_out = output_folder / f"{species}.h5"


rename_dict = {
    "tissues": {
        "Bone Marrow": "Bone_Marrow",
        "Intestine": "Gut",
    },
    "cell_types": {
        "Eosinophil": "eosinophil",
        "Macrophage": "macrophage",
        "Erythrocyte": "erythrocyte",
        "B cell": "B",
        "T cell": "T",
        "Stromal cell": "stromal",
        "Smooth muscle cell": "smooth muscle",
        "Mast cell": "mast",
        "Chondrocyte": "chondrocyte",
        "Fibroblast": "fibroblast",
        #'Radial glia': 'glia',
        "Acinar cell": "acinar",
        "Beta cell": "beta",
        "Delta cell": "delta",
        ("Pancreas", "Epithelial cell"): "ductal",
        "Neuroendocrine cell": "neuroendocrine",
        "Neuron": "neuron",
        "Oligodendrocyte": "oligodendrocyte",
        "Astrocyte": "astrocyte",
        "Radial glia": "glia",
        "GABAergic neuron": "neuron",
        "Endothelial cell": "capillary",
        "Rod photoreceptor": "rod",
        "Melanocyte": "melanocyte",
        "AT2 cell": "AT2 epi",
        "Cumulus cell": "cumulus",
        "Cardiomyocyte": "cardiomyocyte",
        "Stomach Pit cell": "pit",
        "Stomach Parietal cell": "parietal",
        "Secretory cell": "chief",
        ("Stomach", "Epithelial cell"): "mucous",
        "Spermatogonia": "spermatogonia",
        "Spermatocyte": "spermatocyte",
        "Spermatid": "spermatid",
        "Liver sinusoidal endothelial cell": "capillary",
        "Hepatocyte": "hepatocyte",
        "Secretory cell": "secretory",
        "Proximal tubule cell": "proximal tubule epi",
        "Collecting duct intercalated cell": "collecting duct epi",
        "Enterocyte": "enterocyte",
        "Goblet cell": "goblet",
        ("Gut", "Epithelial cell"): "epithelial",
        "Lens fiber cell": "corneal",
        ("Ovary", "Cumulus cell"): "cumulus",
        "Epidermal stem cell": "epidermal",
        ("Skin", "Epithelial cell"): "epidermal",
        ("Bladder", "Epithelial cell"): "epithelial",
    },
}

# FIXME: Myeloid cell and leukocyte are coarse cell types, move them there
celltype_tissue_blacklist = {
    "Bone_Marrow": [
        "Radial glia",
    ],
    "Pancreas": [
        "Endothelial cell",
        "Enterocyte",
        "Radial glia",
        "Liver sinusoidal endothelial cell",
    ],
    "Spleen": [
        "AT2 cell",
    ],
    "Brain": [
        "Epithelial cell",
    ],
    "Lung": [
        "Collecting duct intercalated cell",
        "Epithelial cell",
        "Enterocyte",
        "Cardiomyocyte",
    ],
    "Heart": [
        "Liver sinusoidal endothelial cell",
    ],
    "Stomach": [
        "Collecting duct intercalated cell",
    ],
    "Testis": [
        "Testis unknown cell",
        "Cumulus cell",  # support the ovary
    ],
    "Liver": [],
    "Muscle": [
        "Cardiomyocyte",
        "Oligodendrocyte",
    ],
    "Kidney": [
        "Liver sinusoidal endothelial cell",
        "Spermatocyte",
        "Spermatogonia",
        "Spermatid",
        "Epithelial cell",
        "Epidermal stem cell",
        "Secretory cell",  # this should be a coarse cell type: it's unclear which type they meant
    ],
    "Gut": [
        "Stomach Pit cell",
        "Liver sinusoidal endothelial cell",
    ],
    "Eye": [
        "Cardiomyocyte",
        "Collecting duct intercalated cell",
        "Oligodendrocyte",
        "Cumulus cell",
        "Epidermal stem cell",
        "Epithelial cell",
    ],
    "Ovary": [
        "Spermatogonia",
    ],
    "Bladder": [
        "Collecting duct intercalated cell",
        "Spermatid",
        "Spermatocyte",
        "Goblet cell",
        "Spermatogonia",
        "Hepatocyte",
    ]
}

coarse_cell_types = [
    "Myeloid cell",
    "Leukocyte",
]

subannotation_kwargs = {
    "markers": {},
    "bad_prefixes": [],
    "skip_subannotation": True,
}

# We store an organism-wide complete ordering of cell types, and each tissue
# will cherry pick the necessary
celltype_order = [
    (
        "immune",
        [
            #'HSC',
            #'hematopoietic',
            "neutrophil",
            #'basophil',
            "mast",
            "eosinophil",
            #'granulocytopoietic',
            #'granulocyte',
            #'promonocyte',
            #'myeloid',
            #'monocyte',
            #'alveolar macrophage',
            "macrophage",
            #'dendritic',
            #'Langerhans',
            #'megakaryocyte-erythroid',
            #'proerythroblast',
            #'erythroblast',
            #'erythroid',
            "erythrocyte",
            #'precursor B',
            #'late pro-B',
            #'immature B',
            "B",
            #'plasma cell',
            "T",
            #'Treg',
            #'NKT',
            #'NK',
            #'plasmacytoid',
            "glia",
        ],
    ),
    (
        "epithelial",
        [
            'epithelial',
            'goblet',
            #'brush',
            #'crypt',
            'enterocyte',
            "proximal tubule epi",
            #'distal tubule epi',
            #'podocyte',
            #'Henle limb epi',
            "collecting duct epi",
            "AT2 epi",
            #'club',
            #'ciliated',
            "ductal",
            "acinar",
            #'keratinocyte',
            #'basal',
            "melanocyte",
            "pit",
            "parietal",
            "chief",
            "mucous",
            "secretory",
            "corneal",
            "epidermal",
        ],
    ),
    (
        "endothelial",
        [
            #'arterial',
            #'venous',
            #'coronary',
            #'fenestrated',
            "capillary",
            #'lymphatic',
        ],
    ),
    (
        "mesenchymal",
        [
            "fibroblast",
            #'alveolar fibroblast',
            #'endocardial',
            #'ventricular',
            #'stellate',
            "chondrocyte",
            "cardiomyocyte",
            "smooth muscle",
            #'vascular smooth muscle',
            #'pericyte',
            #'mesangial',
            "stromal",  # bone marrow
        ],
    ),
    (
        "other",
        [
            "neuron",
            "astrocyte",
            "oligodendrocyte",
            "rod",
            "neuroendocrine",
            #'alpha',
            "beta",
            #'PP',
            "delta",
            "hepatocyte",
            "spermatogonia",
            "spermatocyte",
            "spermatid",
            "cumulus",
        ],
    ),
]


if __name__ == "__main__":
    # Remove existing compressed atlas file if present
    if os.path.isfile(fn_out):
        os.remove(fn_out)

    compressed_atlas = {}

    print(f"Load whole animal data ({species})")
    adata = anndata.read(atlas_data_fn)

    print("Exclude developmental stages")
    adata = adata[~adata.obs["tissue"].isin(["St48", "St66", "St54", "St59"])]

    print("Rename tissues")
    adata.obs["tissue"] = adata.obs["tissue"].map(
        lambda x: rename_dict["tissues"].get(x, x),
    )

    print("Renormalize to cptt throughout")
    adata = adata.raw.to_adata()
    adata.X.data = np.exp(adata.X.data) - 1

    tissues = list(adata.obs["tissue"].value_counts().index)
    for it, tissue in enumerate(tissues):
        print(tissue)
        adata_tissue = adata[adata.obs["tissue"] == tissue]

        # Ignore cell types with too few cells... probably misannotations
        ct_counts = adata_tissue.obs["anno"].value_counts()
        cell_types_abundant = list(ct_counts.index[ct_counts >= 5])
        adata_tissue = adata_tissue[
            adata_tissue.obs["anno"].isin(
                cell_types_abundant,
            )
        ]

        # Fix cell type annotations
        adata_tissue.obs["cellType"] = fix_annotations(
            adata_tissue,
            "anno",
            "frog",
            tissue,
            rename_dict,
            coarse_cell_types,
            blacklist=celltype_tissue_blacklist,
            subannotation_kwargs=subannotation_kwargs,
        )

        # Correction might declare some cells as untyped/low quality
        # they have an empty string instead of an actual annotation
        if (adata_tissue.obs["cellType"].isin(["", "unknown"])).sum() > 0:
            idx = ~adata_tissue.obs["cellType"].isin(["", "unknown"])
            adata_tissue = adata_tissue[idx]

        celltypes = get_celltype_order(
            adata_tissue.obs["cellType"].value_counts().index,
            celltype_order,
        )

        print(f"Add data to celltype group: {tissue}")
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
            np.zeros(len(celltypes), np.int64),
            index=celltypes,
        )
        for celltype in celltypes:
            idx = adata_tissue.obs["cellType"] == celltype
            Xidx = adata_tissue[idx].X
            avg_ge[celltype] = np.asarray(Xidx.mean(axis=0))[0]
            frac_ge[celltype] = np.asarray((Xidx > 0).mean(axis=0))[0]
            ncells_ge[celltype] = idx.sum()

        compressed_atlas[tissue] = {
            "features": genes,
            "celltype": {
                "avg": avg_ge,
                "frac": frac_ge,
                "ncells": ncells_ge,
            },
        }

    print("Get gene list")
    genes = adata.var_names.copy()

    print("Add gene annotations")
    gene_annos = collect_gene_annotations(anno_fn, genes)

    print("Store compressed atlas to file")
    store_compressed_atlas(
        fn_out,
        compressed_atlas,
        tissues,
        gene_annos,
        celltype_order,
    )
