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
    load_config,
    root_repo_folder,
    output_folder,
    get_tissue_data_dict,
    subannotate,
    normalise_counts,
    fix_annotations,
    get_celltype_order,
    compress_tissue,
    collect_feature_annotations,
    store_compressed_atlas,
    )


if __name__ == '__main__':

    species = 'i_pulchra'


    config = load_config(species)
    atlas_data_folder = root_repo_folder / 'data' / 'full_atlases' / 'RNA' / species
    fn_out = output_folder / f'{species}.h5'

    # Remove existing compressed atlas file if present
    if os.path.isfile(fn_out):
        os.remove(fn_out)

    for measurement_type in config["measurement_types"]:
        compressed_atlas = {}
        config_mt = config[measurement_type]
        tissues = config_mt["tissues"]
        celltype_order = config_mt["cell_annotations"]["celltype_order"]
        for tissue in tissues:
            print(tissue)

            print("Read full atlas")
            adata_tissue = anndata.read(atlas_data_folder / (species + '.h5ad'))

            print("Normalise")
            adata_tissue = normalise_counts(adata_tissue, config_mt['normalisation'])

            print("Correct cell annotations")
            adata_tissue = fix_annotations(
                adata_tissue, config_mt['cell_annotations']['column'],
                species,
                tissue,
                config_mt['cell_annotations']['rename_dict'],
                config_mt['cell_annotations']['coarse_cell_types'],
                blacklist=config_mt['cell_annotations']['blacklist'],
            )

            print("Compress atlas")
            compressed_atlas[tissue] = compress_tissue(
                adata_tissue, celltype_order,
            )

        print('Feature annotations')
        feature_annos = collect_feature_annotations(
                config_mt['feature_annotation'],
                adata_tissue.var_names,
                measurement_type,
        )

        print('Store compressed atlas')
        store_compressed_atlas(
                fn_out,
                compressed_atlas,
                tissues,
                feature_annos,
                celltype_order,
        )
