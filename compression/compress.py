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
    filter_cells,
    )


if __name__ == '__main__':

    species_list = [
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

        # Multi-organ species
        'h_sapiens',
        'm_musculus',
        'm_myoxinus',
        'd_melanogaster',
        'x_laevis',
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

        for measurement_type in config["measurement_types"]:
            compressed_atlas = {}
            config_mt = config[measurement_type]
            tissues = config_mt["tissues"]
            celltype_order = config_mt["cell_annotations"]["celltype_order"]

            if "path_global" in config_mt:
                print(f"Read full atlas")
                adata = anndata.read(config_mt["path_global"])

            if "path_metadata_global" in config_mt:
                print("Read global metadata separately")
                meta = pd.read_csv(config_mt["path_metadata_global"], sep='\t', index_col=0).loc[adata.obs_names]

                if ("filter_cells_global" in config_mt) and ("metadata" in config_mt["filter_cells_global"]):
                    column, value = config_mt["filter_cells_global"]["metadata"]
                    meta = meta.loc[meta[column] == value]


            for tissue in tissues:
                print(tissue)

                if "path_metadata_global" in config_mt:
                    meta_tissue = meta.loc[meta['tissue'] == tissue]

                if "path_global" not in config_mt:
                    print(f"Read full atlas for {tissue}")
                    adata_tissue = anndata.read(config_mt["path"][tissue])
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
