import os
import pathlib
import gzip
import numpy as np
import pandas as pd
import h5py
import re
import scanpy as sc

# root_repo_folder = pathlib.Path(__file__).parent.parent
# output_folder = root_repo_folder / 'data' / 'atlas_approximations'

input_directory = "../data/full_atlases/drosophila_melanogaster"
output_folder = "./compression_output"

def get_tissue_data_dict(species, rename_dict=None):
    result = []

    # for x in input_folder:
    #     if x has '.h5ad'
    #         fns.append(x)

    filenames = os.listdir(input_directory)
    fns = [x for x in filenames if '.h5ad' in x]
    print(filenames)

    # use regex to search a match string
    # my_name = r'^joanna$'

    pattern = r'(?:s|r)_fca_biohub_(.*?)_10x\.h5ad'
    for filename in fns:
        print(filename)
        tissue = re.search(pattern, filename).group(1)
        result.append({
            'tissue': tissue,
            # 'filename_count': filename,
            # 'filename_meta': filename.replace('Matrix_', 'Metadata_')[:-3],
        })
    print(result)
    # assigning new value to result
    # result = pd.DataFrame(result).set_index('tissue')
    # Order tissues alphabetically
    # result = result.sort_index()
    return result

output = get_tissue_data_dict('d.melanogaster')
print(output)