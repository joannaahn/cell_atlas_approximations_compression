# Cell Atlas Approximations - compression
Cell atlases such as Tabula Muris and Tabula Sapiens are multi-organ single cell omics data sets describing entire organisms. A cell atlas approximation is a lossy and lightweight compression of a cell atlas that can be streamed via the internet.

This project enables biologists, doctors, and data scientist to quickly find answers for questions such as:

- *What types of cells populate the human heart?*
- *What is the expression of a specific gene across cell types in C elegans?*
- *What are the marker genes of a specific cell type in mouse pancreas*?
- *What fraction of cells (of a specific type) express a gene of interest?*

To answer these questions, we construct approximations from all kinds of cell atlases using the code in this repo.

**NOTE:** At this stage, this repo is designed for internal use and made public for transparency. End users are expected to interact with the atlas approximations using the [API](https://atlasapprox.readthedocs.io) or - when we are done building it - the web interface.


## Proposing a new atlas approximation
If you would like to propose an additional cell atlas to be added to the list of approximations, please open an [issue](https://github.com/fabilab/cell_atlas_approximations_compression/issues/new/choose) and specify:

- a hyperlink to a publication where the atlas is described
- a hyperlink to a data source where the table of counts can be downloaded
- any information you might have about cell type annotation
- the reason why you would like to have that atlas added here

Thank you for your time.
