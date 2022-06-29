# Requirements
`pip install scanpy pyarrow`

# Nanostring CosMx Data 
The file `nanostring_cosmx.py` downloads the data and creates a two output files: 
- basic adata object (`nanostring_lung5_rep2.h5ad`) 
- Single gene points table (`nanostring_lung5_rep2_single_genes.parquet`)


The dataset is stored in `./data_lung5_rep2` and structured as:
- tables:
    - counts: `Lung5_Rep2_exprMat_file.csv`
    - obs: `Lung5_Rep2_metadata_file.csv`
    - fov metadata: `Lung5_Rep2_fov_positions_file.csv`
- images: `CellComposite/`
- regions: `CellLabels/`
- points: `Lung5_Rep2_tx_file.csv`

In addition, the data contains:
- cell outlines in `CellOverlay/` (`.jpp`)
- compartment labels in `CellOverlay/` (`.tif`).