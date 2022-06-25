import os
import numpy as np
import pandas as pd
import scanpy as sc

from anndata import AnnData
from pathlib import Path
from scipy.sparse import csr_matrix

nanostring_dir = Path().resolve() / "nanostring"
assert nanostring_dir.exists()
path = nanostring_dir / "data_lung5_rep2"

if not path.exists():
    url = "https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/Lung5_Rep2/Lung5_Rep2+SMI+Flat+data.tar.gz"
    ugly_name = "Lung5_Rep2/Lung5_Rep2-Flat_files_and_images/"
    os.system(f"wget -P {nanostring_dir} {url}")
    os.system(
        f"tar -xzf {nanostring_dir/'Lung5_Rep2+SMI+Flat+data.tar.gz'} -C {nanostring_dir}"
    )
    os.system(f"mv {nanostring_dir/ugly_name} {path}")
    os.system(f"rm -r {nanostring_dir/'Lung5_Rep2'}")
    os.system(f"rm {nanostring_dir/'Lung5_Rep2+SMI+Flat+data.tar.gz'}")

counts_file = "Lung5_Rep2_exprMat_file.csv"
meta_file = "Lung5_Rep2_metadata_file.csv"
fov_file = "Lung5_Rep2_fov_positions_file.csv"

# Gene counts per cell
fov_key = "fov"
cell_id_key = "cell_ID"
counts = pd.read_csv(path / counts_file, header=0, index_col=cell_id_key)
counts.index = counts.index.astype(str).str.cat(
    counts.pop(fov_key).astype(str).values, sep="_"
)

# Obs per cell
obs = pd.read_csv(path / meta_file, header=0, index_col=cell_id_key)
obs[fov_key] = pd.Categorical(obs[fov_key].astype(str))
obs[cell_id_key] = obs.index.astype(np.int64)
obs.rename_axis(None, inplace=True)
obs.index = obs.index.astype(str).str.cat(obs[fov_key].values, sep="_")
common_index = obs.index.intersection(counts.index)

adata = AnnData(
    csr_matrix(counts.loc[common_index, :].values),
    dtype=counts.values.dtype,
    obs=obs.loc[common_index, :],
    uns={"spatial": {}},
)

# Vars
adata.var_names = counts.columns

# Global and local coords
adata.obsm["spatial"] = adata.obs[["CenterX_local_px", "CenterY_local_px"]].values
adata.obsm["spatial_fov"] = adata.obs[["CenterX_global_px", "CenterY_global_px"]].values
adata.obs.drop(columns=["CenterX_local_px", "CenterY_local_px"], inplace=True)

# Scalefactors per fov
for fov in adata.obs[fov_key].cat.categories:
    adata.uns["spatial"][fov] = {
        "images": {},
        "scalefactors": {"tissue_hires_scalef": 1, "spot_diameter_fullres": 1},
    }

# Metadata  per fov
if fov_file is not None:
    fov_positions = pd.read_csv(path / fov_file, header=0, index_col=fov_key)
    for fov, row in fov_positions.iterrows():
        if str(fov) in adata.uns["spatial"]:
            adata.uns["spatial"][str(fov)]["metadata"] = row.to_dict()
        else:
            print(f"FOV '{str(fov)}' not in `adata.uns['spatial']`")


sc.write(path / "nanostring_lung5_rep2.h5ad", adata)
df_single_genes = pd.read_csv(
    path / "Lung5_Rep2_tx_file.csv", header=0, index_col=fov_key
)
pd.DataFrame.to_parquet(
    df_single_genes, path / "nanostring_lung5_rep2_single_genes.parquet"
)
