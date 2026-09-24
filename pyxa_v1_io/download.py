##
import os
import subprocess
import zipfile
from pathlib import Path

# from https://huggingface.co/datasets/Stellaromics/demo (BSD 3-Clause)
# xsmall/: a 100 x 100 x 100 um cube cropped from the full demo region (small/), ~8 MB

BASE_URL = "https://huggingface.co/datasets/Stellaromics/demo/resolve/main/xsmall"
FILES = [
    "cell_assigned_gene_v1.csv",
    "cell_by_gene_v1.csv",
    "cell_metadata_v1.csv",
    "segmentation_geometries_v1.parquet",
    "mosaic_3d.ome.zarr.zip",
]

data_dir = Path(__file__).resolve().parent / "data"
os.makedirs(data_dir, exist_ok=True)

##
# download the data
for filename in FILES:
    command = f"curl -L -C - -o {data_dir / filename} {BASE_URL}/{filename}"
    subprocess.run(command, shell=True, check=True)

##
# unzip the DAPI mosaic (a zipped OME-Zarr store containing mosaic_3d.ome.zarr/)
with zipfile.ZipFile(data_dir / "mosaic_3d.ome.zarr.zip") as zf:
    zf.extractall(data_dir)
os.remove(data_dir / "mosaic_3d.ome.zarr.zip")
