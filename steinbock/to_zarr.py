
import os
os.environ['USE_PYGEOS'] = '0'
import json
import numpy as np
import scanpy as sc
import shutil
from pathlib import Path
import spatialdata as sd
import imageio.v3 as iio
import pandas as pd
from glob import glob

print(f"os.getcwd() = {os.getcwd()}")
data_dir = Path().resolve() / "data"
data_dir.mkdir(parents=True, exist_ok=True)
assert data_dir.exists()
raw = data_dir / "steinbock"
raw.mkdir(parents=True, exist_ok=True)
processed = data_dir / "processed"
processed.mkdir(parents=True, exist_ok=True)

raw_images = { f.stem.split('.')[0]: iio.imread(f)  for f in (raw / "ome").glob("*.tiff") } 
marker_list = pd.read_csv(raw / "panel.csv")
image_list = pd.read_csv(raw / "images.csv")

raw_deepcell = { f.stem.split('.')[0]: iio.imread(f)  for f in (raw / "masks_deepcell").glob("*.tiff") } 
raw_ilastik = { f.stem.split('.')[0]: iio.imread(f)  for f in (raw / "masks_ilastik").glob("*.tiff") } 


images = { k: sd.Image2DModel.parse(v, dims=("c", "y", "x"))  for k, v in raw_images.items() }

# create labels for each sample and both segmentation methods
deepcell = { f'{k}_deepcell': sd.Labels2DModel.parse(v, dims=("y", "x"))  for k, v in raw_deepcell.items() }
ilastik = { f'{k}_ilastik': sd.Labels2DModel.parse(v, dims=("y", "x"))  for k, v in raw_ilastik.items() }

adata = sc.read_h5ad(raw / "cells.h5ad")

adata.obs = (
    adata.obs
    .assign(cell_id=lambda df: df.index.to_series().apply(lambda r: int(r.split(' ')[1]) ))
    .assign(library_id=lambda df: df['image'].apply(lambda r: f"/labels/{r.rstrip('.tiff')}_deepcell" )
    )
)

table = sd.TableModel.parse(
    adata = adata,
    region = [ f'/labels/{n}' for n in list(deepcell.keys()) + list(ilastik.keys()) ],
    region_key = 'library_id',
    instance_key = 'cell_id',
)

sdata = sd.SpatialData(
    table=table,
    images=images,
    labels={**deepcell, **ilastik}
)
 
path_write = processed / "data.zarr"

if path_write.exists():
    shutil.rmtree(path_write)

sdata.write(path_write)
