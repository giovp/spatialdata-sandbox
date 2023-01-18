##
import pyarrow as pa
import os
import anndata as ad
os.environ['USE_PYGEOS'] = '0'
import json
import numpy as np
import scanpy as sc
import shutil
from pathlib import Path
import spatialdata as sd
import imageio.v3 as iio
import pandas as pd


print(f"os.getcwd() = {os.getcwd()}")
data_dir = Path().resolve() / "data"
data_dir.mkdir(parents=True, exist_ok=True)
assert data_dir.exists()
raw = data_dir / "raw"
raw.mkdir(parents=True, exist_ok=True)
processed = data_dir / "processed"
processed.mkdir(parents=True, exist_ok=True)


image = iio.imread(raw / 'lymphnode.tif')
segmentation = iio.imread(raw / 'segmentation.tif')
annotation = pd.read_csv(raw / 'annotations.csv', index_col=0)
markerlist = pd.read_csv(raw / 'markerlist.csv', index_col=0)
expression = pd.read_csv(raw / 'expression.csv', index_col=0)

# construct the spatial data object (ln stands for lymph node)
images = { 'ln': sd.Image2DModel.parse(image) }
labels = { 'ln': sd.Labels2DModel.parse(segmentation,  dims=("y", "x")) }

# creates the anndata object for the spatial data table
adata = ad.AnnData(
    X=expression.values, 
    obs=annotation,
    var=pd.DataFrame(expression.columns, columns=['marker']).set_index('marker'),
    dtype=expression.values.dtype
)

adata.obs['library_id'] = f'/labels/ln'
table = sd.TableModel.parse(adata, region='labels/ln', instance_key='cell_id')

sdata = sd.SpatialData(table=table, labels=labels, images=images)
sdata.images['ln'].coords['c'] = np.array(markerlist['marker'].tolist())


path_write = processed / "data.zarr"

if path_write.exists():
    shutil.rmtree(path_write)

sdata.write(path_write)
