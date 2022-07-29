##
import shutil

import json
import numpy as np
import scanpy as sc
import os
import imageio
import shutil
from pathlib import Path
import spatialdata as sd
import ngff_tables_prototype as viz
import anndata as ad

##
print(f"os.getcwd() = {os.getcwd()}")
path = Path().resolve() / "mibitof"
assert path.exists()
path_read = path / "data"
path_write = path / "data.zarr"

##
libraries = ["point8", "point16", "point23"]

table = ad.concat(
    [ad.read(path_read / f"{lib}_table.h5ad") for lib in libraries],
    label="fov",
    keys=libraries,
)
##

sdata = sd.SpatialData(
    adata=table,
    regions={lib: imageio.imread(path_read / f"{lib}_labels.png") for lib in libraries},
    images={lib: imageio.imread(path_read / f"{lib}_image.png") for lib in libraries},
    tables_region=libraries,
    tables_region_key="fov",
    tables_instance_key='cell_id',
)
print(sdata)

##
if path_write.exists():
    shutil.rmtree(path_write)
sdata.to_zarr(path_write)
print('done')
##

