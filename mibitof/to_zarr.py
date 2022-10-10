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
import anndata as ad

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("mibitof"):
    path /= "mibitof"
    assert path.exists()
path_read = path / "data"
path_write = path / "data.zarr"

##
libraries = ["point8", "point16", "point23"]

table_list = []
for lib in libraries:
    table = ad.read(path_read / f"{lib}_table.h5ad")
    table.obs["library_id"] = f'/labels/{lib}'
    table.obs["cell_id"] = np.arange(len(table))
    table_list.append(table)

table = ad.concat(
    table_list,
    keys=libraries,
)
table.uns["mapping_info"] = {
    "regions": [f'/labels/{lib}' for lib in libraries],
    "regions_key": "library_id",
    "instance_key": "cell_id",
}
##

sdata = sd.SpatialData(
    table=table,
    labels={lib: imageio.imread(path_read / f"{lib}_labels.png") for lib in libraries},
    images={
        lib: np.moveaxis(imageio.imread(path_read / f"{lib}_image.png"), 2, 0)
        for lib in libraries
    },
    images_axes={lib: ("c", "y", "x") for lib in libraries},
    labels_axes={lib: ("y", "x") for lib in libraries},
    transformations={(f"/images/{lib}", lib): None for lib in libraries}
    | {(f"/labels/{lib}", lib): None for lib in libraries},
    coordinate_systems=[
        {
            "name": lib,
            "axes": [
                {"name": "c", "type": "channel"},
                {"name": "x", "type": "space", "unit": "micrometer"},
                {"name": "y", "type": "space", "unit": "micrometer"},
            ],
        }
        for lib in libraries
    ],
)
print(sdata)

##
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")
print(f'view with "python -m spatialdata view data.zarr"')
##
sdata = sd.SpatialData.read(path_write)
print(sdata)
print('read')