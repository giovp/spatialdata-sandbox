##
import numpy as np
import anndata as ad
import shutil
from pathlib import Path
import spatialdata as sd
import re
import os
from tqdm import tqdm
import xarray as xr

# from spatialdata_io import (
#     circles_anndata_from_coordinates,
#     table_update_anndata,
# )

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("visium"):
    path /= "visium"
    assert path.exists()
path_read = path / "data"
path_write = path / "data.zarr"

libraries = [re.sub(".h5ad", "", i) for i in os.listdir(path_read / "tables")][:2]

##
table_list = []
images = {}
points = {}
for lib in tqdm(libraries, desc="loading visium libraries"):
    # prepare table
    table = ad.read_h5ad(path_read / "tables" / f"{lib}.h5ad")
    table.var_names_make_unique()
    table.obs["annotating"] = f"/shapes/{lib}"
    table.obs["library"] = lib
    table.obs["visium_spot_id"] = np.arange(len(table))
    table_list.append(table)

    # setup
    lib_keys = list(table.uns["spatial"].keys())
    assert len(lib_keys) == 1
    lib_key = lib_keys[0]

    # prepare image transformation
    scale_factors = np.array(
        [1.0]
        + [1 / table.uns["spatial"][lib_key]["scalefactors"]["tissue_hires_scalef"]] * 2
    )
    transform = sd.Scale(scale=scale_factors)

    # prepare image
    img = table.uns["spatial"][lib_key]["images"]["hires"]
    assert img.dtype == np.float32 and np.min(img) >= 0.0 and np.max(img) <= 1.0
    scaled = (img * 255).astype(np.uint8)
    scaled = sd.Image2DModel.parse(scaled, transform=transform, dims=("y", "x", "c"))
    images[lib] = scaled

    # prepare circles
    radius = table.uns["spatial"][lib_key]["scalefactors"]["spot_diameter_fullres"] / 2
    shape_regions = sd.ShapesModel.parse(
        coords=table.obsm["spatial"],
        shape_type="Circle",
        shape_size=radius,
        instance_key="visium_spot_id",
        instance_values=np.arange(len(table)),
    )
    points[lib] = shape_regions


table = ad.concat(
    table_list,
    label="library",
    keys=libraries,
)

del table.obsm["spatial"]
adata = sd.TableModel.parse(
    table,
    region=[f"/shapes/{lib}" for lib in libraries],
    region_key="annotating",
    instance_key="visium_spot_id",
)

sdata = sd.SpatialData(
    table=table,
    images=images,
    points=points,
)
print(sdata)

##
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")
print(f'view with "python -m spatialdata view data.zarr"')
sdata = sd.SpatialData.read(path_write)
print(sdata)
print("read")
