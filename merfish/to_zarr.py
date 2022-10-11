##
import json
import numpy as np
import scanpy as sc
import anndata as ad
import imageio
import shutil
from pathlib import Path
import spatialdata as sd
from spatialdata._core import Polygons
import imageio.v3 as iio
import xarray as xr

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("merfish"):
    path /= "merfish"
    assert path.exists()
##
path_read = path / "data" / "processed"
path_write = path / "data.zarr"

##
cells = sc.read_h5ad(path_read / "cells.h5ad")
img = xr.DataArray(iio.imread(path_read / "image.png"), dims=('y', 'x'))
adata = sc.read_h5ad(path_read / "single_molecule.h5ad")
single_molecule = ad.AnnData(shape=(len(adata), 0))
single_molecule.obsm["spatial"] = adata.X
single_molecule.obsm["spatial_type"] = adata.obsm["cell_type"]

j = json.load(open(path_read / "image_transform.json", "r"))
image_translation = np.array([j["translation_y"], j["translation_x"]])
image_scale_factors = np.array([j["scale_factor_y"], j["scale_factor_x"]])

expression = cells.copy()
del expression.obsm["region_radius"]
del expression.obsm["spatial"]
expression.uns["mapping_info"] = {
    "regions": "/points/cells",
    "regions_key": "regions_id",
    "instance_key": "cell_id",
}
expression.obs["regions_id"] = "/points/cells"
expression.obs["cell_id"] = np.arange(len(expression))

regions = ad.AnnData(shape=(len(cells.obsm["spatial"]), 0))
regions.obsm["region_radius"] = cells.obsm["region_radius"]
regions.obsm["spatial"] = cells.obsm["spatial"]
regions.obs["cell_id"] = np.arange(len(regions))

##
adata_polygons = Polygons.anndata_from_geojson(path_read / "anatomical.geojson")

##
translation = sd.Translation(translation=image_translation)
scale = sd.Scale(scale=image_scale_factors)
# these are equivalent
# composed = sd.compose_transformations(scale, translation)
composed = sd.compose_transformations(scale, translation).to_affine()

##
sdata = sd.SpatialData(
    table=expression,
    points={"cells": regions, "single_molecule": single_molecule},
    images={"rasterized": img},
    polygons={"anatomical": adata_polygons},
    transformations={
        ("/images/rasterized", "global"): composed,
    },
)
print(sdata)
##
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")
print(f'view with "python -m napari_spatialdata view data.zarr"')
##
sdata = sd.SpatialData.read(path_write)
print(sdata)
print('read')