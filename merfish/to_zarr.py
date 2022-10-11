##
import json
import numpy as np
import scanpy as sc
import shutil
from pathlib import Path
import spatialdata as sd
import imageio.v3 as iio
import xarray as xr
from spatialdata_io import (
    points_anndata_from_coordinates,
    polygons_anndata_from_geojson,
    circles_anndata_from_coordinates,
    table_update_anndata,
)

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
img = xr.DataArray(iio.imread(path_read / "image.png"), dims=("y", "x"))
adata = sc.read_h5ad(path_read / "single_molecule.h5ad")

##
single_molecule = points_anndata_from_coordinates(coordinates=adata.X, points_types=adata.obsm["cell_type"])

expression = cells.copy()
del expression.obsm["region_radius"]
del expression.obsm["spatial"]
table_update_anndata(
    adata=expression,
    regions="/points/cells",
    regions_key="regions_id",
    instance_key="cell_id",
    regions_values="/points/cells",
    instance_values=np.arange(len(cells)),
)
xy = cells.obsm["spatial"]
regions = circles_anndata_from_coordinates(
    coordinates=xy,
    radii=cells.obsm["region_radius"],
    instance_key="cell_id",
    instance_values=np.arange(len(xy)),
)

adata_polygons = polygons_anndata_from_geojson(path_read / "anatomical.geojson")

##
j = json.load(open(path_read / "image_transform.json", "r"))
image_translation = np.array([j["translation_y"], j["translation_x"]])
image_scale_factors = np.array([j["scale_factor_y"], j["scale_factor_x"]])

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
print("read")
