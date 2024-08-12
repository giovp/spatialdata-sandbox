# this script helps a user read his seqfish data, which are formatted in a newer version that introduced some changes to the initial format
# the labels are not parsed because the download of the labels data took too long so I could not test it
# see discussion here: https://github.com/scverse/spatialdata-io/issues/193

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from dask_image.imread import imread
from spatialdata import SpatialData
from spatialdata.models import (
    Image2DModel,
    Labels2DModel,
    PointsModel,
    ShapesModel,
    TableModel,
)
from spatialdata.transformations import Identity

from spatialdata_io._constants._constants import SeqfishKeys as SK

__all__ = ["seqfish"]


from spatialdata_io import seqfish
from pathlib import Path

path = Path("/Users/macbook/Downloads")

prefix = "TestSet"


def get_cell_file() -> str:
    return f"{prefix}_{SK.CELL_COORDINATES}{SK.CSV_FILE}"


def get_count_file() -> str:
    # SK.COUNTS_FILE is "CxG", now it is CellxGene
    return f"{prefix}_CellxGene{SK.CSV_FILE}"


def get_dapi_file() -> str:
    # it was OME_TIFF, now it's TIFF
    return f"{prefix}_{SK.DAPI}{SK.TIFF_FILE}"


def get_cell_mask_file() -> str:
    return f"{prefix}_{SK.CELL_MASK_FILE}{SK.TIFF_FILE}"


def get_transcript_file() -> str:
    # SK.TRANSCRIPT_COORDINATES is "TranscriptCoordinates", now it is TranscriptList
    return f"{prefix}_TranscriptList{SK.CSV_FILE}"


cell_file = get_cell_file()
count_matrix = get_count_file()
adata = ad.read_csv(path / count_matrix, delimiter=",")
cell_info = pd.read_csv(path / cell_file, delimiter=",")
adata.obsm[SK.SPATIAL_KEY] = cell_info[[SK.CELL_X, SK.CELL_Y]].to_numpy()
adata.obs[SK.AREA] = np.reshape(cell_info[SK.AREA].to_numpy(), (-1, 1))
region = f"cells"
adata.obs[SK.REGION_KEY] = region
adata.obs[SK.INSTANCE_KEY_TABLE] = adata.obs.index.astype(int)

scale_factors = [2, 2, 2, 2]

images = {
    f"image": Image2DModel.parse(
        imread(path / get_dapi_file()),
        dims=("c", "y", "x"),
        scale_factors=scale_factors,
        transformations={"global": Identity()},
    )
}

# labels = {
#     f"labels": Labels2DModel.parse(
#         imread(path / get_cell_mask_file()).squeeze(),
#         dims=("y", "x"),
#         scale_factors=scale_factors,
#         transformations={"global": Identity()},
#     )
# }
labels = {}

# there are 4 z values in the transcript files, if needed
df = pd.read_csv(path / get_transcript_file(), delimiter=",")
print(df.z.value_counts())

points = {
    f"transcripts": PointsModel.parse(
        df,
        # adding the z coordinate which is now available
        coordinates={"x": SK.TRANSCRIPTS_X, "y": SK.TRANSCRIPTS_Y, "z": "z"},
        feature_key=SK.FEATURE_KEY.value,
        # 'cell' column is not available anymore
        # instance_key=SK.INSTANCE_KEY_POINTS.value,
        transformations={"global": Identity()},
    )
}

adata.obs[SK.REGION_KEY] = adata.obs[SK.REGION_KEY].astype("category")
adata.obs = adata.obs.reset_index(drop=True)
table = TableModel.parse(
    adata,
    region=f"cells",
    region_key=SK.REGION_KEY.value,
    instance_key=SK.INSTANCE_KEY_TABLE.value,
)

shapes = {
    f"cells": ShapesModel.parse(
        adata.obsm[SK.SPATIAL_KEY],
        geometry=0,
        radius=np.sqrt(adata.obs[SK.AREA].to_numpy() / np.pi),
        index=adata.obs[SK.INSTANCE_KEY_TABLE].copy(),
        transformations={"global": Identity()},
    )
}

sdata = SpatialData(images=images, labels=labels, points=points, tables={'table': table}, shapes=shapes)

print('writing... ', end='')
sdata.write("/Users/macbook/Downloads/seqfish.zarr")
print('DONE')
sdata = SpatialData.read("/Users/macbook/Downloads/seqfish.zarr")

from napari_spatialdata import Interactive
Interactive(sdata)

