##
# data preparation
import os
import spatialdata as sd
from napari_spatialdata import Interactive
from spatialdata_io import xenium, visium

SPATIALDATA_SANDBOX_PATH = "../../spatialdata-sandbox"
XENIUM_RAW_DATA_PATH = os.path.join(
    SPATIALDATA_SANDBOX_PATH, "xenium_io/data/xenium/outs"
)
XENIUM_SDATA_PATH = os.path.join(SPATIALDATA_SANDBOX_PATH, "xenium_rep1_io/data.zarr")
VISIUM_SDATA_PATH = os.path.join(SPATIALDATA_SANDBOX_PATH, "visium_associated_xenium_io/data.zarr")
LANDMARKS_SDATA_PATH = os.path.join(
    SPATIALDATA_SANDBOX_PATH, "xenium_io/landmarks.zarr"
)

assert os.path.isdir(XENIUM_SDATA_PATH)
assert os.path.isdir(VISIUM_SDATA_PATH)

xenium_sdata = sd.read_zarr(XENIUM_SDATA_PATH)
visium_sdata = sd.read_zarr(VISIUM_SDATA_PATH)

##
# tiles for deep learning

merged = sd.SpatialData(
    images={
        "xenium": xenium_sdata.images["morphology_mip"],
        "visium": visium_sdata.images["CytAssist_FFPE_Human_Breast_Cancer_full_image"],
    },
    shapes={
        "CytAssist_FFPE_Human_Breast_Cancer": visium_sdata.shapes[
            "CytAssist_FFPE_Human_Breast_Cancer"
        ]
    },
    table=visium_sdata.table,
)

min_coordinate = [12790, 12194]
max_coordinate = [15100, 14221]
cropped = merged.query.bounding_box(
    min_coordinate=min_coordinate,
    max_coordinate=max_coordinate,
    axes=["y", "x"],
    target_coordinate_system="aligned",
)
# cropped.add_image('xenium_full', xenium_sdata.images['morphology_mip'])
# cropped.add_image('visium_full', visium_sdata.images['CytAssist_FFPE_Human_Breast_Cancer_full_image'])
# Interactive(cropped)
##
from spatialdata.dataloader.datasets import ImageTilesDataset
from tqdm import tqdm
from spatialdata.transformations import get_transformation, remove_transformation
from spatialdata import transform

circles = merged['CytAssist_FFPE_Human_Breast_Cancer']
t = get_transformation(circles, 'aligned')

transformed_circles =transform(circles, t)
visium_circle_diameter = 2 * transformed_circles.radius.iloc[0]
dataset = ImageTilesDataset(
    sdata=cropped,
    regions_to_images={"CytAssist_FFPE_Human_Breast_Cancer": "xenium"},
    tile_dim_in_units=visium_circle_diameter,
    tile_dim_in_pixels=32,
    target_coordinate_system="aligned",
)
##
tiles = []
for i, tile in enumerate(tqdm(dataset)):
    tiles.append(tile)
##
from spatialdata import SpatialData

full_image = visium_sdata.images["CytAssist_FFPE_Human_Breast_Cancer_full_image"]
all_transformations = get_transformation(full_image, get_all=True)
if "global" in all_transformations:
    remove_transformation(full_image, "global")

tiles_sdata = SpatialData(
    images={"full": full_image}
    | {f"tile_{region_name}_{region_index}": tile for (tile, region_name, region_index) in tiles}
)

##
# Interactive([merged, cropped])
Interactive([tiles_sdata, cropped])
##

# from torch.utils.data import DataLoader
#
# loader = DataLoader(dataset, batch_size=256, num_workers=8)


