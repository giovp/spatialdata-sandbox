##
# data preparation
import os
import spatialdata as sd
from spatialdata._core.transformations import Sequence
from napari_spatialdata import Interactive
from spatialdata_io import xenium, visium

SPATIALDATA_SANDBOX_PATH = "../../spatialdata-sandbox"
XENIUM_RAW_DATA_PATH = os.path.join(
    SPATIALDATA_SANDBOX_PATH, "xenium_io/data/xenium/outs"
)
VISIUM_RAW_DATA_PATH = os.path.join(SPATIALDATA_SANDBOX_PATH, "xenium_io/data/visium")
XENIUM_SDATA_PATH = os.path.join(SPATIALDATA_SANDBOX_PATH, "xenium_io/data.zarr")
VISIUM_SDATA_PATH = os.path.join(SPATIALDATA_SANDBOX_PATH, "xenium_io/data_visium.zarr")
LANDMARKS_SDATA_PATH = os.path.join(
    SPATIALDATA_SANDBOX_PATH, "xenium_io/landmarks.zarr"
)

assert os.path.isdir(XENIUM_RAW_DATA_PATH)
assert os.path.isdir(VISIUM_RAW_DATA_PATH)


ALREADY_IN_ZARR = True

if ALREADY_IN_ZARR:
    xenium_sdata = sd.read_zarr(XENIUM_SDATA_PATH)
else:
    # here we read the raw data into spatialdata objects
    print("reading the xenium data... ", end="")
    xenium_sdata = xenium(XENIUM_RAW_DATA_PATH, cells_as_shapes=True)
    print("done")

    # for large datasets, like in this case, saving the data to zarr will automatically improve the performance of
    # downstream processing thanks to file backing and chunked data storage
    print("saving the xenium data to zarr... ", end="")
    xenium_sdata.write(XENIUM_SDATA_PATH, overwrite=True)
    print("done")

print(xenium_sdata)


if ALREADY_IN_ZARR:
    visium_sdata = sd.read_zarr(VISIUM_SDATA_PATH)
else:
    print("reading the visium data... ", end="")
    visium_sdata = visium(VISIUM_RAW_DATA_PATH)
    print("done")

    print("saving the visium data to zarr... ", end="")
    visium_sdata.write(VISIUM_SDATA_PATH, overwrite=True)
    print("done")

print(visium_sdata)

merged = sd.SpatialData(
    images={
        "xenium": xenium_sdata.images["morphology_mip"],
        "visium": visium_sdata.images["CytAssist_FFPE_Human_Breast_Cancer_full_image"],
    }
)

# let's save them to Zarr so that we can use in the future if needed
ALREADY_IN_ZARR = True
if ALREADY_IN_ZARR:
    landmarks_sdata = sd.read_zarr(LANDMARKS_SDATA_PATH)
else:
    landmarks_sdata = sd.SpatialData(shapes=merged.shapes)
    landmarks_sdata.write(LANDMARKS_SDATA_PATH, overwrite=True)


from spatialdata._core._transform_elements import align_elements_using_landmarks

affine = align_elements_using_landmarks(
    references_coords=landmarks_sdata.shapes["xenium_landmarks"],
    moving_coords=landmarks_sdata.shapes["visium_landmarks"],
    reference_element=merged.images["xenium"],
    moving_element=merged.images["visium"],
    reference_coordinate_system="global",
    moving_coordinate_system="global",
    new_coordinate_system="aligned",
)

from spatialdata._core._spatialdata_ops import (
    get_transformation,
    set_transformation,
    remove_transformation,
)

for element in [
    visium_sdata.images["CytAssist_FFPE_Human_Breast_Cancer_full_image"],
    visium_sdata.shapes["CytAssist_FFPE_Human_Breast_Cancer"],
]:
    transformation = get_transformation(element, "global")
    new_transformation = Sequence([transformation, affine])
    set_transformation(element, new_transformation, "aligned")

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
from spatialdata._dl.datasets import ImageTilesDataset
from tqdm import tqdm

dataset = ImageTilesDataset(
    sdata=cropped,
    regions_to_images={"CytAssist_FFPE_Human_Breast_Cancer": "xenium"},
    tile_dim_in_units=150,
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


