##
import spatialdata as sd
import shutil
import anndata as ad
from pathlib import Path
import numpy as np
from tqdm import tqdm
import time
import pyarrow as pa
import pyarrow.parquet as pq
import imageio.v3 as iio

path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("nanostring_cosmx"):
    path /= "nanostring_cosmx"
    assert path.exists()

path_read = path / "data" / "data_lung5_rep2"
path_write = path / "data.zarr"
path_write_small = path / "data_small.zarr"

##
# read table
table = ad.read(path_read / "nanostring_lung5_rep2.h5ad")
categories = table.obs.fov.cat.categories.to_list()

##
pp = pq.read_table(path_read / "nanostring_lung5_rep2_single_genes.parquet")
print(f"read {len(pp)} points")

##
# let's reduce the dataset size
categories = categories[:3]

##
list_of_images = []
list_of_labels = []
list_of_points = []
for fov in tqdm(categories, desc="images, labels, points"):
    sfov = fov.zfill(3)

    image = iio.imread(path_read / "CellComposite" / f"CellComposite_F{sfov}.jpg")
    labels = iio.imread(path_read / "CellLabels" / f"CellLabels_F{sfov}.tif")

    list_of_images.append(sd.Image2DModel.parse(np.flipud(image), dims=("y", "x", "c")))
    list_of_labels.append(sd.Labels2DModel.parse(np.flipud(labels), dims=("y", "x")))

    # subsetting points
    df = pp.filter(pa.compute.equal(pp.column("fov"), int(fov))).to_pandas()
    xy = df[["x_local_px", "y_local_px"]].to_numpy()
    points = sd.PointsModel.parse(coords=xy)
    list_of_points.append(points)

##
list_of_circles = []
list_of_circles_transforms = []
for fov in tqdm(categories, desc="circles"):
    cells = table[table.obs.fov == fov]
    xy = cells.obsm["spatial"]
    radii = np.sqrt(cells.obs["Area"].to_numpy() / np.pi)
    circles = sd.ShapesModel.parse(
        coords=xy,
        shape_type="Circle",
        shape_size=radii,
    )
    list_of_circles.append(circles)

table = table[table.obs.fov.isin(categories)]
table.obs["fov_path"] = table.obs["fov"].apply(lambda x: f"/shapes/cells{x}")
table = sd.TableModel.parse(
    table,
    region=[f"/shapes/cells{fov}" for fov in categories],
    region_key="fov_path",
    instance_key="cell_ID",
)
##
# name the coordinate systems after the sample name
# TODO: make a API for doing this kind of operations, maybe an argument in the parser to specify the coordinate system name
def rename_coordinate_syestem(element: sd.SpatialElement, new_name: str):
    t = sd.get_transform(element)
    t.output_coordinate_system._name = new_name
    new_element = sd.set_transform(element, t)
    return new_element

for fov in categories:
    for list_of_elements in [list_of_images, list_of_labels, list_of_circles, list_of_points]:
        element = list_of_elements[categories.index(fov)]
        new_element = rename_coordinate_syestem(element, fov)
        list_of_elements[categories.index(fov)] = new_element
##
sdata = sd.SpatialData(
    labels={fov: labels for fov, labels in zip(categories, list_of_labels)},
    images={fov: image for fov, image in zip(categories, list_of_images)},
    points={
        f"points{fov}": points_subset
        for fov, points_subset in zip(categories, list_of_points)
    },
    shapes={f"cells{fov}": circles for fov, circles in zip(categories, list_of_circles)},
    # transformations={(f"/images/{fov}", fov): None for fov in categories}
    # | {(f"/labels/{fov}", fov): None for fov in categories}
    # | {(f"/points/points{fov}", fov): None for fov in categories}
    # | {(f"/points/cells{fov}", fov): None for fov in categories},
    table=table,
)
print(sdata)

##
if path_write.exists():
    shutil.rmtree(path_write)

# takes 45 seconds, 3.9 GB written to disk (raw data is 5.1 GB)
start = time.time()
sdata.write(path_write)
print(f"saving .zarr file: {time.time() - start}")
print(f'view with "python -m spatialdata view data.zarr"')

##
# sdata = sd.SpatialData.read(path_write, filter_table=True)
sdata = sd.SpatialData.read(path_write)
print(sdata)
print("read")
