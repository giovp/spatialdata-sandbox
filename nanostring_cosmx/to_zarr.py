##
import spatialdata as sd
import shutil
import anndata as ad
from pathlib import Path
import pandas as pd
import numpy as np
import imageio
from tqdm import tqdm
import time
import json
from ome_zarr.io import parse_url
import zarr
from ome_zarr.io import parse_url
from ome_zarr.writer import write_image, write_labels
import zarr
from anndata.experimental import write_elem

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

# read points
# takes 20 seconds
start = time.time()
points_df = pd.read_parquet(path_read / "nanostring_lung5_rep2_single_genes.parquet")
print(f'read_parquet(): {time.time() - start}')

##
points_df.reset_index(inplace=True, drop=False)
# takes 50 seconds
start = time.time()
points_df["fov"] = pd.Categorical(points_df["fov"].astype(str))
print(f'created categorical col: {time.time() - start}')
points_df.reset_index(inplace=True, drop=False)

##
cols = ["y_local_px", "x_local_px"]
# takes 10 seconds
start = time.time()
xy = points_df[["x_local_px", "y_local_px"]].to_numpy()
print(f'xy: {time.time() - start}')

##
# takes 80 seconds
start = time.time()
points = ad.AnnData(
    shape=(len(xy), 0),
    obsm={"spatial": xy},
    obs=points_df[points_df.columns.difference(cols)],
    dtype=np.float_,
)
print(f'creating anndata: {time.time() - start}')
##
list_of_images = []
list_of_labels = []
list_of_points = []
for fov in tqdm(categories, desc="fovs"):
    sfov = fov.zfill(3)

    image = imageio.imread(path_read / "CellComposite" / f"CellComposite_F{sfov}.jpg")
    labels = imageio.imread(path_read / "CellLabels" / f"CellLabels_F{sfov}.tif")

    list_of_images.append(np.flipud(image))
    list_of_labels.append(labels)
    # write_image(
    #     image=image,
    #     group=image_group,
    #     axes=["c", "y", "x"],
    #     scaler=None,
    # )
    # print(f"Written image `{fov}` .")
    #
    # write_labels(
    #     labels=labels,
    #     group=image_group,
    #     name=fov,
    #     axes=["y", "x"],
    #     scaler=None,
    # )
    # print(f"Written labels `{fov}` .")

    points_subset = points[points.obs.fov == fov].copy()
    # write_table_points(
    #     group=image_group,
    #     adata=points[points.obs.fov == fov].copy(),
    # )
    # print(f"Written points `{fov}` .")
    list_of_points.append(points_subset)
##
sdata = sd.SpatialData(
    labels={
        fov: labels for fov, labels in zip(categories, list_of_labels)
    },
    images={
        fov: image for fov, image in zip(categories, list_of_images)
    },
    points={
        fov: points_subset
        for fov, points_subset in zip(categories, list_of_points)
    },
)
print(sdata)
##
if path_write.exists():
    shutil.rmtree(path_write)

# takes 45 seconds, 3.9 GB written to disk (raw data is 5.1 GB)
start = time.time()
sdata.write(path_write)
print(f'saving .zarr file: {time.time() - start}')

print(f'view with "python -m spatialdata view data.zarr"')

##
my_labels={
    fov: labels for fov, labels in zip(categories, list_of_labels) if fov in ['1', '2']
}
my_images={
    fov: image for fov, image in zip(categories, list_of_images) if fov in ['1', '2']
}
my_points={
    fov: points_subset
    for fov, points_subset in zip(categories, list_of_points) if fov in ['1', '2']
}
sdata = sd.SpatialData(labels=my_labels, images=my_images, points=my_points)

if path_write_small.exists():
    shutil.rmtree(path_write_small)

# takes 1 second, 0.1 GB written to disk
start = time.time()
sdata.write(path_write_small)
print(f'saving .zarr file: {time.time() - start}')

print(f'view with "python -m spatialdata view data_small.zarr"')

print("done")
