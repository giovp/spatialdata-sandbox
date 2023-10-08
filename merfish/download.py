##
# the best is to run this .py file as a script, if you want you can also to convert it to a jupyter notebook with
# jupytext --to notebook merfish.py
# When runnning the notebook the working directory should be the project root folder otherwise the data paths will be
# not found (use os.chdir(<project root path>)

##
import pandas as pd
import numpy as np
import anndata as ad
import json
import matplotlib
import matplotlib.pyplot as plt
import os
import shutil
import math
import datashader
import imageio
from pathlib import Path
import sys
from PIL import Image
sys.path.insert(1, os.path.join(sys.path[0], Path(__file__).parent.parent.resolve()))

from utils import download

plt.style.use("dark_background")

# PLOT = True
PLOT = False

##
print(f"os.getcwd() = {os.getcwd()}")
merfish_dir = Path().resolve() / "data"
merfish_dir.mkdir(parents=True, exist_ok=True)
assert merfish_dir.exists()
raw = merfish_dir / "raw"
raw.mkdir(parents=True, exist_ok=True)
processed = merfish_dir / "processed"
processed.mkdir(parents=True, exist_ok=True)

##
def my_download(url):
    return
    # aria2c asks for relative paths
    download(
        url,
        (raw / os.path.basename(url)).relative_to(Path().resolve()),
        desc="data",
    )


for url in [
    "https://github.com/spacetx-spacejam/data/blob/master/gene_lists/MERFISH_genes.csv",
    "https://s3.amazonaws.com/starfish.data.spacetx/spacejam2/MERFISH_Allen_VISp/Allen_MERFISH_spots_with_anatomy.csv",
    "https://s3.amazonaws.com/starfish.data.spacetx/spacejam2/MERFISH_Allen_VISp/fixed_1001844875.csv",
    "https://raw.githubusercontent.com/spacetx-spacejam/data/master/annotations/Allen_MERFISH_Layers.geojson",
]:
    my_download(url)


##
class BoundingBox:
    pass


bb = BoundingBox()
bb.x0 = 1154
bb.x1 = 3172
bb.y0 = 4548
bb.y1 = 6566


##
# extract single molecule (points) data
df = pd.read_csv(raw / "Allen_MERFISH_spots_with_anatomy.csv")
plt.figure(figsize=(10, 10))
points_bb = bb
df = df[
    (df["x_um"] > points_bb.x0)
    & (df["x_um"] < points_bb.x1)
    & (df["y_um"] > points_bb.y0)
    & (df["y_um"] < points_bb.y1)
]
xy = df[["x_um", "y_um"]].to_numpy()
a_points = ad.AnnData(xy)
# a_points.obsm["spatial"] = xy
a_points.obsm["cell_type"] = df["layer"].to_numpy()
# a_points.obs.columns = ["0"]
a_points.write_h5ad(processed / "single_molecule.h5ad")

##
if PLOT:
    plt.figure()
    plt.scatter(xy[:, 0], xy[:, 1], s=1)
    plt.gca().set_aspect("equal")
    plt.show()

##
def plot_single_molecule_anndata(adata: ad.AnnData, ax=None):
    xy = adata.X
    c = adata.obsm["cell_type"]
    if ax is None:
        plt.figure(figsize=(10, 10))
        cax = plt.gca()
    else:
        cax = ax
    all_types = list(set(c))
    cax.scatter(xy[:, 0], xy[:, 1], s=1, c=[all_types.index(cc) for cc in c])
    cax.set_aspect("equal")
    if ax is None:
        plt.show()


if PLOT:
    plot_single_molecule_anndata(a_points)

##
# extract single cell masks data
df = pd.read_csv(raw / "fixed_1001844875.csv")
df.drop(columns=[df.columns[0], df.columns[1]], inplace=True)
df["radius"] = df["area"].apply(lambda x: math.sqrt(x / math.pi))
genes = df.columns.tolist()
genes = genes[: genes.index("area")]
xy = df[["x_um", "y_um"]].to_numpy()

a_cells = ad.AnnData(X=df[genes])
a_cells.obsm["spatial"] = xy
a_cells.obsm["region_radius"] = df["radius"].to_numpy()
a_cells.write_h5ad(processed / "cells.h5ad")

##
def plot_shape_masks_anndata(adata: ad.AnnData, ax=None):
    if ax is None:
        plt.figure(figsize=(10, 10))
        cax = plt.gca()
    else:
        cax = ax
    xy = adata.obsm["spatial"]
    radius = adata.obsm["region_radius"]
    patches = []
    for (x, y), r in zip(xy, radius):
        if bb.x0 - r < x < bb.x1 + r and bb.y0 - r < y < bb.y1 + r:
            patch = matplotlib.patches.Circle(
                (x, y), r, color=np.append(np.random.rand(3), 0.3)
            )
            patches.append(patch)
    p = matplotlib.collections.PatchCollection(patches, match_original=True)
    cax.add_collection(p)
    cax.set(xlim=(bb.x0, bb.x1), ylim=(bb.y0, bb.y1))
    cax.set_aspect("equal")
    if ax is None:
        plt.show()


if PLOT:
    plot_shape_masks_anndata(a_cells)

##

# generate the raster image
df = pd.read_csv(raw / "Allen_MERFISH_spots_with_anatomy.csv")
df["datashader"] = np.array([1] * len(df))
raster_w = 600
raster_h = 600
cvs = datashader.Canvas(plot_width=raster_w, plot_height=raster_h)
# agg = cvs.points(df, x="x_um", y="y_um", agg=datashader.any())
agg = cvs.points(df, x="x_um", y="y_um", agg=datashader.count())
# img = datashader.tf.shade(agg)
# raster = img.to_numpy()
raster = agg.values
raster = raster.astype(np.float64)
raster /= raster.max()
# raster = np.flipud(raster)
# raster = np.log(1 + raster)

if PLOT:
    plt.figure()
    len(df)
    # plt.hist(raster.flatten())
    plt.hist(raster)
    # plt.hist(raster.flatten()[raster.flatten() >1000], bins=1000)
    plt.show()

##
# let's manually adjust the levels to make the image brighter in this example
raster = np.clip(raster, a_min=0.0, a_max=0.2)
raster *= 5
#
# plt.figure()
# plt.imshow(raster, origin='lower')
# plt.xlim([400, 425])
# plt.ylim([175, 200])
# plt.show()


##
if PLOT:
    plt.figure()
    plt.imshow(
        raster,
        extent=(df["x_um"].min(), df["x_um"].max(), df["y_um"].min(), df["y_um"].max()),
        origin="lower",
    )
    plt.scatter(df["x_um"], df["y_um"], s=1, alpha=0.01)
    plt.show()

##
min_x = df["x_um"].min()
min_y = df["y_um"].min()
max_x = df["x_um"].max()
max_y = df["y_um"].max()

raster_bb = BoundingBox()
raster_bb.x0 = int((bb.x0 - min_x) / (max_x - min_x) * raster_w)
raster_bb.x1 = int((bb.x1 - min_x) / (max_x - min_x) * raster_w)
raster_bb.y0 = int((bb.y0 - min_y) / (max_y - min_y) * raster_h)
raster_bb.y1 = int((bb.y1 - min_y) / (max_y - min_y) * raster_h)
print(
    f"raster_bb.x0 = {raster_bb.x0}, raster_bb.x1 = {raster_bb.x1}, raster_bb.y0 = {raster_bb.y0}, raster_bb.y1 = {raster_bb.y1}"
)

# raster_crop = np.flipud(np.flipud(raster)[raster_bb.y0 : raster_bb.y1, raster_bb.x0 : raster_bb.x1])
raster_crop = raster[raster_bb.y0 : raster_bb.y1, raster_bb.x0 : raster_bb.x1]
print(raster_crop.shape)

if PLOT:
    plt.figure()
    plt.imshow(raster_crop, origin="lower")
    plt.show()
##
translation = np.array([bb.x0, bb.y0])
# assert bb.x1 - bb.x0 == bb.y1 - bb.y0
# assert raster_w == raster_h
scale_factor_x = (max_x - min_x) / raster_w
scale_factor_y = (max_y - min_y) / raster_h
scale_factors = np.array([scale_factor_x, scale_factor_y])
# wrong
# scale_factor = (bb.x1 - bb.x0) / raster_w
##
raster_crop = (raster_crop * 255).astype(np.uint8)
imageio.imwrite(processed / "image.png", raster_crop)
# np.save(os.path.join(output_dir, "image"), raster_crop)

d = {}
d["translation_x"] = float(translation[0])
d["translation_y"] = float(translation[1])
d["scale_factor_x"] = scale_factor_x
d["scale_factor_y"] = scale_factor_y
with open(processed / "image_transform.json", "w") as outfile:
    json.dump(d, fp=outfile)

##
print(f"translation = {translation}, scale_factor = {scale_factors}")
x0 = translation[0]
y0 = translation[1]
x1 = translation[0] + raster.shape[1] * scale_factor_x
y1 = translation[1] + raster.shape[0] * scale_factor_y
extent = (x0, x1, y0, y1)
print(extent)

if PLOT:
    plt.figure()
    plt.imshow(raster_crop, extent=extent, origin="lower")
    plt.show()
##


# first translation, then scaling
def plot_raster(
    raster: np.ndarray, translation: np.array, scale_factors: np.array, ax=None
):
    assert len(translation) == 2
    # grayscale image or rgb/rgba
    assert (
        len(raster.shape) == 2 or len(raster.shape) == 3 and raster.shape[2] in [3, 4]
    )
    if ax is None:
        plt.figure(figsize=(10, 10))
        cax = plt.gca()
    else:
        cax = ax
    x0 = translation[0]
    y0 = translation[1]
    x1 = translation[0] + raster.shape[1] * scale_factors[0]
    y1 = translation[1] + raster.shape[0] * scale_factors[1]
    extent = (x0, x1, y0, y1)
    cax.imshow(
        raster_crop,
        extent=extent,
        cmap=plt.cm.get_cmap("gray"),
        alpha=0.4,
        origin="lower",
    )
    if ax is None:
        plt.show()


if PLOT:
    plot_raster(raster, translation=translation, scale_factors=scale_factors)
##
# polygon information
layers = json.load(open(raw / "Allen_MERFISH_Layers.geojson", "r"))
shutil.copyfile(raw / "Allen_MERFISH_Layers.geojson", processed / "anatomical.geojson")
layers

##
# there is a bug here and things are not plotted aligned, but things are aligned in the napari viewer, that is our
# real goal so we will fix this plot when implementing more general plotting functions
if PLOT:
    plt.figure(figsize=(10, 10))
    ax = plt.gca()
    plot_single_molecule_anndata(a_points, ax)
    plot_shape_masks_anndata(a_cells, ax)
    plot_raster(raster, translation=translation, scale_factors=scale_factors, ax=ax)
    plt.show()

##
brain_layers = {}
for layer in layers["geometries"]:
    assert layer["type"] == "Polygon"
    name = layer["name"]
    coordinates = np.array(layer["coordinates"])
    coordinates = np.squeeze(coordinates, 0)
    brain_layers[name] = coordinates

# old way of storing polygon information
# a_polygon = ad.AnnData(None, obs=list(brain_layers.keys()))
# a_polygon.obs.columns = ["layer"]
# a_polygon.obs["vertices"] = list(brain_layers.values())
# # temporary, inefficient arbitrary way of storing the coordinates
# a_polygon.obs["vertices"] = a_polygon.obs["vertices"].apply(lambda x: repr(x))
# a_polygon.write_h5ad(os.path.join(output_dir, "polygons.h5ad"))

##
if PLOT:
    plt.figure()
    for layer, coordinates in brain_layers.items():
        plt.plot(coordinates[:, 0], coordinates[:, 1])
    plt.gca().set_aspect("equal")
    plt.legend(
        [layer for layer in brain_layers.keys()],
        loc="upper center",
        bbox_to_anchor=(0.5, 1.25),
        ncol=3,
    )
    plt.tight_layout()
    plt.show()
