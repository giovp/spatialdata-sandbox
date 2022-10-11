import matplotlib.pyplot as plt
import anndata as ad
import shutil
import numpy as np
from pathlib import Path
import spatialdata as sd
import xarray as xr
from spatialdata_io import (
    points_anndata_from_coordinates,
    circles_anndata_from_coordinates,
    table_update_anndata,
)
##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("toy"):
    path /= "toy"
    assert path.exists()

path_write = path / "data.zarr"

##
# image
def draw_circle_on_image(im: np.ndarray, x: int, y: int, r: int):
    yy, xx = np.mgrid[: im.shape[0], : im.shape[1]]
    circle = (xx - x) ** 2 + (yy - y) ** 2
    return np.logical_or(im, circle < r * r).astype(float)


x = np.zeros((200, 100), dtype=float)
x = draw_circle_on_image(im=x, x=30, y=50, r=5)
x = draw_circle_on_image(im=x, x=50, y=20, r=10)
x = draw_circle_on_image(im=x, x=70, y=150, r=15)
x[:2, :] = 1.0
x[-2:, :] = 1.0
x[:, :2] = 1.0
x[:, -2:] = 1.0

plt.figure()
plt.imshow(x, origin="lower")
plt.show()

##
# circles
plt.figure()
x_max = 500
y_max = 800
translation_x = 100
translation_y = 50
scale_factor = 3
a_x, a_y = (30 * scale_factor + translation_x, 50 * scale_factor + translation_y)
b_x, b_y = (50 * scale_factor + translation_x, 20 * scale_factor + translation_y)
c_x, c_y = (70 * scale_factor + translation_x, 150 * scale_factor + translation_y)
points = np.array(
    [
        [translation_x, translation_y],
        [x_max - translation_x, translation_y],
        [x_max - translation_x, y_max - 3 * translation_y],
        (a_x, a_y),
        (b_x, b_y),
        (c_x, c_y),
    ]
)
sizes = [100, 400, 900, 50, 50, 50]


def f(t, begin_end):
    d = 0.25
    segments = np.arange(0, 1.1, d)
    n = np.argmax(t < segments) - 1
    low = segments[n]
    up = segments[n + 1]
    begin, end = begin_end[n]
    s = (t - low) / (up - low)
    return begin * s + end * (1 - s)


def get_points(begin_end):
    tt = np.linspace(0, 1, 50)
    lt = []
    for t in tt:
        lt.append(f(t, begin_end))
    xy = np.array(lt)
    return xy


begin_end = np.array(
    [
        [(1, 1), (x_max - 2, 1)],
        [(x_max - 2, 1), (x_max - 2, y_max - 2)],
        [(x_max - 2, y_max - 2), (1, y_max - 2)],
        [(1, y_max - 2), (1, 1)],
    ]
)
xy_border = get_points(begin_end)
points = np.concatenate((points.astype(float), xy_border))
sizes = sizes + [50] * len(xy_border)

ax = plt.gca()
ax.imshow(
    x,
    origin="lower",
    extent=[
        translation_x,
        x_max - translation_x,
        translation_y,
        y_max - 3 * translation_y,
    ],
)
ax.scatter(points[:, 0], points[:, 1], s=sizes, c="r")
ax.set_aspect("equal")
ax.set(xlim=(0, x_max), ylim=(0, y_max))

# points
begin_end = np.array(
    [
        [points[0, :].tolist(), (a_x, a_y)],
        [(a_x, a_y), (b_x, b_y)],
        [(b_x, b_y), (c_x, c_y)],
        [(c_x, c_y), points[2, :].tolist()],
    ]
)
xy = get_points(begin_end)

ax.scatter(xy[:, 0], xy[:, 1], s=1)
plt.show()
##
a_circles = circles_anndata_from_coordinates(coordinates=points, radii=np.sqrt(np.array(sizes) / np.pi))
a_points = points_anndata_from_coordinates(coordinates=xy)

##
transformation = sd.compose_transformations(
    sd.Scale(scale=np.array([scale_factor, scale_factor])),
    sd.Translation(translation=np.array([translation_x, translation_y])),
)
sdata = sd.SpatialData(
    images={"image": xr.DataArray(x, dims=("y", "x"))},
    points={"points": a_points, "circles": a_circles},
    transformations={
        ("/images/image", "global"): transformation,
    },
)

print(sdata)
##
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")

sdata = sd.SpatialData.read(path_write)
print(sdata)
print('read')
# viewer = load_to_napari_viewer(
#     file_path=output_fpath,
#     groups=["circles/circles_table", "points/points_table"],
#     # groups=["circles/circles_table", "tables/regions_table", "points/points_table"],
# )
# napari.run()
