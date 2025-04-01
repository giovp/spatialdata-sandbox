##
import spatialdata as sd
import spatialdata_plot
import geopandas
import anndata
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from spatialdata.datasets import blobs_annotating_element

##
# create small dataset, write and read to zarr, print it and plot it
sdata = blobs_annotating_element("blobs_polygons")

# full IO roundtrip using the spatialdata APIs (i.e. all batteries included)  
sdata.write("data.zarr", overwrite=True)
sdata = sd.read_zarr("data.zarr")

print(sdata)
# note: the instance_id column is from the AnnData tabular annotation
# below we will show how one can read and write the table
(
    sdata.pl.render_images("blobs_image", cmap="gray", channel=0)
    .pl.render_shapes("blobs_polygons", color="instance_id")
    .pl.show()
)
plt.show()

## save image to png
data = sdata["blobs_image"].data.compute()

data_normalized = (255 * (data - data.min()) / (data.max() - data.min())).astype(
    np.uint8
)
data_normalized = np.transpose(data_normalized, (1, 2, 0))
image = Image.fromarray(data_normalized)

image.save("image.png")

## save polygons to geojson
shapes = sdata["blobs_polygons"]
shapes.to_file("shapes.geojson", driver="GeoJSON")

## save tabular annotations to AnnData
sdata["table"].write_zarr("table.zarr")

## create SpatialData object from scratch, print it and plot it
image = Image.open("image.png")
shapes = geopandas.GeoDataFrame.from_file("shapes.geojson")
table = anndata.read_zarr("table.zarr")

parsed_image = sd.models.Image2DModel.parse(np.array(image), dims=("y", "x", "c"))
parsed_shapes = sd.models.ShapesModel.parse(shapes)
parsed_table = sd.models.TableModel.parse(table)

sdata = sd.SpatialData.init_from_elements(
    {
        "blobs_image": parsed_image,
        "blobs_polygons": parsed_shapes,
        "table": parsed_table,
    }
)
print(sdata)
sdata = sd.read_zarr("data.zarr")
(
    sdata.pl.render_images("blobs_image", cmap="gray", channel=0)
    .pl.render_shapes("blobs_polygons", color="instance_id")
    .pl.show()
)
plt.show()

## alternative 1: write SpatialData object incrementally, then read it, print it and
# plot it
sdata = sd.SpatialData()
sdata.write("incremental.zarr", overwrite=True)

sdata["blobs_image"] = parsed_image
sdata.write_element("blobs_image")

sdata["blobs_polygons"] = parsed_shapes
sdata.write_element("blobs_polygons")

sdata["table"] = parsed_table
sdata.write_element("table")

sdata = sd.read_zarr("incremental.zarr")
print(sdata)
sdata = sd.read_zarr("data.zarr")
(
    sdata.pl.render_images("blobs_image", cmap="gray", channel=0)
    .pl.render_shapes("blobs_polygons", color="instance_id")
    .pl.show()
)
plt.show()

## alternative 2: as above, but using the spatialdata_io APIs as a shorthand and only
# reading it in memory
from spatialdata_io import generic

image = generic("image.png", data_axes=("y", "x", "c"))
shapes = generic("shapes.geojson")
# we already read the table before

sdata = sd.SpatialData.init_from_elements(
    {"blobs_image": image, "blobs_polygons": shapes, "table": table}
)
print(sdata)

## alternative 3: now using the spatialdata_io to even write to disk incrementally
from spatialdata_io import generic_to_zarr
import shutil
import os

if os.path.isdir('generic_to_zarr.zarr'):
    shutil.rmtree("generic_to_zarr.zarr")

generic_to_zarr(
    input="image.png",
    output="generic_to_zarr.zarr",
    name="blobs_image",
    data_axes=("y", "x", "c"),
    coordinate_system="global",
)
generic_to_zarr(
    input="shapes.geojson",
    output="generic_to_zarr.zarr",
    name="blobs_polygons",
    coordinate_system="global",
)
print(sdata)

## alternative 4: the above is equivalent to the CLI usage below
# python -m spatialdata_io generic \
#     --input "image.png" \
#     --output "generic_to_zarr.zarr" \
#     --name "blobs_image" \
#     --data-axes "yxc" \
#     --coordinate-system "global"
#
# python -m spatialdata_io generic \
#     --input "shapes.geojson" \
#     --output "generic_to_zarr.zarr" \
#     --name "blobs_polygons" \
#     --coordinate-system "global"

