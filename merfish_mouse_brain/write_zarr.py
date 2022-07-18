##
# the best is to run this .py file as a script, if you want you can also to convert it to a jupyter notebook with
# jupytext --to notebook merfish.py
# When runnning the notebook the working directory should be the project root folder otherwise the data paths will be
# not found (use os.chdir(<project root path>)

##
import anndata as ad
from pathlib import Path
import pandas as pd
import imageio
import json
from ome_zarr.io import parse_url
import zarr
from ome_zarr.io import parse_url
from ome_zarr.writer import write_image, write_labels
import zarr
from anndata.experimental import write_elem
from ngff_tables_prototype.writer import write_table_regions, write_table_points

path = Path().resolve() / "merfish_mouse_brain"
path_read = path / "processed"
path_write = path / "data.zarr"

cells = ad.read(path_read / "cells.h5ad")
transcripts = ad.read(path_read / "single_molecule.h5ad")
transcripts.obs["cell_type"] = pd.Categorical(transcripts.obsm["cell_type"])
image = imageio.imread(path_read / "image.png")
with open(path_read / "image_transform.json", "r") as f:
    image_transform = json.load(f)
with open(path_read / "anatomical.geojson", "r") as f:
    anatomical_region = json.load(f)

store = parse_url(path_write, mode="w").store
root = zarr.group(store=store)

transform_translation = {
    "type": "translation",
    "translation": [
        image_transform["translation_y"],
        image_transform["translation_x"],
    ],
}
transform_scale = {
    "type": "scale",
    "scale": [
        image_transform["scale_factor_y"],
        image_transform["scale_factor_x"],
    ],
}
coordinate_transformations = [[transform_scale, transform_translation]]

write_image(
    image=image,
    group=root,
    coordinate_transformations=coordinate_transformations,
    axes=["y", "x"],
    scaler=None,
)

tables_group = root.create_group(name="tables")
write_table_regions(
    group=tables_group,
    adata=cells,
)
points_group = root.create_group(name="points")
write_table_points(
    group=tables_group,
    adata=transcripts,
)
