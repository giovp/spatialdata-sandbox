import anndata as ad
from pathlib import Path
import pandas as pd
import numpy as np
import imageio
import json
from ome_zarr.io import parse_url
import zarr
from ome_zarr.io import parse_url
from ome_zarr.writer import write_image, write_labels
import zarr
from anndata.experimental import write_elem
from ngff_tables_prototype.writer import write_table_regions, write_table_points

path = Path().resolve()
path_read = path / "data" / "data_lung5_rep2"
path_write = path / "data.zarr"

# read table
table = ad.read(path_read / "nanostring_lung5_rep2.h5ad")

# read points
points_df = pd.read_parquet(path_read / "nanostring_lung5_rep2_single_genes.parquet")
points_df.reset_index(inplace=True, drop=False)
points_df["fov"] = pd.Categorical(points_df["fov"].astype(str))
points_df.reset_index(inplace=True, drop=False)
cols = ["y_local_px", "x_local_px"]
points = ad.AnnData(
    points_df[["y_local_px", "x_local_px"]].to_numpy(),
    obs=points_df[points_df.columns.difference(cols)],
    dtype=np.float_,
)

store = parse_url(path_write, mode="w").store
root = zarr.group(store=store)

for fov in table.obs.fov.cat.categories:
    sfov = fov.zfill(3)

    image = imageio.imread(path_read / "CellComposite" / f"CellComposite_F{sfov}.jpg")
    labels = imageio.imread(path_read / "CellLabels" / f"CellLabels_F{sfov}.tif")
    image_group = root.create_group(name=fov)

    write_image(
        image=image,
        group=image_group,
        axes=["c", "y", "x"],
        scaler=None,
    )
    print(f"Written image `{fov}` .")

    write_labels(
        labels=labels,
        group=image_group,
        name=fov,
        axes=["y", "x"],
        scaler=None,
    )
    print(f"Written labels `{fov}` .")

    write_table_points(
        group=image_group,
        adata=points[points.obs.fov == fov].copy(),
    )
    print(f"Written points `{fov}` .")


tables_group = root.create_group(name="tables")
write_table_regions(
    group=tables_group,
    adata=table,
    region=table.obs.fov.cat.categories.tolist(),
    region_key="fov",
)
