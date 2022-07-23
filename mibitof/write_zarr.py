import anndata as ad
from pathlib import Path
import imageio
from ome_zarr.io import parse_url
import zarr
from ome_zarr.io import parse_url
from ome_zarr.writer import write_image, write_labels
import zarr
from ngff_tables_prototype.writer import write_table_regions

path = Path().resolve()
path_read = path / "data"
path_write = path / "data.zarr"

libraries = ["point8", "point16", "point23"]

table = ad.concat(
    [ad.read(path_read / f"{lib}_table.h5ad") for lib in libraries],
    label="region_key",
    keys=libraries,
)

store = parse_url(path_write, mode="w").store
root = zarr.group(store=store)

for lib in libraries:
    image = imageio.imread(path_read / f"{lib}_image.png")
    labels = imageio.imread(path_read / f"{lib}_labels.png")

    image_group = root.create_group(name=lib)

    write_image(
        image=image,
        group=image_group,
        axes=["c", "y", "x"],
        scaler=None,
    )
    print(f"Written image `{lib}` .")

    write_labels(
        labels=labels,
        group=image_group,
        name=lib,
        axes=["y", "x"],
        scaler=None,
    )
    print(f"Written labels `{lib}` .")

tables_group = root.create_group(name="tables")
write_table_regions(
    group=tables_group,
    adata=table,
)
