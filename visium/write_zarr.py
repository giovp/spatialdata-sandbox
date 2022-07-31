from typing import Mapping, Any, List, Dict, Union
import anndata as ad
from pathlib import Path
import pandas as pd
import numpy as np
import imageio
import re
from ome_zarr.io import parse_url
import zarr
from ome_zarr.io import parse_url
from ome_zarr.writer import write_image, write_labels
import zarr
from anndata.experimental import write_elem
from ngff_tables_prototype.writer import write_table_regions, write_table_points
from ome_zarr.format import CurrentFormat
from anndata import AnnData
import numpy as np
import pandas as pd
from ome_zarr.io import parse_url
from ome_zarr.writer import write_image, write_labels
import zarr
from anndata.experimental import write_elem
from ome_zarr.writer import _get_valid_axes, _validate_datasets
from ome_zarr.types import JSONDict
import os

# taken from prototype
def write_table_shape(
    group: zarr.Group,
    adata: AnnData,
    table_group_name: str = "shape_table",
    group_type: str = "ngff:shape_region",
    shape_type: str = "circle",
    axes: List[str] = ["y", "x"],
    scale: float = 1.0,
    shape_params: Mapping = {"radius": 1},
    **metadata: Union[str, JSONDict, List[JSONDict]],
):
    write_elem(group, table_group_name, adata)

    fmt = CurrentFormat()
    axes = _get_valid_axes(adata.X.ndim, axes, fmt)
    datasets: List[dict] = []
    datasets.append({"path": str(0)})

    transform_scale = {
        "type": "scale",
        "scale": [
            scale,
            scale,
        ],
    }
    coordinate_transformations = [[transform_scale]]

    fmt.validate_coordinate_transformations(adata.X.ndim, 1, coordinate_transformations)

    if coordinate_transformations is not None:
        for dataset, transform in zip(datasets, coordinate_transformations):
            dataset["coordinateTransformations"] = transform

    # datasets["coordinateTransformations"] = coordinate_transformations

    multiscales = [
        dict(
            version=fmt.version,
            datasets=_validate_datasets(datasets, adata.X.ndim, fmt),
            **metadata,
        )
    ]
    if axes is not None:
        multiscales[0]["axes"] = axes

    table_group = group[table_group_name]
    table_group.attrs["multiscales"] = multiscales
    table_group.attrs["@type"] = group_type
    table_group.attrs["shape_type"] = shape_type
    table_group.attrs["shape_params"] = shape_params


path = Path().resolve()
path_read = path / "data"
path_write = path / "data.zarr"

libraries = [re.sub(".h5ad", "", i) for i in os.listdir(path_read / "tables")]

store = parse_url(path_write, mode="w").store
root = zarr.group(store=store)

table_list = []
for lib in libraries:
    table = ad.read(path_read / "tables" / "ST8059051.h5ad")
    lib_key = list(table.uns["spatial"].keys())[0]
    img = table.uns["spatial"][lib_key]["images"]["hires"]

    image_group = root.create_group(name=lib)

    write_image(
        image=img,
        group=image_group,
        axes=["c", "y", "x"],
        scaler=None,
    )
    print(f"Written image `{lib}` .")

    shape_region = ad.AnnData(table.obsm["spatial"][:, [1, 0]])

    write_table_shape(
        group=image_group,
        adata=shape_region,
        table_group_name="circles",
        axes=["y", "x"],
        scale=table.uns["spatial"][lib_key]["scalefactors"]["tissue_hires_scalef"],
        shape_type="circle",
        shape_params={
            "radius": table.uns["spatial"][lib_key]["scalefactors"][
                "spot_diameter_fullres"
            ]
            / 2
        },
    )
    table.uns.pop("spatial")
    table.var_names_make_unique()
    table_list.append(table)

    print(f"Written shape region `{lib}` .")

table = ad.concat(
    table_list,
    label="library",
    keys=libraries,
)

tables_group = root.create_group(name="tables")
write_table_regions(
    group=tables_group,
    adata=table,
    region=table.obs.library.cat.categories.tolist(),
    region_key="library",
)
