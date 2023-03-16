##
import shutil
from pathlib import Path
import anndata as ad
import imageio.v3 as iio
from spatialdata import SpatialData
from spatialdata.transformations import Identity
from spatialdata.models import Labels2DModel, Image2DModel, TableModel

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("mibitof"):
    path /= "mibitof"
    assert path.exists()
path_read = path / "data"
path_write = path / "data.zarr"

##
libraries = ["point8", "point16", "point23"]

table_list = []
for lib in libraries:
    table = ad.read(path_read / f"{lib}_table.h5ad")
    table.obs["library_id"] = f"{lib}_labels"
    table_list.append(table)

table = ad.concat(
    table_list,
    keys=libraries,
)
table = TableModel.parse(
    adata=table,
    region=[f"{lib}_labels" for lib in libraries],
    region_key="library_id",
    instance_key="cell_id",
)

##
labels = {
    f"{lib}_labels": Labels2DModel.parse(
        iio.imread(path_read / f"{lib}_labels.png"),
        dims=("y", "x"),
        transformations={lib: Identity()},
    )
    for lib in libraries
}
images = {
    f"{lib}_image": Image2DModel.parse(
        iio.imread(path_read / f"{lib}_image.png"),
        dims=("y", "x", "c"),
        transformations={lib: Identity()},
    )
    for lib in libraries
}

##

table.obs['donor'] = table.obs['donor'].astype('category')
table.obs['batch'] = table.obs['batch'].astype('category')

sdata = SpatialData(
    table=table,
    labels=labels,
    images=images,
)
print(sdata)

##
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")
print(f'view with "python -m spatialdata view data.zarr"')
##
sdata = SpatialData.read(path_write)
print(sdata)
print("read")
