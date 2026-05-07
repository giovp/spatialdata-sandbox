##
from spatialdata_io import xenium
import spatialdata as sd

##
from pathlib import Path
import shutil

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("atera_breast_cancer_io"):
    path /= "atera_breast_cancer_io"
    assert path.exists()

path_read = path / "data" / "WTA_Preview_FFPE_Breast_Cancer_outs"
path_write = path / "data.zarr"

##
print("parsing the data... ", end="")
sdata = xenium(
    path=str(path_read),
    n_jobs=8,
    cell_boundaries=True,
    nucleus_boundaries=True,
    morphology_focus=True,
    cells_as_circles=False,
)
print("done")

##
print("writing the data... ", end="")
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")

##
sdata = sd.SpatialData.read("./data.zarr/")
print(sdata)
##
