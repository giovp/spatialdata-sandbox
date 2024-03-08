##
from spatialdata_io import visium
import spatialdata as sd
##
from pathlib import Path
import shutil

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("visium_1.0.0_io"):
    path /= "visium_1.0.0_io"
    assert path.exists()
path_read = path / "data"
path_write = path / "data.zarr"
##
sdata = visium(path_read, fullres_image_file="V1_Adult_Mouse_Brain_image.tif")
##
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")
##
print(f'view with "python -m napari_spatialdata view data.zarr"')
sdata = sd.SpatialData.read(path_write)
print(sdata)
print("read")

# from napari_spatialdata import Interactive
# Interactive(sdata)
