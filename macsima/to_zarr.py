##
from spatialdata_io import macsima
import spatialdata as sd
##
from pathlib import Path
import shutil

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("macsima_io"):
    path /= "macsima_io"
    assert path.exists()
path_read = path / "data/Lung_adc_demo"
path_write = path / "data.zarr"
##
sdata = macsima(path_read, subset=1000)
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
