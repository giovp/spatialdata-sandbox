##
from spatialdata_io import cosmx
import spatialdata as sd

##
from pathlib import Path
import shutil

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("cosmx_io"):
    path /= "cosmx_io"
    assert path.exists()
path_read = path / "data/data_lung5_rep2"

# variation on the dataset in download.py, from here: https://github.com/scverse/spatialdata-io/issues/272#issue-2836388176
# path_read = "/Users/macbook/ssd/biodata/cosmx/2025 download/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images"
path_write = path / "data.zarr"
##
sdata = cosmx(path_read)
# sdata = cosmx(path_read, transcripts=False)
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
