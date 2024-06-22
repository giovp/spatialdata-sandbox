##
from spatialdata_io import visium_hd
import spatialdata as sd

##
from pathlib import Path
import shutil

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("visium_hd_3.0.0_io"):
    path /= "visium_hd_3.0.0_io"
    assert path.exists()
path_read = path / "data"
path_write = path / "data.zarr"

##
sdata = visium_hd(
    path_read,
    load_all_images=True,
    fullres_image_file="Visium_HD_Mouse_Small_Intestine_tissue_image.btf",
    bin_size=[2, 16]
)

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
