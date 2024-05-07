##
from pathlib import Path
from spatialdata_io.experimental import iss
import shutil
import spatialdata as sd

##
path = Path().resolve()
if not str(path).endswith("iss_io"):
    path /= "iss_io"
    assert path.exists()

dataset = Path("data")

path_read = path / dataset
path_write = path / "data.zarr"

WRITE = True
# WRITE = False

if WRITE:
    ##
    sdata = iss(
        path_read,
        raw_relative_path="iss-demo-raw.tif",
        labels_relative_path="iss-demo-label.tif",
        h5ad_relative_path="iss-demo-anndata.h5ad",
        instance_key='cell_id',
    )

    ##
    print("writing the data... ", end="")
    if path_write.exists():
        shutil.rmtree(path_write)
    sdata.write(path_write)
    print("done")

    ##
print(f'view with "python -m napari_spatialdata view data.zarr"')
sdata = sd.SpatialData.read(path_write)
print("read")
print(sdata)

##
SHOW = False
SHOW = True
if SHOW:
    from napari_spatialdata import Interactive

    interactive = Interactive(sdata, headless=True)
    interactive.add_element(element="region_raw_image", element_coordinate_system="global")
    interactive.add_element(
        element="region_labels_image",
        element_coordinate_system="global",
        view_element_system=True,
    )
    interactive.run()
