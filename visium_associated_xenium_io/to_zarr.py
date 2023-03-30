##
from spatialdata_io import visium
import spatialdata as sd
import shutil

##
from pathlib import Path


##
def main():
    ##
    path = Path().resolve()
    # luca's workaround for pycharm
    if not str(path).endswith("visium_associated_xenium_io"):
        path /= "visium_associated_xenium_io"
        assert path.exists()
    path_read = path / "data"
    path_write = path / "data.zarr"
    ##
    print("reading the data... ", end="")
    sdata_visium = visium(str(path_read))
    print("done")

    print("writing the data... ", end="")
    if path_write.exists():
        shutil.rmtree(path_write)
    sdata_visium.write(str(path_write))
    print("done")
    ##
    print(f'view with "python -m napari_spatialdata view data.zarr"')
    print("read")
    sdata = sd.SpatialData.read(path_write)
    print(sdata)
    ##


if __name__ == "__main__":
    main()
