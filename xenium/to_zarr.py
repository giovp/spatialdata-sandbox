##
from spatialdata_io import convert_xenium_to_ngff, read_visium
import spatialdata as sd

##
from pathlib import Path
import shutil

##

# since we use multiprocessing we need if __name__ == "__main__", otherwise we get this error https://stackoverflow.com/questions/55057957/an-attempt-has-been-made-to-start-a-new-process-before-the-current-process-has-f
# we also enclose the code in a function, so that we can call it from the tests
def main():
    ##
    path = Path().resolve()
    # luca's workaround for pycharm
    if not str(path).endswith("xenium"):
        path /= "xenium"
        assert path.exists()
    path_read = path / "data"
    path_write = path / "data.zarr"
    ##
    convert_xenium_to_ngff(
        str(path_read),
        str(path_write),
        # skip_nucleus_boundaries=True,
        # skip_cell_boundaries=True,
        # skip_points=True,
        skip_table_and_shapes=True,
        skip_image_morphology=True,
        skip_image_morphology_mip=True,
        skip_image_morphology_focus=True,
    )
    print("done")
    ##
    print(f'view with "python -m spatialdata view data.zarr"')
    print("read")
    sdata = sd.SpatialData.read("./data.zarr/")
    print(sdata)
    ##


if __name__ == "__main__":
    main()
