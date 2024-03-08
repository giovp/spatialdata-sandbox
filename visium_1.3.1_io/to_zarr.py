# The (only) public dataset with spaceranger version 1.3.1 has a structure which is very different from the usual Visium
# dataset. This version is not supported.
# Differences:
# - no tissue positions file
# - the spatial folder has many sufolders, one for each sample
# - each subfolder contains the scalefactors, and the hires and lowres images, but not the tissue positions

# ##
# from spatialdata_io import visium
# import spatialdata as sd
# ##
# from pathlib import Path
# import shutil
#
# ##
# path = Path().resolve()
# # luca's workaround for pycharm
# if not str(path).endswith("visium_1.3.1_io"):
#     path /= "visium_1.3.1_io"
#     assert path.exists()
# path_read = path / "data"
# path_write = path / "data.zarr"
# ##
# sdata = visium(path_read, fullres_image_file=None)
# ##
# if path_write.exists():
#     shutil.rmtree(path_write)
# sdata.write(path_write)
# print("done")
# ##
# print(f'view with "python -m napari_spatialdata view data.zarr"')
# sdata = sd.SpatialData.read(path_write)
# print(sdata)
# print("read")
#
# from napari_spatialdata import Interactive
# Interactive(sdata)
