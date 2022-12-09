from spatialdata import SpatialData
from napari_spatialdata import Interactive


def test_merfish():
    sdata = SpatialData.read("merfish/data.zarr")
    Interactive(sdata, headless=True)


def test_mibitof():
    sdata = SpatialData.read("mibitof/data.zarr")
    Interactive(sdata, headless=True)


def test_nanostring_cosmx():
    sdata = SpatialData.read("nanostring_cosmx/data.zarr")
    Interactive(sdata, headless=True)


def test_toy():
    sdata = SpatialData.read("toy/data.zarr")
    Interactive(sdata, headless=True)


def test_visium():
    sdata = SpatialData.read("visium/data.zarr")
    Interactive(sdata, headless=True)


def test_visium2():
    sdata = SpatialData.read("visium2/data.zarr")
    Interactive(sdata, headless=True)


def test_xenium():
    sdata = SpatialData.read("xenium/data.zarr")
    Interactive(sdata, headless=True)
