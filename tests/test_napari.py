import pytest
from spatialdata import SpatialData
from napari_spatialdata import Interactive


def test_merfish():
    sdata = SpatialData.read("merfish/data.zarr")
    Interactive(sdata, headless=True)


def test_mibitof():
    sdata = SpatialData.read("mibitof/data.zarr")
    Interactive(sdata, headless=True)


def test_cosmx_io():
    sdata = SpatialData.read("cosmx_io/data.zarr")
    Interactive(sdata, headless=True)


def test_toy():
    sdata = SpatialData.read("toy/data.zarr")
    Interactive(sdata, headless=True)


def test_visium():
    sdata = SpatialData.read("visium/data.zarr")
    Interactive(sdata, headless=True)


def test_visium_io():
    sdata = SpatialData.read("visium_2.0.0_1_io/data.zarr")
    Interactive(sdata, headless=True)


@pytest.mark.skip(reason="Skipping because skipping the corrsponding to_zarr()")
def test_xenium_io():
    sdata = SpatialData.read("xenium_io/data.zarr")
    Interactive(sdata, headless=True)


if __name__ == "__main__":
    test_merfish()
