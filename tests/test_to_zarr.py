import os
import pathlib

import pytest

path = pathlib.Path(__file__).parent.parent.resolve()


@pytest.mark.skip(reason="large dataset")
def test_cosmx_io():
    os.chdir(path / "cosmx_io")
    import cosmx_io.to_zarr

def test_mcmicro_io():
    os.chdir(path / "mcmicro_io")
    import mcmicro_io.to_zarr

def test_merfish():
    os.chdir(path / "merfish")
    import merfish.to_zarr


def test_mibitof():
    os.chdir(path / "mibitof")
    import mibitof.to_zarr


def test_steinbock_io():
    os.chdir(path / "steinbock_io")
    import steinbock_io.to_zarr

def test_toy():
    os.chdir(path / "toy")
    import toy.to_zarr


def test_visium():
    os.chdir(path / "visium")
    import visium.to_zarr


def test_visium_io():
    os.chdir(path / "visium_io")
    import visium_io.to_zarr


@pytest.mark.skip(reason="large dataset")
def test_xenium_io():
    os.chdir(path / "xenium_io")
    import xenium_rep1_io.to_zarr
