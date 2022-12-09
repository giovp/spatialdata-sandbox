import os
import pathlib

path = pathlib.Path(__file__).parent.parent.resolve()

def test_merfish():
    os.chdir(path / "merfish")
    import merfish.to_zarr


def test_mibitof():
    os.chdir(path / "mibitof")
    import mibitof.to_zarr


def test_nanostring_cosmx():
    os.chdir(path / "nanostring_cosmx")
    import nanostring_cosmx.to_zarr


def test_toy():
    os.chdir(path / "toy")
    import toy.to_zarr


def test_visium():
    os.chdir(path / "visium")
    import visium.to_zarr


def test_visium2():
    os.chdir(path / "visium2")
    import visium2.to_zarr


def test_xenium():
    os.chdir(path / "xenium")
    from xenium.to_zarr import main

    main()
