import os


def test_merfish():
    os.chdir("../merfish")
    import merfish.to_zarr


def test_mibitof():
    os.chdir("../mibitof")
    import mibitof.to_zarr


def test_nanostring_cosmx():
    os.chdir("../nanostring_cosmx")
    import nanostring_cosmx.to_zarr


def test_toy():
    os.chdir("../toy")
    import toy.to_zarr


def test_visium():
    os.chdir("../visium")
    import visium.to_zarr


def test_visium2():
    os.chdir("../visium2")
    import visium2.to_zarr


def test_xenium():
    os.chdir("../xenium")
    from xenium.to_zarr import main

    main()
