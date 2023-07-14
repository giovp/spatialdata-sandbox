import shutil
from pathlib import Path

from spatialdata_io import merscope


def main():
    sample_path = Path(__file__).parent / "data" / "vizgen_sample"
    zarr_path = Path(__file__).parent / "data.zarr"

    sdata = merscope(sample_path, vpt_outputs=sample_path / "vpt-run")

    if zarr_path.exists():
        shutil.rmtree(zarr_path)
    sdata.write(zarr_path)

    print(f"Saved merscope data:\n{sdata}")


def viz():
    from spatialdata import read_zarr
    from napari_spatialdata import Interactive

    sdata = read_zarr(Path(__file__).parent / "data.zarr")
    Interactive(sdata).run()


if __name__ == "__main__":
    main()
    # viz()
