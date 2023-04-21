import shutil
from pathlib import Path

from spatialdata_io import merfish


def main():
    sample_path = Path(__file__).parent / "data" / "vizgen_sample"
    zarr_path = Path(__file__).parent / "data.zarr"

    sdata = merfish(sample_path, vpt_outputs=sample_path / "vpt-run")

    if zarr_path.exists():
        shutil.rmtree(zarr_path)
    sdata.write(zarr_path)

    print(f"Saved merfish data:\n{sdata}")


if __name__ == "__main__":
    main()
