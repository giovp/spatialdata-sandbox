from pathlib import Path

from spatialdata_io import merfish


def main():
    sample_path = Path(__file__).parents[1] / "data" / "vizgen_sample"
    zarr_path = Path(__file__).parent / "data.zarr"

    if zarr_path.exists():
        print(f"Zarr folder already existing at {zarr_path}. Exiting.")
        return

    sdata = merfish(sample_path, vpt_outputs=sample_path / "vpt-run")
    sdata.write(zarr_path)
    
    print(f"Saved merfish data:\n{sdata}")


if __name__ == "__main__":
    main()
