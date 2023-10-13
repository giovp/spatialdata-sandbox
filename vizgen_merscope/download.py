import subprocess
from pathlib import Path


def main():
    URL = Path("https://zenodo.org/record/7852005/files/vizgen_sample.zip")

    zip_path = Path(__file__).parent / "data" / URL.name
    folder_path = Path(__file__).parent / "data" / URL.stem
    zip_path.parent.mkdir(exist_ok=True)

    command = f"curl {URL} --output '{zip_path}'"
    subprocess.run(command, shell=True, check=True)
    subprocess.run(f"unzip -f {zip_path} -d {folder_path}", shell=True, check=True)


if __name__ == "__main__":
    main()
