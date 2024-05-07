##
import os
from pathlib import Path
import subprocess

urls = [
    "https://cf.10xgenomics.com/samples/xenium/1.9.0/Xenium_V1_hHeart_nondiseased_section_FFPE/Xenium_V1_hHeart_nondiseased_section_FFPE_he_image.ome.tif",
    "https://cf.10xgenomics.com/samples/xenium/1.9.0/Xenium_V1_hHeart_nondiseased_section_FFPE/Xenium_V1_hHeart_nondiseased_section_FFPE_he_imagealignment.csv",
    "https://cf.10xgenomics.com/samples/xenium/1.9.0/Xenium_V1_hHeart_nondiseased_section_FFPE/Xenium_V1_hHeart_nondiseased_section_FFPE_outs.zip",
]

##
# download the data
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("xenium_1.9.0_human_heart_io"):
    path /= "xenium_1.9.0_human_heart_io"
    assert path.exists()

path = path / "data"

##
# download the data
for url in urls:
    filename = Path(url).name
    os.makedirs(path, exist_ok=True)
    command = f"curl -o {path/filename} {url}"
    subprocess.run(command, shell=True, check=True)

# ##
# unzip the data
subprocess.run(
    f"unzip -o {path}/Xenium_V1_hHeart_nondiseased_section_FFPE_outs.zip -d {path}/",
    shell=True,
    check=True,
)