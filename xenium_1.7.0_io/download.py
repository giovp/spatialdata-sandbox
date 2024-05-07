##
import os
from pathlib import Path
import subprocess

urls = [
    "https://cf.10xgenomics.com/samples/xenium/1.7.0/Xeniumranger_V1_hSkin_Melanoma_Add_on_FFPE/Xeniumranger_V1_hSkin_Melanoma_Add_on_FFPE_outs.zip",
    "https://cf.10xgenomics.com/samples/xenium/1.7.0/Xeniumranger_V1_hSkin_Melanoma_Add_on_FFPE/Xeniumranger_V1_hSkin_Melanoma_Add_on_FFPE_he_image.ome.tif",
    "https://cf.10xgenomics.com/samples/xenium/1.7.0/Xeniumranger_V1_hSkin_Melanoma_Add_on_FFPE/Xeniumranger_V1_hSkin_Melanoma_Add_on_FFPE_he_imagealignment.csv"
]

##
# download the data
for url in urls:
    filename = Path(url).name
    os.makedirs("data", exist_ok=True)
    command = f"curl -o {'data/' + filename} {url}"
    subprocess.run(command, shell=True, check=True)

##
# unzip the data
subprocess.run(
    f"unzip -o data/Xeniumranger_V1_hSkin_Melanoma_Add_on_FFPE_outs.zip -d data/",
    shell=True,
    check=True,
)
