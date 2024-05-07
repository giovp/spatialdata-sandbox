##
import os
from pathlib import Path
import subprocess

urls = [
    "https://cf.10xgenomics.com/samples/xenium/1.5.0/Xenium_V1_hKidney_cancer_section/Xenium_V1_hKidney_cancer_section_outs.zip",
    "https://cf.10xgenomics.com/samples/xenium/1.5.0/Xenium_V1_hKidney_cancer_section/Xenium_V1_hKidney_cancer_section_he_image.ome.tif",
    "https://cf.10xgenomics.com/samples/xenium/1.5.0/Xenium_V1_hKidney_cancer_section/Xenium_V1_hKidney_cancer_section_he_imagealignment.csv",
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
    f"unzip -o data/Xenium_V1_hKidney_cancer_section_outs.zip -d data/",
    shell=True,
    check=True,
)
