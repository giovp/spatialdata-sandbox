##
import os
from pathlib import Path
import subprocess

# from https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-skin

urls = [
    "https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_Human_Skin_FFPE/Xenium_Prime_Human_Skin_FFPE_outs.zip"
]
print('current working directory', os.getcwd())
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
    f"unzip -o data/Xenium_Prime_Human_Skin_FFPE_outs.zip -d data/",
    shell=True,
    check=True,
)
