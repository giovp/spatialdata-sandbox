##
import os
from pathlib import Path
import subprocess

# from https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-skin

urls = [
    "curl -O https://...zip"
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
    f"unzip -o data/....zip -d data/",
    shell=True,
    check=True,
)
