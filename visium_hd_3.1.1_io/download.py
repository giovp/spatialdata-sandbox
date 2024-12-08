#!/usr/bin/env python3
##
import os
import subprocess
from pathlib import Path
import shutil

import sys

sys.path.insert(1, os.path.join(sys.path[0], Path(__file__).parent.parent.resolve()))

from utils import download

##

# from https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-human-lung-cancer-fixed-frozen
# Visium HD Spatial Gene Expression Library, Human Lung Cancer (Fixed Frozen)

urls = [
    "https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Lung_Cancer_Fixed_Frozen/Visium_HD_Human_Lung_Cancer_Fixed_Frozen_cloupe_008um.cloupe",
    "https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Lung_Cancer_Fixed_Frozen/Visium_HD_Human_Lung_Cancer_Fixed_Frozen_feature_slice.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Lung_Cancer_Fixed_Frozen/Visium_HD_Human_Lung_Cancer_Fixed_Frozen_spatial.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Lung_Cancer_Fixed_Frozen/Visium_HD_Human_Lung_Cancer_Fixed_Frozen_binned_outputs.tar.gz"
]


##
os.makedirs("data", exist_ok=True)
for url in urls:
    print(url)
    name = Path(url).name
    download(url, os.path.join("data", name), name)

files = [
    "Visium_HD_Human_Lung_Cancer_Fixed_Frozen_spatial.tar.gz",
    "Visium_HD_Human_Lung_Cancer_Fixed_Frozen_binned_outputs.tar.gz",
]
for file in files:
    subprocess.run(f"tar -xzf data/{file} -C data", shell=True, check=True)
    subprocess.run(f"rm data/{file}", shell=True, check=True)

for dir in list(Path("data/binned_outputs").glob("*")):
    if not(Path("data") / dir.name).exists():
        shutil.move(dir,"data")
subprocess.run(f"rm -r data/binned_outputs", shell=True, check=True)
