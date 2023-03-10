##
import os
from pathlib import Path
import subprocess

url = "https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip"

##
# download the data
filename = Path(url).name
os.makedirs("data", exist_ok=True)
command = f"curl -o {'data/' + filename} {url}"
subprocess.run(command, shell=True, check=True)

##
# unzip the data
subprocess.run(f'unzip data/{filename} -d data/xenium', shell=True, check=True)
