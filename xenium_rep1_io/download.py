##
import os
from pathlib import Path
import subprocess

url = "https://s3-us-west-2.amazonaws.com/10x.files/samples/xenium/1.0.2/Xenium_V1_FFPE_Human_Breast_IDC_Big_1/Xenium_V1_FFPE_Human_Breast_IDC_Big_1_outs.zip"

##
# download the data
filename = Path(url).name
os.makedirs("data", exist_ok=True)
command = f"curl -o {'data/' + filename} {url}"
subprocess.run(command, shell=True, check=True)

##
# unzip the data
subprocess.run(f'unzip data/{filename} -d data/xenium', shell=True, check=True)
