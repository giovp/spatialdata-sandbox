##
import os
from pathlib import Path
from tqdm import tqdm

url = "https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip"

##
# download the data
filename = Path(url).name
os.makedirs("data", exist_ok=True)
command = f"curl -O {url} --output {'data/' + filename}"
os.system(command)

##
# unzip the data
os.system(f'unzip data/{filename} -d data/xenium')
