import os
from pathlib import Path
from tqdm import tqdm

urls = [
    "https://github.com/labsyspharm/mcmicro/suites/9547753268/artifacts/454126381",
]

os.makedirs("data", exist_ok=True)
os.makedirs("data/steinbock", exist_ok=True)
for url in tqdm(urls, desc="downloading"):
    command = f"curl {url} --output 'data/steinbock/{Path(url).name}'"
    os.system(command)

# os.chdir("data/xenium")
# os.system("unzip Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip")
