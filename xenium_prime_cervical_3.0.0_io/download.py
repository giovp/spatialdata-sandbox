# Source: 10x Genomics Xenium Prime Cervical Cancer FFPE dataset
# https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-cervical-cancer
# License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
##
import os
from pathlib import Path
import subprocess

urls = [
    "https://s3-us-west-2.amazonaws.com/10x.files/samples/xenium/3.0.0/Xenium_Prime_Cervical_Cancer_FFPE/Xenium_Prime_Cervical_Cancer_FFPE_outs.zip",
    "https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_Cervical_Cancer_FFPE/Xenium_Prime_Cervical_Cancer_FFPE_he_image.ome.tif",
    "https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_Cervical_Cancer_FFPE/Xenium_Prime_Cervical_Cancer_FFPE_he_imagealignment.csv",
    "https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_Cervical_Cancer_FFPE/Xenium_Prime_Cervical_Cancer_FFPE_gene_panel.json",
]

##
# download the data
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("xenium_prime_cervical_3.0.0_io"):
    path /= "xenium_prime_cervical_3.0.0_io"
    assert path.exists()

path = path / "data"

##
# download the data
for url in urls:
    filename = Path(url).name
    os.makedirs(path, exist_ok=True)
    command = f"curl -o {path/filename} {url}"
    subprocess.run(command, shell=True, check=True)

##
# prepare files for spatialdata-io ingestion
outs_path = path / "Xenium_Prime_Cervical_Cancer_FFPE_outs"
subprocess.run(
    f"unzip -o {path}/Xenium_Prime_Cervical_Cancer_FFPE_outs.zip -d {outs_path}",
    shell=True,
    check=True,
)
subprocess.run(
    f"mv {path}/Xenium_Prime_Cervical_Cancer_FFPE_he_image.ome.tif {outs_path}/",
    shell=True,
    check=True,
)
subprocess.run(
    f"mv {path}/Xenium_Prime_Cervical_Cancer_FFPE_he_imagealignment.csv {outs_path}/",
    shell=True,
    check=True,
)
subprocess.run(
    f"mv {path}/Xenium_Prime_Cervical_Cancer_FFPE_gene_panel.json {outs_path}/",
    shell=True,
    check=True,
)
