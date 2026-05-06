##
import os
from pathlib import Path
import subprocess

# from https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-cervical-cancer

urls = [
    "https://s3-us-west-2.amazonaws.com/10x.files/samples/atera/dev/WTA_Preview_FFPE_Cervical_Cancer/WTA_Preview_FFPE_Cervical_Cancer_outs.zip",
    "https://s3-us-west-2.amazonaws.com/10x.files/samples/atera/dev/WTA_Preview_FFPE_Cervical_Cancer/WTA_Preview_FFPE_Cervical_Cancer_xe_outs.zip",
    "https://cf.10xgenomics.com/samples/atera/dev/WTA_Preview_FFPE_Cervical_Cancer/WTA_Preview_FFPE_Cervical_Cancer_gene_groups.csv",
    "https://cf.10xgenomics.com/samples/atera/dev/WTA_Preview_FFPE_Cervical_Cancer/WTA_Preview_FFPE_Cervical_Cancer_cell_groups.csv",
    "https://cf.10xgenomics.com/samples/atera/dev/WTA_Preview_FFPE_Cervical_Cancer/WTA_Preview_FFPE_Cervical_Cancer_he_image.ome.tif",
    "https://cf.10xgenomics.com/samples/atera/dev/WTA_Preview_FFPE_Cervical_Cancer/WTA_Preview_FFPE_Cervical_Cancer_he_alignment.csv",
    "https://cf.10xgenomics.com/samples/atera/dev/WTA_Preview_FFPE_Cervical_Cancer/WTA_Preview_FFPE_Cervical_Cancer_keypoints.csv",
]
print('current working directory', os.getcwd())
##
# download the data
os.makedirs("data", exist_ok=True)
for url in urls:
    filename = Path(url).name
    command = f"curl -L -o {'data/' + filename} {url}"
    subprocess.run(command, shell=True, check=True)

##
# unzip the data
subprocess.run(
    f"unzip -o data/WTA_Preview_FFPE_Cervical_Cancer_outs.zip -d data/",
    shell=True,
    check=True,
)
subprocess.run(
    f"unzip -o data/WTA_Preview_FFPE_Cervical_Cancer_xe_outs.zip -d data/",
    shell=True,
    check=True,
)
