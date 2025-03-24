##
import os
from pathlib import Path
import subprocess

# Singular Genomics sample data

urls = [
    "https://singular-public-repo.s3.us-west-1.amazonaws.com/g4x_tutorial_dataset.zip"
]

# download the data
for url in urls:
    filename = Path(url).name
    os.makedirs("data", exist_ok=True)
    command = f"curl -o {'data/' + filename} {url}"
    subprocess.run(command, shell=True, check=True)

##
# unzip the data
subprocess.run(
    "unzip -o data/g4x_tutorial_dataset.zip -d data/",
    shell=True,
    check=True,
)
