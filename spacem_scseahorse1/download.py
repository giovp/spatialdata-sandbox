#!/usr/bin/env python3
import os
import subprocess

URL = "https://s3.embl.de/spatialdata/raw_data/20220121_ScSeahorse1.small.zip"

os.chdir(os.path.dirname(__file__))
command = f"curl {URL} --output 'data.zip'"
subprocess.run(command, shell=True, check=True)
os.makedirs('data', exist_ok=True)
os.rename('data.zip', 'data/data.zip')
subprocess.run("unzip -o data/data.zip -d data/", shell=True, check=True)
subprocess.run("mv data/spatialdata.zarr data/data.zarr", shell=True, check=True)