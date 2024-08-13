#!/usr/bin/env python3
import os
import subprocess

URL = "https://s3.embl.de/spatialdata/raw_data/20220121_ScSeahorse1.small.zip"

os.chdir(os.path.dirname(__file__))
command = f"curl {URL} --output 'data.zip'"
subprocess.run(command, shell=True, check=True)
subprocess.run("unzip -o data.zip", shell=True, check=True)
subprocess.run("mv spatialdata.zarr data.zarr", shell=True, check=True)
