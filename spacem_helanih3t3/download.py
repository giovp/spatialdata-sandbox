#!/usr/bin/env python3
import os
import subprocess

URL = "TODO",

os.chdir(os.path.dirname(__file__))
command = f"curl {url} --output 'data.zip'"
subprocess.run(command, shell=True, check=True)
subprocess.run("unzip -f -o data.zip", shell=True, check=True)
subprocess.run("mv spatialdata.zarr data.zarr", shell=True, check=True)
