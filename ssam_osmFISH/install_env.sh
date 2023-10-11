#!/bin/bash

source ~/.bash_profile
echo "-------- remove old env --------"
conda deactivate
set -e
conda remove -n ome --all -y

echo "-------- create new env --------"
mamba create -n ome python==3.10 -y
conda activate ome

# echo "-------- install spatialdata-io --------"
# cd spatialdata-io
# pip install -e .

# echo "-------- install spatialdata-plot --------"
# cd ../spatialdata-plot
# pip install -e .

# echo "-------- install napari-spatialdata --------"
# cd ../napari-spatialdata
# pip install numpy magicgui qtpy numba anndata scikit-image squidpy loguru jupyter-black
# pip install PyQt6
# pip install -e . --no-deps

# echo "-------- install spatialdata --------"
# cd ../spatialdata
# pip install -e ".[all]"
# pip install torch torchvision black pre-commit tqdm scanpy pyarrow imagecodecs loguru pytest psutil napari-ome-zarr myst_nb jupyterlab napari-matplotlib apache-airflow readfcs

# echo "-------- fix a problem with M1 chip: https://github.com/spatial-image/multiscale-spatial-image/issues/77 --------"
# pip uninstall numcodecs -y
# mamba install numcodecs -y
# pip uninstall vispy -y
# mamba install vispy -y
# # other conflicting packages
# mamba install -c conda-forge xarray xarray-spatial xarray-datatree typing_extensions scikit-image dask scipy -y

# # echo "-------- fix the pygeos/shapely incompatibility (I don't think this works) --------"
# pip uninstall geopandas shapely pygeos -y
# mamba install -c conda-forge shapely pygeos geopandas -y

mamba install -c conda-forge -c anaconda "anndata>=0.9.1" jupyterlab bump2version dask-image geopandas h5py imagecodecs "ipython>=8.6.0" joblib loguru magicgui matplotlib myst-nb napari networkx numba numpy pre-commit pyarrow pygeos pytest pytest-cov qtpy rich scanpy scikit-image scikit-learn "shapely>=2.0.1" sphinx-autodoc-typehints "sphinx-book-theme>=1.0.0" sphinx-copybutton sphinx-design "sphinx>=4.5" sphinx_rtd_theme "sphinxcontrib-bibtex>=1.0.0" squidpy tqdm "typing_extensions>=4.0.0" xarray xarray-schema "xarray-spatial>=0.3.5" zarr pytorch::pytorch torchvision torchaudio -c pytorch -y

# torch
# "typing_extensions>=4.0.0,<4.6.0"

pip install "ome_zarr>=0.7.0" "multiscale_spatial_image>=0.11.2" "matplotlib_scalebar" "readfcs" "spatial_image>=0.3.0" "PyQt6" "jupyter-black" "napari-matplotlib"

echo "-------- installing spatialdata and related libraries --------"
cd spatialdata
pip install -e . --no-deps

cd ../spatialdata-io
pip install -e . --no-deps

cd ../spatialdata-plot
pip install -e . --no-deps

cd ../napari-spatialdata
pip install -e . --no-deps

echo "-------- test --------"
cd ../spatialdata-sandbox/mibitof
python -m spatialdata peek data.zarr
python -m napari_spatialdata view data.zarr
