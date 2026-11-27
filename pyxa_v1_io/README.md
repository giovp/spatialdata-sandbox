# Stellaromics Pyxa (v1 output files)

The `xsmall` subset of the [Stellaromics/demo](https://huggingface.co/datasets/Stellaromics/demo) dataset on the
Hugging Face Hub (BSD 3-Clause): a 100 × 100 × 100 µm cube (187 cells, ~23k transcripts) with the four Pyxa output
files and a DAPI mosaic (multiscale OME-Zarr), ~8 MB in total. It is the same subset used in the spatialdata-io CI
tests for the experimental `pyxa` reader.

`download.py` fetches it into `data/`, and `to_zarr.py` converts it with `spatialdata_io.experimental.pyxa()`. The
full demo region (`small/`, ~290 MB) has the same layout.
