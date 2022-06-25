```
wget https://cell2location.cog.sanger.ac.uk/tutorial/mouse_brain_visium_wo_cloupe_data.zip
unzip mouse_brain_visium_wo_cloupe_data.zip
```

This will download mouse brain Visium data from the [cell2location paper](https://doi.org/10.1038/s41587-021-01139-4).
There are five Visium slides, located in `rawdata/ST8059048`, `rawdata/ST8059049`, `rawdata/ST8059050`, `rawdata/ST8059051`, `rawdata/ST8059052`.
Use scanpy's `read_visium` function to read them, e.g. `scanpy.read_visium("rawdata/ST8059048")`.
This will return an AnnData object with gene counts. Spatial coordinates are stored in `.obsm["spatial"]`, images are stored in `.uns["spatial"]['spaceranger100_count_30458_ST8059048_mm10-3_0_0_premrna']['images']`, `.uns["spatial"]['spaceranger100_count_30458_ST8059049_mm10-3_0_0_premrna']['images']` etc.
