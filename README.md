![SpatialData banner](https://github.com/scverse/spatialdata/blob/main/docs/_static/img/spatialdata_horizontal.png?raw=true)

# spatialdata-sandbox

The respository contains scripts to convert common formats into the `SpatialData` storage format. To learn more about the storage format, please see the [SpatialData design document](https://spatialdata.scverse.org/en/latest/design_doc.html). You can download the datasets [here](https://spatialdata.scverse.org/en/latest/tutorials/notebooks/datasets/README.html).

## How to add datasets

* Directory per dataset, named after the technology
* Single python executable `download.py` that dump all available data in `data` folder.
* Single python executable `write_zarr.py` that saves all available data in zarr file `data.zarr`.
* Be explicit about where the components are, e.g.
  * Table
  * Points
  * Images
  * Regions
* Have a readme.md with prose
