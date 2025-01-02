##
import spatialdata as sd
from napari_spatialdata import Interactive
from spatialdata_io import visium_hd

sdata = sd.read_zarr("data.zarr")
print(sdata)

##
# test rasterize bins after converting the data to Zarr
rasterized = sd.rasterize_bins(
    sdata,
    bins="Visium_HD_Human_Lung_Cancer_Fixed_Frozen_square_016um",
    table_name="square_016um",
    col_key="array_col",
    row_key="array_row",
    value_key=None,
    return_region_as_labels=True,
)
sdata['rasterized'] = rasterized
sd.rasterize_bins_link_table_to_labels(sdata=sdata, table_name="square_016um", rasterized_labels_name="rasterized")

##
# Interactive(sdata)

##
# test rasterize bins before converting the data to Zarr
sdata_labels = visium_hd('data', annotate_table_by_labels=True, bin_size=[16])
print(sdata_labels)

##
# Interactive(sdata_labels)
