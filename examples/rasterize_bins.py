##
from spatialdata import read_zarr, rasterize_bins

f = '../visium_hd_3.0.0_io/data.zarr'

sdata = read_zarr(f)

##
sdata['square_002um'].X = sdata['square_002um'].X.tocsc()
while True:
    var_names = ["Epcam", "Cd8a", "Cd4"]
    rasterize_bins(
        sdata,
        "Visium_HD_Mouse_Small_Intestine_square_002um",
        "square_002um",
        "array_col",
        "array_row",
        value_key=var_names,
    )
    print('done')
