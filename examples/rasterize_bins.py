##
import time
from spatialdata import read_zarr, rasterize_bins

f = '../visium_hd_3.0.0_io/data.zarr'

sdata = read_zarr(f)

##
# bins = '002'
bins = '016'
sdata[f'square_{bins}um'].X = sdata[f'square_{bins}um'].X.tocsc()
var_names = ["Epcam", "Cd8a", "Cd4"]
start = time.time()
rasterized = rasterize_bins(
    sdata,
    f"Visium_HD_Mouse_Small_Intestine_square_{bins}um",
    f"square_{bins}um",
    "array_col",
    "array_row",
    # value_key=var_names,
)
print(f"rasterize: {time.time() - start}")
print('done')

from napari_spatialdata import Interactive
sdata['rasterized'] = rasterized
interactive = Interactive(sdata, headless=True)
interactive.add_element(element='rasterized', element_coordinate_system='global', view_element_system=True)
interactive.run()

##
