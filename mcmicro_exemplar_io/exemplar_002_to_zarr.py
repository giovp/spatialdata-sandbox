##
import os
from spatialdata_io import mcmicro
import shutil
import spatialdata as sd

##
f = 'data/exemplar-002'
assert os.path.isdir(f)
sdata = mcmicro(f)
print(sdata)

outfile = 'exemplar_002.zarr'
if os.path.exists(outfile):
    shutil.rmtree(outfile)
sdata.write(outfile)

sdata = sd.read_zarr(outfile)

# from napari_spatialdata import Interactive
# Interactive(sdata)
