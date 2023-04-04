##
import os
from spatialdata_io import mcmicro
import shutil

##
f = 'data/exemplar-001'
assert os.path.isdir(f)
sdata = mcmicro(f, dataset_id='exemplar-001')

outfile = 'exemplar_001.zarr'
if os.path.exists(outfile):
    shutil.rmtree(outfile)
sdata.write(outfile)

from napari_spatialdata import Interactive
Interactive(sdata)