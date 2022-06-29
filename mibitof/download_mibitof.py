import squidpy as sq
import imageio
from pathlib import Path

PATH = Path().resolve()
SAVE = PATH / "data"
SAVE.mkdir(parents=True, exist_ok=True)
adata = sq.datasets.mibitof()

libraries = list(adata.uns["spatial"].keys())

for lib in libraries:
    img = adata.uns["spatial"][lib]["images"]["hires"]
    imageio.imwrite(SAVE / f"{lib}_image.png", img)

    seg = adata.uns["spatial"][lib]["images"]["segmentation"]
    imageio.imwrite(SAVE / f"{lib}_labels.png", seg)

    adata[adata.obs.library_id == lib].copy().write(SAVE / f"{lib}_table.h5ad")
