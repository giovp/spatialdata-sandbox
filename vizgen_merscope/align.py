from pathlib import Path

import geopandas
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

path = Path(__file__).parents[0] / "data" / "vizgen_sample"
images_dir = path / "images"

microns_to_pixels = np.genfromtxt(images_dir / "micron_to_mosaic_pixel_transform.csv")

image = np.stack(
	[
		np.array(Image.open(images_dir / name))
		for name in [
			"mosaic_DAPI_z3.tif",
			"mosaic_stain1_z3.tif",
			"mosaic_stain2_z3.tif",
		]
	],
	axis=2,
)

transcripts = pd.read_csv(path / "detected_transcripts.csv")
coords = np.concatenate(
	[transcripts[["global_x", "global_y"]].values, np.ones((len(transcripts), 1))],
	axis=1,
)

pixels_coords = (coords @ microns_to_pixels.T)[:, :2]
pixels_coords = pixels_coords[
	np.random.choice(len(pixels_coords), 10_000, replace=False)
]  # Subsample 10k transcripts

geo_df = geopandas.read_parquet(path / "vpt-run" / "cellpose_micron_space.parquet")

### Plot image + transcripts + boundaries
plt.figure(figsize=(12, 12))
plt.imshow(image.clip(0, 40_000) / 40_000)
plt.scatter(pixels_coords[:, 0], pixels_coords[:, 1], s=1, c="k")

for geom in geo_df[geo_df.ZIndex == 0].Geometry:
	x, y = geom.geoms[0].exterior.coords.xy
	x, y = (
		np.stack([x, y], axis=1) @ microns_to_pixels[:2, :2].T + microns_to_pixels[:2, 2]
	).T
	plt.plot(x, y, "w", linewidth=0.5)

plt.show()
plt.savefig("cropped_images_aligned.png")

