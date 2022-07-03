from pathlib import Path
from anndata import AnnData
import pandas as pd

mcmicro_dir = Path().resolve() / "data"
modality_dirs = [x for x in mcmicro_dir.glob("*") if x.is_dir()]
feature_dir_name = "quantification"

for x in modality_dirs:
  pilot_name = x.stem
  feature_dir = x / feature_dir_name
  csv_file = [*feature_dir.glob("*.csv")][0]
  
  # Data contains intensities, cell center location (global) and shape measurements. Intensities are treated as counts
  data_df = pd.read_csv(csv_file, header=0, index_col='CellID')
  # Delete copy of X_centroid and Y_centroid
  data_df.drop(columns=['column_centroid', 'row_centroid'], inplace=True)
  
  adata = AnnData(data_df.loc[:,:"X_centroid"].values[:,:-1], 
                  dtype=data_df.values.dtype, 
                  obs=data_df.loc[:,"X_centroid":], 
                  uns={"spatial":{}},
                  )
  
  # Var names set to be equal to the column names corresponding to antibody intensities
  split_index = data_df.columns.to_list().index("X_centroid")
  adata.var_names = data_df.columns[:split_index]
  
  adata.obsm["spatial"] = adata.obs[["X_centroid","Y_centroid"]]
  adata.uns["spatial"] = {
                          "images": {},
                          "segmentations": {},
                          }
	
  adata.write(x / f"{pilot_name}.h5ad")
