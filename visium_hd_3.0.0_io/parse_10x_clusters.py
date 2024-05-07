##
from anndata import AnnData
import pandas as pd
from pathlib import Path

parent_dir = Path(__file__).parent

f = parent_dir / 'data/square_016um/analysis/clustering/gene_expression_graphclust/clusters.csv'
assert f.exists()

df = pd.read_csv(f)
df
pass

##
