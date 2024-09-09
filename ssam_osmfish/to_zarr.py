import ssam
import pandas as pd
import numpy as np
import h5py
import scanpy as sc
import spatialdata as sd
import anndata as ad
import dask.dataframe as dd

##Before run this file, SSAM should be runned by to_zarr_with_plots.ipynb
ds = ssam.SSAMDataset("data/processed/ssamdataset-osmFISH.zarr")
analysis = ssam.SSAMAnalysis(ds, ncores=40, verbose=True)

pixel_per_um = 15.3846 # from BioRxiv paper
um_per_pixel = 1.0 / pixel_per_um

blacklists = ['Cnr1_Hybridization4', 'Plp1_Hybridization4', 'Vtn_Hybridization4',
              'Klk6_Hybridization5', 'Lum_Hybridization9', 'Tbr1_Hybridization11']

with h5py.File("./data/raw/mRNA_coords_raw_counting.hdf5", 'r') as f:
    # Extract gene names (keys)
    gene_names = list(f.keys())
    # Extract x, y coordinates for each gene and store them along with the gene name
    data_list = []
    for gene in gene_names:
        if gene in blacklists:
            continue
        xx, yy = f[gene][:].T
        gene = gene.split("_")[0]
        
        if gene == 'Tmem6':
            gene = 'Tmem2'
        elif gene == 'Kcnip':
            gene = 'Kcnip2'

        for x, y in zip(xx, yy):
            data_list.append([gene, x * um_per_pixel, y * um_per_pixel])

# Create a pandas DataFrame
df = pd.DataFrame(data_list, columns=['gene', 'x', 'y']).set_index('gene')

# Convert the DataFrame to a Dask DataFrame
dask_df = dd.from_pandas(df, npartitions=2) 
single_molecule = sd.models.PointsModel.parse(dask_df)

analysis.compute_cell_by_gene_matrix(df)
cell_type = []
for i in range(len(ds.cell_by_gene_matrix)):
    cell_type.append(int(np.unique(ds.watershed_celltype_maps[ds.watershed_segments == i]))+1)
domain = []

unique_domains  = [np.unique(ds.inferred_domains[ds.watershed_segments == i], return_counts=True) for i in range(len(ds.cell_by_gene_matrix))]

for unique_values, counts in unique_domains:
    max_count_index = np.argmax(counts)
    most_frequent_value = unique_values[max_count_index]
    domain.append(most_frequent_value)
    
table = sc.AnnData(ds.cell_by_gene_matrix)

table.var_names = ds.genes
table.obs["region"] = ["watershed_segments" for i in range(len(ds.cell_by_gene_matrix))]
table.obs['cell_id'] = np.arange(len(ds.cell_by_gene_matrix))+1
table.obs['cell_type'] = np.array(cell_type)
table.obs['domain'] = np.array(domain)+1

expression = sd.models.TableModel.parse(
    adata=table,
    region="watershed_segments",
    region_key='region',
    instance_key="cell_id",
)

celltype_maps = sd.models.Image2DModel.parse(ds.celltype_maps, dims=( 'x', 'y' ,'c'))
filtered_celltype_maps = sd.models.Image2DModel.parse(ds.filtered_celltype_maps, dims=('x', 'y','c'))
watershed_segments = sd.models.Labels2DModel.parse(ds.watershed_segments+1, dims=('x', 'y'))
inferred_domains = sd.models.Labels2DModel.parse(np.squeeze(ds.inferred_domains)+1, dims=('x', 'y'))

sdata = sd.SpatialData(points={"single_molecule": single_molecule}, images = {'celltype_maps': celltype_maps, 'filtered_celltype_maps': filtered_celltype_maps},labels={'watershed_segments': watershed_segments,'inferred_domains': inferred_domains, 'inferred_domains_cells': inferred_domains_cells},table=expression)
sdata.write('data.zarr') 