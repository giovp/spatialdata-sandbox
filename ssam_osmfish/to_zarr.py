import ssam
import pandas as pd
import numpy as np
import h5py
import scanpy as sc
import spatialdata as sd
import anndata as ad

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

analysis.compute_cell_by_gene_matrix(df)
new_element = [0 for i in range(len(ds.genes))]
cell_by_gene_matrix = np.insert(ds.cell_by_gene_matrix, 0, new_element, axis=0)

cell_type = [0]
for i in range(len(cell_by_gene_matrix)-1):
    cell_type.append(int(np.unique(ds.watershed_celltype_maps[ds.watershed_segments == i]))+1)

main_domain = [0]
secondary_domain = [0]
additional_domains = [0]

unique_domains = [np.unique(ds.inferred_domains[ds.watershed_segments == i]) for i in range(len(cell_by_gene_matrix) - 1)]

for k in unique_domains:
    if len(k) == 1:
        main_domain.append(k[0] + 1)
        secondary_domain.append(k[0] + 1)
        additional_domains.append(0)
    elif len(k) == 2:
        main_domain.append(k[0] + 1)
        secondary_domain.append(k[1] + 1)
        additional_domains.append(0)
    elif len(k) > 2:
        main_domain.append(k[0] + 1)
        secondary_domain.append(k[1] + 1)
        additional_domains.append(list(k[2:] + 1))
     

table = sc.AnnData(cell_by_gene_matrix)

table.var_names = ds.genes
table.obs["region"] = ["watershed_segments" for i in range(len(cell_by_gene_matrix))]
table.obs['cell_id'] = np.arange(len(cell_by_gene_matrix))+1
table.obs['cell_type'] = np.array(cell_type)
table.obs['main_domain'] = np.array(main_domain)
table.obs['secondary_domain'] = np.array(secondary_domain)

expression = sd.models.TableModel.parse(
    adata=table,
    region="watershed_segments",
    region_key='region',
    instance_key="cell_id",
)

celltype_maps = sd.models.Image2DModel.parse(ds.celltype_maps[::-1], dims=('x', 'y', 'c'))
filtered_celltype_maps = sd.models.Image2DModel.parse(ds.filtered_celltype_maps[::-1], dims=('x', 'y', 'c'))
watershed_segments = sd.models.Labels2DModel.parse(ds.watershed_segments[::-1]+1, dims=('x', 'y'))
inferred_domains = sd.models.Labels2DModel.parse(np.squeeze(ds.inferred_domains[::-1]+1), dims=('x', 'y'))
inferred_domains_cells = sd.models.Labels2DModel.parse(np.squeeze(ds.inferred_domains_cells[::-1]+1), dims=('x', 'y'))

sdata = sd.SpatialData(images = {'celltype_maps': celltype_maps, 'filtered_celltype_maps': filtered_celltype_maps},labels={'watershed_segments': watershed_segments,'inferred_domains': inferred_domains, 'inferred_domains_cells': inferred_domains_cells},table=expression)
sdata.write('data.zarr') 