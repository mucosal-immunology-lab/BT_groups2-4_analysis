# Here we describe the process for preparing a pseudobulk dataset
# from the Integrated Human Lung Cell Database (full).
# Original dataset downloaded from:
# https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293
# and saved as "hlca.h5ad"

import numpy as np
import pyreadr
import pandas as pd
import anndata as ad
import h5py
import gc
import tables
import scanpy as sc
import decoupler as dc
from psutil import virtual_memory as mem
from scipy.sparse import csr_matrix

# Import the AnnData object
adata_full = sc.read_h5ad('hlca.h5ad')
adata = adata_full[adata_full.obs.disease.isin(['COVID-19', 'chronic rhinitis', 'cystic fibrosis', 'hypersensitivity pneumonitis', 'normal', 'pneumonia']) 
              & ~adata_full.obs.smoking_status.isin(['active', 'former', 'hist of marijuana use']) 
              & ~adata_full.obs.cell_type.isin(['native cell']) 
              & ~adata_full.obs.tissue.isin(['nose'])]

# Save the data (this is because you may need to restart python to reduce memory utilisation)
adata.write('hlca_lung.h5ad', compression = 'gzip')

# (Restart python and) import the AnnData object
adata = sc.read_h5ad('hlca_lung.h5ad')

# Get pseudo-bulk profile
pdata = dc.get_pseudobulk(adata, sample_col='study', groups_col='cell_type', mode='mean', min_cells=10, min_counts=100)
pdata_df = pdata.to_df()
gene_name = pdata.var.feature_name.astype('string')
pdata_df.columns = gene_name
pdata_df.to_csv('hlca_lung_pseudobulk_mean.csv')