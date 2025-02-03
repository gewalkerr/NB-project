#!/usr/bin/env python

# Catch any errors
import sys
import os

log_file = 'c2l_errors_FINAL2.log'
sys.stdout = open(log_file, 'w') # Logs go to directory where script is
sys.stderr = open(log_file, 'w')

# Packages
import scanpy as sc
import squidpy as sq
import numpy as np
import cell2location as c2l
import matplotlib
import scipy
from matplotlib import rcParams
import matplotlib.pyplot as plt
import scipy.sparse

# Settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor="white")
sc.settings.figdir = '/project/data/gew123/Spatial/Deconvolution'

# Figure settings
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs
rcParams['figure.figsize'] = 6,6

# Set up results directory and paths to results folders
results = '/project/data/gew123/Spatial/Deconvolution/Results_final'

ref_run_name = f'{results}/Reference_signatures'
run_name = f'{results}/Models'

# Read in reference and spatial objects
## Reference has already been filtered and formatted for the model
adata_ref_fil = sc.read_h5ad('/project/data/gew123/Spatial/Objects/Reference_c2l_filtered_FINAL.h5ad')
adata_st = sc.read_h5ad('/project/data/gew123/Spatial/Objects/spatial_clustered_unlabelled.h5ad')

# Get rid of MT genes in spatial
adata_st.var["MT_gene"] = [gene.startswith("MT-") for gene in adata_st.var_names]
adata_st.obsm["MT"] = adata_st[:, adata_st.var["MT_gene"].values].X.toarray()
adata_st = adata_st[:, ~adata_st.var["MT_gene"].values]


# (1) Generating reference gene signatures

# Prepare anndata for the regression model
c2l.models.RegressionModel.setup_anndata(adata=adata_ref_fil,
                        # 10X reaction / sample / batch
                        batch_key='Source',
                        # cell type, covariate used for constructing signatures
                        labels_key='annotations_final',
                       )

# Create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref_fil)

# view anndata_setup as a sanity check
print(mod.view_anndata_setup())

print('Training the model....')
# Training the model
mod.train(max_epochs=250)
print('Model training completed. Saving model.')

# Exporting the model to results file
adata_ref_fil = mod.export_posterior(
    adata_ref_fil, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/reference_signatures.h5ad"
adata_ref_fil.write(adata_file)

print('Script ran successfully.')