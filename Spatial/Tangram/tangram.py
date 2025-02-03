import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
from anndata import AnnData
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import skimage
import seaborn as sns
import tangram as tg

sc.logging.print_header()
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor="white", dpi_save=400)

signature_genes = ['BASP1', 'BSG', 'CD147', 'CCDC125', 'CD9', 'DLG2', 'FNBP1', 'FRMD3', 'GABRB3',
              'GNB2L1', 'RACK1', 'HAPLN4', 'HEBP2', 'HSD17B12', 'IGSF10', 'IL11RA', 'IQCE', 'KCNQ3', 'TOX2']

# Read in objects
spatial_object_path = "/project/data/gew123/Spatial/Annotated_objects/spatial_annotated_final.h5ad"
reference_object_path = "/project/data/gew123/Processed_h5ad/Jansky_processed.h5ad"
tangram_object_path = "/project/data/gew123/Spatial/Imputation/mapped_tangram_object.h5ad"

adata_sc = sc.read_h5ad(reference_object_path)
adata_st = sc.read_h5ad(spatial_object_path)

# Get overlapping markers
markers = list(set.intersection(set(adata_sc.var_names), set(adata_st.var_names)))
len(markers)

# Filter genes
tg.pp_adatas(adata_sc, adata_st, genes=markers)

# Map cells to space
ad_map = tg.map_cells_to_space(
    adata_sc,
    adata_st,
    mode="cells",
    density_prior="uniform", # For single cell resolution
    num_epochs=500,
    device="cpu", # For CPU
)

# Save result
ad_map.write(tangram_object_path)