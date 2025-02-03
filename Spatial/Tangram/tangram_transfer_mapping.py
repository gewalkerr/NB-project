#! usr/bin/env python

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

# Read in spatial and sc objects
spatial_object_path = "/project/data/gew123/Spatial/Imputation/spatial_mapped.h5ad"
reference_object_path = "/project/data/gew123/Processed_h5ad/Jansky_processed.h5ad"

print('Reading in spatial and sc datasets...')
adata_st = sc.read_h5ad(spatial_object_path)
adata_sc = sc.read_h5ad(reference_object_path)
print('Done.')

# Read in mapping results
print("Reading in mapping results...")
ad_map = sc.read_h5ad("/project/data/gew123/Spatial/Imputation/mapped_tangram_object.h5ad")
print('Done.')


print('Projecting gene expression to spatial data...')
ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=adata_sc)
print('Done')

print('Saving new spatial object...')
adata_ge.write('/project/data/gew123/Spatial/Imputation/adata_ge.h5ad')
print('Done, script ran successfully.')






