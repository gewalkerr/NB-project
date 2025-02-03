#! usr/bin/env R

# Loading libraries
library(nicheDE)
library(Seurat)
library(ggplot2)
library(Matrix)

# Setting file paths
setwd('/project/data/gew123/Spatial/NicheDE')
reference_dataset <- "/project/data/gew123/Spatial/Annotated_objects/atlas_annotated_seurat.rds"
spatial_dataset <- "/project/data/gew123/Spatial/Annotated_objects/spatial_annotated_seurat.rds"
# Change to selected tissue for analysis
tissue_to_analyze <- 'tissue_01'

### STEP 1: Generate Average Expression Profile from reference ###
# Read in reference
atlas <- readRDS(file = reference_dataset)
# Checking cell types
unique_cell_types <- levels(atlas$annotations_final)

# Set idents
Idents(object = atlas) <- atlas$annotations_final
# Drop levels
atlas$annotations_final = droplevels(atlas$annotations_final)

# Generate library_matrix from reference
library_matrix <- CreateLibraryMatrixFromSeurat(atlas, assay = 'originalexp')

### STEP 2: Read in spatial data and subset to single tissue ###
# Read in spatial object
spatial <- readRDS(file = spatial_dataset)

# Subset
tissue = subset(x = spatial, subset = tissue == tissue_to_analyze)

# Get spatial counts matrix (INPUT VARIABLE) 
counts = tissue@assays$originalexp@counts
# Saving gene and cell names
spots = colnames(counts)
genes = rownames(counts)
# Format to sparse
counts = as(counts,'sparseMatrix')
# Ensure proper row/col names
colnames(counts) = spots
rownames(counts) = genes

# Extracting cell centroids from global coordinates
cell_centroids = data.frame(y = tissue$CenterX_global_px, x = tissue$CenterY_global_px, cell = colnames(tissue))
centroid_data = list("centroids" = CreateCentroids(cell_centroids))
coords = CreateFOV(coords = centroid_data, type = c("centroids"), assay = "RNA")
# Add coords to metadata
tissue[["global"]] = coords

# Get coordinates df
centroids = GetTissueCoordinates(tissue)
# Set cells as rownames
rownames(centroids) = centroids$cell
centroids <- centroids[,-3]

### STEP 3: Build and subsequently fill a deconvolution matrix ###
# Initialize empty deconvolution matrix
cell_types = Idents(atlas)
cell_types = droplevels(cell_types)
unique_celltypes = unique(cell_types)
lib_sizes = rep(NA,length(unique_celltypes))

# Get cell specific library sizes and CTs
cell_libs = colSums(atlas)
# Get cell types of cells from reference 
cell_ct = Idents(atlas)

# Initialize the matrix 
deconv_matrix = as.data.frame(matrix(0, nrow(centroids), length(rownames(library_matrix))))

# Set colnames as cell types
colnames(deconv_matrix) = rownames(library_matrix)
# Set rownames as spots (cells)
rownames(deconv_matrix) = rownames(centroids)

# Confirming colnames
colnames(deconv_matrix)

# Filling deconvolution matrix
deconv_matrix$'Neuroblasts' <- ifelse(tissue$fine_annotations == "Neuroblasts", 1, 0)
deconv_matrix$'Cytotoxic (CD8+) T cells' <- ifelse(tissue$fine_annotations == "Cytotoxic (CD8+) T cells", 1, 0)
deconv_matrix$'Macrophages/monocytes' <- ifelse(tissue$fine_annotations == "Mesenchymal cells", 1, 0)
deconv_matrix$'Mesenchymal cells' <- ifelse(tissue$fine_annotations == "B cells", 1, 0)
deconv_matrix$'B cells' <- ifelse(tissue$fine_annotations == "B cells", 1, 0)
deconv_matrix$'Bridge cells' <- ifelse(tissue$fine_annotations == "Bridge cells", 1, 0)
deconv_matrix$'Endothelial cells' <- ifelse(tissue$fine_annotations == "Endothelial cells", 1, 0)
deconv_matrix$'Connecting progenitor cells' <- ifelse(tissue$fine_annotations == "Connecting progenitor cells", 1, 0)
deconv_matrix$'Naïve T cells' <- ifelse(tissue$fine_annotations == "Naïve T cells", 1, 0)
deconv_matrix$'Late neuroblasts' <- ifelse(tissue$fine_annotations == "Late neuroblasts", 1, 0)
deconv_matrix$'Plasma cells' <- ifelse(tissue$fine_annotations == "Plasma cells", 1, 0)
deconv_matrix$'Proliferating T cells' <- ifelse(tissue$fine_annotations == "Proliferating T cells", 1, 0)
deconv_matrix$'Myofibroblasts' <- ifelse(tissue$fine_annotations == "Myofibroblasts", 1, 0)
deconv_matrix$'SCPs' <- ifelse(tissue$fine_annotations == "SCPs", 1, 0)
deconv_matrix$'Adrenal cortex' <- ifelse(tissue$fine_annotations == "Adrenal cortex", 1, 0)
deconv_matrix$'Chromaffin cells' <- ifelse(tissue$fine_annotations == "Chromaffin cells", 1, 0)
deconv_matrix$'NK cells' <- ifelse(tissue$fine_annotations == "NK cells", 1, 0)
deconv_matrix$'Mast cells' <- ifelse(tissue$fine_annotations == "Mast cells", 1, 0)
deconv_matrix$'Cytotoxic (CD8+/IGKC+) T cells' <- ifelse(tissue$fine_annotations == "Cytotoxic (CD8+/IGKC+) T cells", 1, 0)
deconv_matrix$'Helper (CD4+) T cells' <- ifelse(tissue$fine_annotations == "Helper (CD4+) T cells", 1, 0)
deconv_matrix <- as.matrix(deconv_matrix) # (INPUT VARIABLE)

### STEP 4: Build a NicheDE Object ###
# Get centroids again
centriods = GetTissueCoordinates(tissue)
# Set rownames as cells
rownames(centroids) = centroids$cell
# Remove cell column
centroids <- centroids[,-3]

# Ensuring columns of deconv matrix match rows of reference
colnames(deconv_est) == rownames(library_matrix) 

# Extract transposed spatial counts matrix
spatial_counts <- t(tissue@assays$originalexp@counts)

# CreateNicheDEObject takes 4 input variables:
    # 1) 'spatial_counts' : spatial raw counts matrix (cells x genes)
    # 2) 'centroids' : cell (spot) coordinates matrix
    # 3) 'library_matrix' : average expression profile matrix from reference
    # 4) 'deconv_matrix' : deconvolution matrix of cell type predictions

NDE_obj <- CreateNicheDEObject(
  spatial_counts, 
  centroids, 
  library_matrix, 
  deconv_matrix, 
  sigma = c(10, 100, 300, 1000, 10000, 25000),
  Int = T
)



