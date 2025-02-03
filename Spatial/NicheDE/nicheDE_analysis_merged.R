#! usr/bin/env R

# Loading libraries
library(nicheDE)
library(Seurat)
library(ggplot2)
library(Matrix)

# Setting file paths
setwd('/project/data/gew123/Spatial/NicheDE')
nde_object_dir <- "/project/data/gew123/Spatial/NicheDE/NDE_objects/"


# Read in object
nde_object_path <- paste0(nde_object_dir, "NDE_merged.rds")
NDE_obj = readRDS(nde_object_path)

# Calculate effective niche size
print("Calculating effective niche size:")
NDE_obj = CalculateEffectiveNicheLargeScale(NDE_obj,batch_size = 1000, cutoff = 0.05)
print("Done.") 


# Run NicheDE
print("Running NicheDE....")
# Start cluster
cl <- parallel::makeCluster(3, outfile = "outfile.txt")
# Register cluster
doParallel::registerDoParallel(cl)

NDE_obj = niche_DE(NDE_obj,
                num_cores=4,
                outfile="outfile.txt",
                C = 150, # Default
                M = 10, 
                gamma = 0.8,
                print = T, # Tracks progress of nicheDE in outfile
                Int = T, # Default, Boolean for if data is count data
                batch = T, # Include indicator for each batch
                self_EN = T, # Should we consider self interactions of cells?
                G=1) 

doParallel::stopImplicitCluster()


print("Saving NDE object...")
saveRDS(NDE_obj, nde_object_path)
print("Results saved. Script ran successfully!")


