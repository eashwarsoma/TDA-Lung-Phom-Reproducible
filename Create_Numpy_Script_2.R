#Script Number 2
#Reads in RDS Segmentation File
#Creates npy files for Python

#gets a list of all the segmentations made in the Segment.R file
segment.rad.slice.list <- readRDS("radiomics.segments.rds")
segment.radg.slice.list <- readRDS("radiogenomics.segments.rds")


#gets a list of all Radiomics Lung Nifti Segmentation File paths
files.rad.nift <- list.files(path="./Radiomics_NIFTII/", full.names = TRUE, all.files	
                             = FALSE)

#gets a list of all Radiogenomics Lung Nifti Segmentation File paths
files.radg.nift <- list.files(path="./Radiogenomics_NIFTII/", full.names = TRUE, all.files	
                              = FALSE)

#Don't load dplyr with numpy  
library(reticulate)

np = import("numpy")


for (i in 1:length(segment.rad.slice.list)) {
  
  names1 <- gsub(".dcm-1.nii.gz", "" , files.rad.nift[[i]]) 
  
  names2 <- gsub("./Radiomics_NIFTII//", "" , names1)
  
  
  file.name <- paste("./Numpy_Arrays/", names2, ".npy", sep = "")
  
  np$save(file.name, r_to_py(segment.rad.slice.list[[i]]))
  
}



for (i in 1:length(segment.radg.slice.list)) {
  
  names1 <- gsub(".dcm-1.nii.gz", "" , files.radg.nift[[i]]) 
  
  names2 <- gsub("./Radiogenomics_NIFTII//", "" , names1)
  
  
  file.name <- paste("./Numpy_Arrays/", names2, ".npy", sep = "")
  
  np$save(file.name, r_to_py(segment.radg.slice.list[[i]]))
  
}