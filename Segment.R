library(dplyr)

#File to make segments from DICOM Data
#Before running this script, make sure you have downloaded the DICOM images from google drive
#This code will use the NIFTII segments to create the cubical tumor region of interest
#We will save an RDS file
#And importantly, save numpy files to compute homology in python



#Function to get segmented slices from DICOM image and nifti segment
segment.slicer.adjusted <- function (dicom, segment) {
  
  dicom3d <- create3D(dicom)
  
  #Storing Header Data
  header <- dicom$hdr
  
  #Created a new corrected array with HU units
  hu.array <- array(dim = dim(dicom3d))
  for (i in 1:dim(dicom3d)[3]) {
    
    #Figuring out which index of values stores slope and intercept
    slop.ind <- which(header[[i]]$name == "RescaleSlope")
    int.ind <- which(header[[i]]$name == "RescaleIntercept")
    
    #Actually grabbing that value
    slop <- as.numeric(header[[i]]$value[slop.ind])
    int <- as.numeric(header[[i]]$value[int.ind])
    
    #Grabbing Slice i and performing linear rescale
    #Linear Rescaling to convert raw value to HU
    #Throw back Tuesday to 7th grade and y = mx + b
    hu.array[,,i] <- (dicom3d[,,i]*slop) + int
    
  }
  
  
  #Now that the standardized hu.array has been created
  #We have to correlate the nifti segment VOI with the 
  #Appropriate Slices and length + width indices
  #The ipp attr gives the appropriate slice location in 3d space
  #We have to match this to the nifti 3d space location
  attr.ipp <- attr(dicom3d, "ipp")
  
  #Getting the raw array data
  segment.arr <- img_data(segment)
  
  #Matrix of all the array indices of the segment that is marked with tumor
  bitmap.simp <- which(segment.arr > 0, arr.ind=TRUE)
  
  #Getting the indices translated to 3d space
  #We are interested particularly in the z coordinate which corresponds to slices
  trans <- voxelToWorld(bitmap.simp, segment)
  
  #Getting the real world coordinates of interest
  zslice.of.interest <- unique(trans[,3])
  
  #Getting the appropriate indices for the dicom
  #Adding a .002% tolerance due to rounding error
  #Adding absolute value as some locations are negative
  dicom.indices.of.interest <- sapply(zslice.of.interest, function(x) 
    which(abs(attr.ipp[,3]) <= abs(1.00002*x) & abs(attr.ipp[,3]) >= abs(.99998*x)))
  
  #Selecting the Z slices
  #Sorting them to get ride order
  hu.array.z.selected <- hu.array[,,sort(dicom.indices.of.interest)]
  
  #Converting Z selected array to list format
  hu.list.z.selected <- lapply(seq(dim(hu.array.z.selected)[3]), function(x) 
    hu.array.z.selected[ , , x])
  
  #Now to Select Appropriate Y parameter
  #Selects only the columns in the max and min range of y
  max.y <- max(bitmap.simp[,2])
  min.y <- min(bitmap.simp[,2])
  
  hu.list.zy.selected <- lapply(hu.list.z.selected, function (x) x[ , min.y:max.y])
  
  #Now to Select Appropriate X parameter
  #Selects only the columns in the max and min range of X
  max.x <- max(bitmap.simp[,1])
  min.x <- min(bitmap.simp[,1])
  
  hu.list.zyx.selected <- lapply(hu.list.zy.selected, function (x) x[min.x:max.x, ])
  
  return(hu.list.zyx.selected)
  
}

#Collecting all Radiomics Files
#gets a list of all Radiomics Lung DICOM Files
files.rad.dicom <- list.dirs(path="./Radiomics_DICOM/", full.names = TRUE,
                             recursive = FALSE)

#gets a list of all Radiomics Lung Nifti Segmentation Files
files.rad.nift <- list.files(path="./Radiomics_NIFTII/", full.names = TRUE, all.files	
                              = FALSE)

#Not all Dicoms have been segmented
#Removing the unsegmented DICOM
true.files.rad.dicom <- vector()
for (i in 1:length(files.rad.nift)) {
  true.files.rad.dicom[[i]] <- gsub(".dcm-1.nii.gz", "" , files.rad.nift[[i]]) %>%
  gsub("./Radiomics_NIFTII//", "" , .) %>% grep(., 
                                                files.rad.dicom, value = TRUE, 
                                               ignore.case = TRUE)
}


#Obtaining all the segments
#This step unfortunately takes some time due to the
#slowness of reading in DICOMS
#Make sure the RDS file is saved so this does not have to be repeated
segment.rad.slice.list <- list()
for (i in 1:length(true.files.rad.dicom)) {
  print(paste("Working on number", i, sep = " "))
  dicom <- readDICOM(true.files.rad.dicom[i])
  segment <- readNIfTI(files.rad.nift[i])
  segment.rad.slice.list[[i]] <- segment.slicer.adjusted(dicom, segment)
}

saveRDS(segment.rad.slice.list, "radiomics.segments.rds")

#Test out on some images
#Change the index 2 to whatever index to visualize the tumor
for (i in 1:length(segment.rad.slice.list[[2]])) {
  image(segment.rad.slice.list[[2]][[i]], col=grey(0:64/64), axes=FALSE, xlab="", ylab="")
}

#Repeating Same Steps for the Radiogenomics set
#gets a list of all Radiomics Lung DICOM Files
files.radg.dicom <- list.dirs(path="./Radiogenomics_DICOM/", full.names = TRUE,
                             recursive = FALSE)

#gets a list of all Radiomics Lung Nifti Segmentation Files
files.radg.nift <- list.files(path="./Radiogenomics_NIFTII/", full.names = TRUE, all.files	
                             = FALSE)


#Removing any mismatched DICOMs and segments 
#Doesn't apply to this data, but code shown here for consistence
true.files.radg.dicom <- vector()
for (i in 1:length(files.radg.nift)) {
  true.files.radg.dicom[[i]] <- gsub(".dcm-1.nii.gz", "" , files.radg.nift[[i]]) %>%
    gsub("./Radiogenomics_NIFTII//", "" , .) %>% grep(., 
                                                  files.radg.dicom, value = TRUE, 
                                                  ignore.case = TRUE)
}

#Obtaining all the segments
segment.radg.slice.list <- list()
for (i in 1:length(true.files.radg.dicom)) {
  print(paste("Working on number", i, sep = " "))
  dicom <- readDICOM(true.files.radg.dicom[i])
  segment <- readNIfTI(files.radg.nift[i])
  segment.radg.slice.list[[i]] <- segment.slicer.adjusted(dicom, segment)
}



saveRDS(segment.radg.slice.list, "radiogenomics.segments.rds")


#Test out on some images
for (i in 1:length(segment.radg.slice.list[[144]])) {
  image(segment.radg.slice.list[[144]][[i]], col=grey(0:64/64), axes=FALSE, xlab="", ylab="")
}



####Importing to Python####
#Going to python
# import numpy
#For loop function
#Don't load numpy with dplyr
#Read back in the RDS files if needed
devtools::unload("dplyr")

segment.rad.slice.list <- readRDS("radiomics.segments.rds")
segment.radg.slice.list <- readRDS("radiogenomics.segments.rds")


library(reticulate)

np = import("numpy")


for (i in 1:length(segment.rad.slice.list)) {
  
  names1 <- gsub(".dcm-1.nii.gz", "" , files.rad.nift[[i]]) 
  
  names2 <- gsub("./Radiomics_Niftii//", "" , names1)
  
  
  file.name <- paste("./Numpy_Arrays/", names2, ".npy", sep = "")
  
  np$save(file.name, r_to_py(segment.rad.slice.list[[i]]))
  
}



for (i in 1:length(segment.radg.slice.list)) {
  
  names1 <- gsub(".dcm-1.nii.gz", "" , files.radg.nift[[i]]) 
  
  names2 <- gsub("./Radiogenomic_Niftii//", "" , names1)
  
  
  file.name <- paste("./Numpy_Arrays/", names2, ".npy", sep = "")
  
  np$save(file.name, r_to_py(segment.radg.slice.list[[i]]))
  
}







