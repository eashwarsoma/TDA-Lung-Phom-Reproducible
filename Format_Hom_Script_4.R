#Script Number 4
#Converts Raw Homology to Feature Curve
#Consolidates Clinical Data and Feature Curve into 1 RDS file
library(dplyr)
library(anytime)


#Getting file names for homolgy
hom.csv <- list.files(path="./Python_Hom/", full.names = TRUE, all.files	
                      = FALSE)

#Getting the actual homology
hom.list <- lapply(hom.csv, read.csv, header = FALSE, nrows = 1000000)

#naming the columns properly
for (i in 1:length(hom.list)) {
  hom.list[[i]] <- hom.list[[i]][c(3,1,2)]
  
  colnames(hom.list[[i]]) <- c("dimension", "birth", "death")
}

#Removing the erraneous inf value
for (i in 1:length(hom.list)) {
  
  ind <- which(hom.list[[i]] == Inf, arr.ind = TRUE)
  
  hom.list[[i]] <- hom.list[[i]][-ind[1], ]
  
}

#Getting absolute min and max
abs.min <- min(unlist(lapply(hom.list, function (x) min(x[,c(2,3)]))))
abs.max <- max(unlist(lapply(hom.list, function (x) max(x[,c(2,3)]))))


#Normalizing all HU unit values to 0 to 1
#Note for normalization:
#Since the final product is a "topological feature curve"
#Which compared # of features vs filtration distance
#It doesn't matter how or where normalization (normalizing DICOM, normalizing here, etc.)
#As long as it is a linear transformation, the total number of topological features
#And their relative position on the filtration x axis will be the same
for (i in 1:length(hom.list)) {
  hom.list[[i]] <- cbind(hom.list[[i]][,1], 
                         ((hom.list[[i]][,c(2,3)] - abs.min)/(abs.max - abs.min)))
  colnames(hom.list[[i]]) <- c("dimension", "birth", "death")
}


#Getting the original scan names to name the phom
#Removing the file part of the name
name1 <- sapply(hom.csv, function (x) gsub("./python_hom_v2//", "" , x))
name2 <- sapply(name1, function (x) gsub(".npy.csv", "" , x))

#Naming the list
names(hom.list) <- name2


#Function to count total features that exist at each threshold value
cubical.total.feat.counter <- function(phom, min, max, res) {
  by <- (max-min)/res
  seq.to.check <- seq(min, max, by)
  mat <- cbind(rep(0, length(seq.to.check)), seq.to.check)
  for (i in 1:length(seq.to.check)) {
    tally <- sum(seq.to.check[i] >= phom[,2] & 
                   seq.to.check[i] <= phom[,3])
    
    mat[i, 1] <- tally
  }
  return(mat[,c(2,1)])
}

#Creating a function to count based on individual dim feat
#Not used, but function is here for posterity
featcounter.vector.cubical.act <- function(phom, min, max, res, feat) {
  
  phom <- phom[which(phom[,1] == feat),]
  
  by <- (max-min)/res
  seq.to.check <- seq(min, max, by)
  mat <- cbind(rep(0, length(seq.to.check)), seq.to.check)
  for (i in 1:length(seq.to.check)) {
    tally <- sum(seq.to.check[i] >= phom[,2] & 
                   seq.to.check[i] <= phom[,3])
    
    mat[i, 1] <- tally
  }
  return(mat[,c(2,1)])
}


#Read in Radiomics Clinic Data
clinic.rad.lung <- read.csv("rad_clinic.csv")

#There is not a perfect match from clinic csv and radiomic data
#Use grep function to search for patient IDs located in the scan file name
#if patient ID is in the list of radiomics files, we keep note of the index
rad.ind <- lapply(clinic.rad.lung$PatientID, function (x) grep(x, 
                                                               hom.csv, value = TRUE, 
                                                               ignore.case = TRUE))
rad.ind <- lapply(rad.ind, function (x) ifelse(length(x) > 0, 1, 0))
rad.ind <- do.call(rbind, rad.ind)

#Keep only the clinical data for which there is radiomic data
#Scan 128 doesn't exist in radiomic dataset, so it was the only 
#clinical data point removed
clinic.rad.lung <- clinic.rad.lung[which(rad.ind > 0), ]

#Repeating for radiogenomics
clinic.radg.lung <- read.csv("radg_clinic.csv")

#There is not a perfect match from clinic csv and radiogenomics data
#Use grep function to search for patient IDs located in the scan file name
#if patient ID is in the list of radiomics files, we keep note of the index
radg.ind <- lapply(clinic.radg.lung$Case.ID, function (x) grep(x, 
                                                               hom.csv, value = TRUE, 
                                                               ignore.case = TRUE))
radg.ind <- lapply(radg.ind, function (x) ifelse(length(x) > 0, 1, 0))
radg.ind <- do.call(rbind, radg.ind)

#Keep only the clinical data for which there is radiomic data
#Many Clinical Data Exist for Patients without scan/segment info 
clinic.radg.lung <- clinic.radg.lung[which(radg.ind > 0), ]



#Creating an overall stage variable
clinic.radg.lung <- clinic.radg.lung %>%
    mutate(
    Overall.Stage = case_when(
      Pathological.T.stage == "Tis" & Pathological.N.stage == "N0" & Pathological.M.stage == "M0" ~ "0",
      Pathological.T.stage == "T1a" & Pathological.N.stage == "N0" & Pathological.M.stage == "M0" ~ "I",
      Pathological.T.stage == "T1b" & Pathological.N.stage == "N0" & Pathological.M.stage == "M0" ~ "I",
      Pathological.T.stage == "T2a" & Pathological.N.stage == "N0" & Pathological.M.stage == "M0" ~ "I",
      Pathological.T.stage == "T2b" & Pathological.N.stage == "N0" & Pathological.M.stage == "M0" ~ "II",
      Pathological.T.stage == "T1a" & Pathological.N.stage == "N1" & Pathological.M.stage == "M0" ~ "II",
      Pathological.T.stage == "T1b" & Pathological.N.stage == "N1" & Pathological.M.stage == "M0" ~ "II",
      Pathological.T.stage == "T2a" & Pathological.N.stage == "N1" & Pathological.M.stage == "M0" ~ "II",
      Pathological.T.stage == "T2b" & Pathological.N.stage == "N1" & Pathological.M.stage == "M0" ~ "II",
      Pathological.T.stage == "T3" & Pathological.N.stage == "N0" & Pathological.M.stage == "M0" ~ "II",
      Pathological.T.stage == "T1a" & Pathological.N.stage == "N2" & Pathological.M.stage == "M0" ~ "IIIa",
      Pathological.T.stage == "T1b" & Pathological.N.stage == "N2" & Pathological.M.stage == "M0" ~ "IIIa",
      Pathological.T.stage == "T2a" & Pathological.N.stage == "N2" & Pathological.M.stage == "M0" ~ "IIIa",
      Pathological.T.stage == "T2b" & Pathological.N.stage == "N2" & Pathological.M.stage == "M0" ~ "IIIa",
      Pathological.T.stage == "T3" & Pathological.N.stage == "N1" & Pathological.M.stage == "M0" ~ "IIIa",
      Pathological.T.stage == "T3" & Pathological.N.stage == "N2" & Pathological.M.stage == "M0" ~ "IIIa",
      Pathological.T.stage == "T4" & Pathological.N.stage == "N0" & Pathological.M.stage == "M0" ~ "IIIa",
      Pathological.T.stage == "T5" & Pathological.N.stage == "N1" & Pathological.M.stage == "M0" ~ "IIIa",
      Pathological.T.stage == "T4" & Pathological.N.stage == "N2" & Pathological.M.stage == "M0" ~ "IIIb",
      Pathological.N.stage == "N3" & Pathological.N.stage == "M0" ~ "IIIb",
      Pathological.M.stage == "M1a" ~ "IV",
      Pathological.M.stage == "M1b" ~ "IV"
    )
  )



#Getting Getting Zero Dim Feature Curve
feat.0.curve <- lapply(hom.list, featcounter.vector.cubical.act, 0, 1, 1000, 0)

#Getting Getting One Dim Feature Curve
feat.1.curve <- lapply(hom.list, featcounter.vector.cubical.act, 0, 1, 1000, 1)

#Getting Getting Two Dim Feature Curve
feat.2.curve <- lapply(hom.list, featcounter.vector.cubical.act, 0, 1, 1000, 2)

#Getting Getting Total Dim Feature Curve
feat.tot.curve <- lapply(hom.list, cubical.total.feat.counter, 0, 1, 1000)

#Ensuring numerical data is numerical
clinic.rad.lung$Survival.time <- as.numeric(as.character(clinic.rad.lung$Survival.time))
clinic.radg.lung$Time.to.Death..days. <- as.numeric(as.character(clinic.radg.lung$Time.to.Death..days.))

#Creating new Clinical Survival Based off CT Date - Date of death
#For Radiogenomics set, survival time was only included for patients who died
#Unlike the radiomics set where time was included for last follow up
#This alternate survival was calculated so that there is consistency between data sets
#Modifying Dead/Alive to Correct number
clinic.radg.lung$mod.surv <- anydate(clinic.radg.lung$Date.of.Last.Known.Alive) - 
                              anydate(clinic.radg.lung$CT.Date)

clinic.radg.lung$Survival.Status <- ifelse(clinic.radg.lung$Survival.Status == "Dead", 1, 0)

#Combining data into one master list with all feature curves and clinical data
#Adding radiomics feature curves and clinic data
complete.lung.results <- list()
for (i in 1:421) {
  complete.lung.results[[i]] <- c(feat.0.curve[i], feat.1.curve[i], 
                                  feat.2.curve[i], feat.tot.curve[i],
                                  clinic.rad.lung[i, ])
}

#Adding radiogenomics feature curves and clinic data
for (i in 422:565) {
  complete.lung.results[[i]] <- c(feat.0.curve[i], feat.1.curve[i], 
                                  feat.2.curve[i], feat.tot.curve[i],
                                  clinic.radg.lung[i-421, ])
}

#Figuring out survival quartiles
surv1 <- sapply(complete.lung.results[1:421], function (x) x$Survival.time)
surv2 <- sapply(complete.lung.results[422:565], function (x) x$mod.surv)
survan <- na.omit(c(surv1, surv2))
quart <- quantile(survan)

#Storing Radiomics names and radiogenomics
#Modifying variables for consistency 
rad.names <- names(complete.lung.results[[1]])
rad.names[1:4] <- c("Zero.Feat.Count", "One.Feat.Count", "Two.Feat.Count", "Tot.Feat.Count")

radg.names <- names(complete.lung.results[[422]])
radg.names[1:4] <- c("Zero.Feat.Count", "One.Feat.Count", "Two.Feat.Count", "Tot.Feat.Count")
radg.names[5] <- "PatientID"
radg.names[41] <- "orig.Survival.time"
radg.names[39] <- "deadstatus.event"
radg.names[46] <- "Survival.time"


#Naming Our variables
for (i in 1:421) {
  names(complete.lung.results[[i]]) <- rad.names
}
for (i in 422:565) {
  names(complete.lung.results[[i]]) <- radg.names
}




#Adding A survival quartile Variable
complete.lung.results <- lapply(complete.lung.results, function (x) 
                                          c(x, ifelse(x[["Survival.time"]] <= quart[2], "surv25",
                                                  ifelse(x[["Survival.time"]] <= quart[3], "surv50",
                                                         ifelse(x[["Survival.time"]] <= quart[4], "surv75",
                                                                ifelse(x[["Survival.time"]] <= quart[5], "surv100", NA))))))


#Name the phom represetnations and survival variable
rad.names <- c(rad.names, "surv.cat")
radg.names <- c(radg.names, "surv.cat")


#Naming Our variables with survival now
for (i in 1:421) {
  names(complete.lung.results[[i]]) <- rad.names
}
for (i in 422:565) {
  names(complete.lung.results[[i]]) <- radg.names
}

#Ensuring survival is treated as a factor
for (i in 1:565) {
  complete.lung.results[[i]][["surv.cat"]] <- factor(complete.lung.results[[i]][["surv.cat"]], 
                                                     c("surv25", "surv50", "surv75", "surv100"))
}

#Saving Result
saveRDS(complete.lung.results, "complete_lung_results.rds")


