#Script Number 5 (Final Script)
#Reads in consolidated data
#Analyzes Survival by discrete and continuos method

library(plyr)

#Reading in Formatted Data and slices
complete.lung.results <- readRDS("complete_lung_results.rds")
segment.slices.rds <- c(readRDS("radiomics.segments.rds"), 
                        readRDS("radiogenomics.segments.rds"))

#Custom Clinic Splitter Function
#Function creates groups based off clinical variable
clinic.splitter.cust <- function(results, phom.rep, clinic.variable) {
  
  
  #Gets all factors for specified variable
  num.fact <- levels(results[[1]][[clinic.variable]])
  
  #Choosing only results with clinical variable and No Na
  clean.results <- lapply(results, function (x) ifelse(is.na(x[[clinic.variable]]) != TRUE, 
                                                       return(x), 
                                                       return(NULL))) %>% compact() 
  #List to store results
  group.list <- list()
  
  for (i in 1:length(num.fact)) {
    
    
    #Get the landscapes or dv
    group <- list()
    for (k in 1:length(clean.results)) {
      temp <- clean.results[[k]]
      if (temp[[clinic.variable]] == as.character(num.fact[i])) {
        group[[k]] <- temp[[phom.rep]]
      } else {
        group[[k]] <- NULL
      }
    }
    group <- group %>% compact()
    
    
    
    #Get the patient names of particular clinic variable
    names <- lapply(clean.results, function (x) ifelse(x[[clinic.variable]] == as.character(num.fact[i]), 
                                                       return(x[["PatientID"]]), 
                                                       return(NULL))) %>% compact() 
    
    #Adding the patient names to each clinical variable
    names(group) <- unlist(names)
    
    
    group.list[[i]]<- group
    
  }
  
  group.list <- lapply(group.list, na.omit)
  names(group.list) <- unlist(lapply(num.fact, function (x) paste(x)))
  return(group.list)
  
}


####Discrete analysis####
#Looking at total feature curve by the survival quartiles
surv <- clinic.splitter.cust(complete.lung.results, "Tot.Feat.Count", "surv.cat")

surv25 <- surv[[1]]
surv50 <- surv[[2]]
surv75 <- surv[[3]]
surv100 <- surv[[4]]

#Function to average feature curves and provide mean and standard error
feat.averager <- function(feat.list) {
  arr <- array(unlist(feat.list), dim = c(dim(feat.list[[1]]), length(feat.list)))
  feat.list.av <- apply(arr, c(1,2), mean)
  return(feat.list.av)
}

#Computing the average of each survival group
surv25.av <- feat.averager(surv25)
surv50.av <- feat.averager(surv50)
surv75.av <- feat.averager(surv75)
surv100.av <- feat.averager(surv100)

#Visualizing the survival groups
par(mfrow=c(2,2))
plot(surv25.av, xlim=c(0,1), ylim=c(0,2000)) + title("Surv25")
plot(surv50.av, xlim=c(0,1), ylim=c(0,2000)) + title("Surv50")
plot(surv75.av, xlim=c(0,1), ylim=c(0,2000)) + title("Surv75")
plot(surv100.av, xlim=c(0,1), ylim=c(0, 2000)) + title("Surv100")
dev.off()

library(pracma)

#Area and sde
surv25areas <- sapply(surv25, function (q) trapz(x = q[,1], y = q[,2]))
surv50areas <- sapply(surv50, function (q) trapz(x = q[,1], y = q[,2]))
surv75areas <- sapply(surv75, function (q) trapz(x = q[,1], y = q[,2]))
surv100areas <- sapply(surv100, function (q) trapz(x = q[,1], y = q[,2]))

par(mfrow=c(2,2))
hist(surv25areas, breaks = c(seq(0, 10000, 200)), ylim = c(0, 70))
hist(surv50areas, breaks = c(seq(0, 10000, 200)), ylim = c(0, 70))
hist(surv75areas, breaks = c(seq(0, 10000, 200)), ylim = c(0, 70))
hist(surv100areas, breaks = c(seq(0, 10000, 200)), ylim = c(0, 70))
dev.off()


#Comparing Survival Groups
wilcox.test(surv25areas, surv100areas, alternative = "two.sided")
wilcox.test(surv25areas, surv75areas, alternative = "two.sided")
wilcox.test(surv25areas, surv50areas, alternative = "two.sided")
wilcox.test(surv50areas, surv100areas, alternative = "two.sided")
wilcox.test(surv50areas, surv75areas, alternative = "two.sided")
wilcox.test(surv75areas, surv100areas, alternative = "two.sided")

#Mean
mean(surv25areas)
mean(surv50areas)
mean(surv75areas)
mean(surv100areas)


####Cox Analysis####
#Survival Time
surv.v <- do.call(rbind, lapply(complete.lung.results, function (x) x[["Survival.time"]]))

#Dead vs Alive Status
dead.alive.v <- do.call(rbind, lapply(complete.lung.results, function (x) x[["deadstatus.event"]]))

#Feature Number
area.v <- do.call(rbind, lapply(complete.lung.results, function (q) trapz(x = q[[1]][,1],
                                                                    y = q[[1]][,2])))

#Scaling area for cox model
area.v <- scale(area.v)

#Pixel Count 
pixelcount.v <- do.call(rbind, lapply(segment.slices.rds, function (x) 
  length(x)*dim(x[[1]])[1]*dim(x[[1]])[1]))
pixelcount.v <- scale(pixelcount.v)

cor(area.v, pixelcount.v)^2
plot(area.v, pixelcount.v)



#Age
age.rad <- do.call(rbind, lapply(complete.lung.results[1:421], function (x) x[["age"]])) 
age.radg <- do.call(rbind, lapply(complete.lung.results[422:565], function (x) 
                                  x[["Age.at.Histological.Diagnosis"]])) 
age.v <- c(age.rad, age.radg)

#Stage, making variable uniform across datasets
stage.v <- do.call(rbind, lapply(complete.lung.results, function (x) x[["Overall.Stage"]]))
stage.v <- case_when(stage.v == 1 ~ "I",
                     stage.v == 2 ~ "II",
                     stage.v == 3 ~ "IIIa",
                     stage.v == 4 ~ "IIIb",
                     TRUE ~ as.character(stage.v))
  


tab <- as.data.frame(cbind(surv.v, dead.alive.v, age.v, area.v, pixelcount.v, stage.v))

names(tab) <- c("surv.v", "dead.alive.v", "age.v", "area.v", "pixelcount.v", "stage.v")

tab$stage.v <- as.factor(tab$stage.v)
tab$surv.v <- as.numeric(as.character(tab$surv.v))
tab$dead.alive.v <- as.numeric(as.character(tab$dead.alive.v))
tab$age.v <- as.numeric(as.character(tab$age.v))
tab$pixelcount.v <- as.numeric(as.character(tab$pixelcount.v))
tab$area.v <- as.numeric(as.character(tab$area.v))

tab$stage.v <- factor(tab$stage.v, levels = c("I", "II", "IIIa", "IIIb", "IV", "0"))




library("survival")
library("survminer")
res.cox <- coxph(Surv(surv.v, dead.alive.v) ~ area.v + stage.v + age.v + pixelcount.v, data = tab)
summary(res.cox)

(res.zph1 <- cox.zph(res.cox))
ggcoxzph(res.zph1)


tab.alt <- tab
tab.alt$area.category <- ifelse(tab$area.v >= 1.2, "high", "low")

fit <- survfit(Surv(surv.v, dead.alive.v) ~ area.category, data = tab.alt)
ggsurvplot(fit, 
           ggtheme = theme_minimal(), data = tab.alt)

tab <- as.data.frame(cbind(surv.v, area.v, age.v, stage.v, dead.alive.v))
names(tab) <- c("surv.v", "area.v", "age.v", "stage.v", "dead.alive.v")

write.csv(tab, "tdasurvival.csv")





