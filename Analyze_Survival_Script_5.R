#Script Number 5 (Final Script)
#Reads in consolidated data
#Analyzes Survival by discrete and continuos method
library(plyr)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(reshape)
library(moments)
library(gridExtra)
library(rstatix)


#Reading in Results
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

#Looking at total feature curve by the survival quartiles
surv <- clinic.splitter.cust(complete.lung.results, "Tot.Feat.Count", "surv.cat")

surv25 <- surv[[1]]
surv50 <- surv[[2]]
surv75 <- surv[[3]]
surv100 <- surv[[4]]

#Function to median feature curves and provide mean and standard error
#Reason median is computed instead of mean is due to exponential distribution of data
feat.medianer <- function(feat.list) {
  arr <- array(unlist(feat.list), dim = c(dim(feat.list[[1]]), length(feat.list)))
  feat.list.av <- apply(arr, c(1,2), median)
  return(feat.list.av)
}

#Computing the average of each survival group
surv25.med <- feat.medianer(surv25)
surv50.med <- feat.medianer(surv50)
surv75.med <- feat.medianer(surv75)
surv100.med <- feat.medianer(surv100)

#Visualizing the survival groups
comb.data <- cbind(data.frame(rbind(surv25.med, surv50.med, surv75.med, surv100.med)), 
                   c(rep("surv25", dim(surv25.med)[1]),
                     rep("surv50", dim(surv50.med)[1]),
                     rep("surv75", dim(surv75.med)[1]),
                     rep("surv100", dim(surv100.med)[1])))
colnames(comb.data) <- c("filtration", "features", "survival_group")
comb.data$survival_group <- factor(comb.data$survival_group, 
                                   levels = c("surv25", "surv50", "surv75", "surv100"))

ggplot(comb.data, aes(filtration, features)) +
  geom_point() + facet_wrap(~survival_group) + theme_bw()

write.csv(comb.data, "./Results/discretesurvgraphs.csv")

#Collecting the moments of each group
surv25mom <- cbind(as.data.frame(t(sapply(surv25, function (q) all.moments(q[,2], 4)))), "surv25")
colnames(surv25mom) <- c("mom0", "mom1", "mom2", "mom3", "mom4", "surv.group")
surv50mom <- cbind(as.data.frame(t(sapply(surv50, function (q) all.moments(q[,2], 4)))), "surv50")
colnames(surv50mom) <- c("mom0", "mom1", "mom2", "mom3", "mom4", "surv.group")
surv75mom <- cbind(as.data.frame(t(sapply(surv75, function (q) all.moments(q[,2], 4)))), "surv75")
colnames(surv75mom) <- c("mom0", "mom1", "mom2", "mom3", "mom4", "surv.group")
surv100mom <- cbind(as.data.frame(t(sapply(surv100, function (q) all.moments(q[,2], 4)))), "surv100")
colnames(surv100mom) <- c("mom0", "mom1", "mom2", "mom3", "mom4", "surv.group")


surv.all.mom <- rbind(surv25mom, surv50mom, surv75mom, surv100mom)
write.csv(surv.all.mom, "./Results/allsurvmoments.csv")


#Creating one super table with all data
all.surv.moments.melt <- melt(surv.all.mom)
colnames(all.surv.moments.melt) <- c("surv.group", "moment", "feature.count")
all.surv.moments.melt.rem.zero <- subset(all.surv.moments.melt, moment != "mom0")
all.surv.moments.melt.rem.zero$surv.group <- factor(all.surv.moments.melt.rem.zero$surv.group, 
                                                    levels = c("surv25", "surv50", 
                                                               "surv75", "surv100"))


#####Analyzing moment 1 Across Survival Groups#####
mom1tab <- cbind(quantile(surv25mom[,2]), quantile(surv50mom[,2]), 
                 quantile(surv75mom[,2]), quantile(surv100mom[,2]))
colnames(mom1tab) <- c("surv25", "surv50", "surv75", "surv100")
mom1tab <- cbind(quantile(surv25mom[,2]), quantile(surv50mom[,2]), 
                 quantile(surv75mom[,2]), quantile(surv100mom[,2]))

#Comparing moment 1 distributions among survival groups
#Kruskal Wallis
mom1.kw <- all.surv.moments.melt.rem.zero %>%
  subset(moment == "mom1") %>%
  kruskal_test(feature.count ~ surv.group)

#Dunn's Test
mom1.dt <- all.surv.moments.melt.rem.zero %>%
  subset(moment == "mom1") %>%
  dunn_test(feature.count ~ surv.group, p.adjust.method = "bonferroni")

mom1comp <- rbind(
cbind(as.matrix(mom1.kw[, -c(2, 4, 6)]), NA, NA, NA),
c("y", "group1", "group2", "stat", "p", "p.adj"),
as.matrix(mom1.dt[, -c(4, 5, 9)])
)


#Writing results to csv
write.csv(mom1tab, "./Results/mom1table.csv")
write.csv(mom1comp, "./Results/mom1comp.csv")


#####Analyzing moment 2 Across Survival Groups#####
mom2tab <- cbind(quantile(surv25mom[,3]), quantile(surv50mom[,3]), 
                 quantile(surv75mom[,3]), quantile(surv100mom[,3]))
colnames(mom2tab) <- c("surv25", "surv50", "surv75", "surv100")
mom2tab <- cbind(quantile(surv25mom[,3]), quantile(surv50mom[,3]), 
                 quantile(surv75mom[,3]), quantile(surv100mom[,3]))

#Comparing moment 2 distributions among survival groups
#Kruskal Wallis
mom2.kw <- all.surv.moments.melt.rem.zero %>%
  subset(moment == "mom2") %>%
  kruskal_test(feature.count ~ surv.group)

#Dunn's Test
mom2.dt <- all.surv.moments.melt.rem.zero %>%
  subset(moment == "mom2") %>%
  dunn_test(feature.count ~ surv.group, p.adjust.method = "bonferroni")

mom2comp <- rbind(
  cbind(as.matrix(mom2.kw[, -c(2, 4, 6)]), NA, NA, NA),
  c("y", "group1", "group2", "stat", "p", "p.adj"),
  as.matrix(mom2.dt[, -c(4, 5, 9)])
)



#Writing results to csv
write.csv(mom2tab, "./Results/mom2table.csv")
write.csv(mom2comp, "./Results/mom2comp.csv")


#####Analyzing moment 3 Across Survival Groups#####
mom3tab <- cbind(quantile(surv25mom[,4]), quantile(surv50mom[,4]), 
                 quantile(surv75mom[,4]), quantile(surv100mom[,4]))
colnames(mom3tab) <- c("surv25", "surv50", "surv75", "surv100")
mom3tab <- cbind(quantile(surv25mom[,4]), quantile(surv50mom[,4]), 
                 quantile(surv75mom[,4]), quantile(surv100mom[,4]))

#Comparing moment 3 distributions among survival groups
#Kruskal Wallis
mom3.kw <- all.surv.moments.melt.rem.zero %>%
  subset(moment == "mom3") %>%
  kruskal_test(feature.count ~ surv.group)

#Dunn's Test
mom3.dt <- all.surv.moments.melt.rem.zero %>%
  subset(moment == "mom3") %>%
  dunn_test(feature.count ~ surv.group, p.adjust.method = "bonferroni")

mom3comp <- rbind(
  cbind(as.matrix(mom3.kw[, -c(2, 4, 6)]), NA, NA, NA),
  c("y", "group1", "group2", "stat", "p", "p.adj"),
  as.matrix(mom3.dt[, -c(4, 5, 9)])
)

#Writing results to csv
write.csv(mom3tab, "./Results/mom3table.csv")
write.csv(mom3comp, "./Results/mom3comp.csv")









#####Analyzing moment 4 Across Survival Groups#####
mom4tab <- cbind(quantile(surv25mom[,5]), quantile(surv50mom[,5]), 
                 quantile(surv75mom[,5]), quantile(surv100mom[,5]))
colnames(mom4tab) <- c("surv25", "surv50", "surv75", "surv100")
mom4tab <- cbind(quantile(surv25mom[,5]), quantile(surv50mom[,5]), 
                 quantile(surv75mom[,5]), quantile(surv100mom[,5]))

#Comparing moment 3 distributions among survival groups
#Kruskal Wallis
mom4.kw <- all.surv.moments.melt.rem.zero %>%
  subset(moment == "mom4") %>%
  kruskal_test(feature.count ~ surv.group)

#Dunn's Test
mom4.dt <- all.surv.moments.melt.rem.zero %>%
  subset(moment == "mom4") %>%
  dunn_test(feature.count ~ surv.group, p.adjust.method = "bonferroni")

mom4comp <- rbind(
  cbind(as.matrix(mom4.kw[, -c(2, 4, 6)]), NA, NA, NA),
  c("y", "group1", "group2", "stat", "p", "p.adj"),
  as.matrix(mom4.dt[, -c(4, 5, 9)])
)

#Writing results to csv
write.csv(mom4tab, "./Results/mom4table.csv")
write.csv(mom4comp, "./Results/mom4comp.csv")












####Multivariate Cox Analysis####
#Function to normalize a vector of values to 0 to 1
normalize <- function(col) {
  return ((col - min(col)) / (max(col) - min(col)))
}

#Survival Time
surv <- do.call(rbind, lapply(complete.lung.results, function (x) x[["Survival.time"]]))

#Dead vs Alive Status
dead.alive <- do.call(rbind, lapply(complete.lung.results, function (x) x[["deadstatus.event"]]))

#Gathering all the moments of the data
mom <- do.call(rbind, lapply(complete.lung.results, function (q) 
  all.moments(q[[1]][,2], order.max = 4)))

mom1 <- mom[,2]
mom2 <- mom[,3]
mom3 <- mom[,4]
mom4 <- mom[,5]

#Pixel Count 
pixelcount <- do.call(rbind, lapply(segment.slices.rds, function (x) 
  length(x)*dim(x[[1]])[1]*dim(x[[1]])[1]))


#Age
age.rad <- do.call(rbind, lapply(complete.lung.results[1:421], function (x) x[["age"]])) 
age.radg <- do.call(rbind, lapply(complete.lung.results[422:565], function (x) 
  x[["Age.at.Histological.Diagnosis"]])) 
age <- c(age.rad, age.radg)

#Stage, making variable uniform across datasets
stage <- do.call(rbind, lapply(complete.lung.results, function (x) x[["Overall.Stage"]]))
stage <- case_when(stage == 1 ~ "I",
                   stage == 2 ~ "II",
                   stage == 3 ~ "IIIa",
                   stage == 4 ~ "IIIb",
                     TRUE ~ as.character(stage))

#Adding Sex
sex.r <- do.call(rbind, lapply(complete.lung.results[1:421], function (x) as.character(x[["gender"]])))
sex.rg <- do.call(rbind, lapply(complete.lung.results[422:565], function (x) as.character(x[["Gender"]])))
sex <- c(sex.r, sex.rg)
sex <- case_when(sex == "male" ~ "Male",
                 sex == "female" ~ "Female", 
                 TRUE ~ sex)



#Consolidating into single table
tab <- as.data.frame(cbind(surv, dead.alive, mom1, mom2, mom3, 
                           mom4, pixelcount, age, stage, sex))

names(tab) <- c("surv", "dead.alive", "mom1", "mom2", "mom3", "mom4", 
                "pixelcount", "age", "stage", "sex")

tab$stage <- factor(tab$stage, levels = c("I", "II", "IIIa", "IIIb", "IV", "0"))
tab$surv <- as.numeric(as.character(tab$surv))
tab$dead.alive <- as.numeric(as.character(tab$dead.alive))
tab$age <- as.numeric(as.character(tab$age))
tab$pixelcount <- as.numeric(as.character(tab$pixelcount))
tab$mom1 <- as.numeric(as.character(tab$mom1))
tab$mom2 <- as.numeric(as.character(tab$mom2))
tab$mom3 <- as.numeric(as.character(tab$mom3))
tab$mom4 <- as.numeric(as.character(tab$mom4))


#Removing Stage 0 for table 1
write.csv(tab[-which(tab$stage == "0"),], "./Results/tab1data.csv")


#Scaling Values 0 to 50...more comparable with age
#Scale could be normalized to anything, 
#Just changes the hazard ratio value but does not impact significance
#Or meaning of results
tab$mom1 <- normalize(mom[,2]) * 50
tab$mom2 <- normalize(mom[,3]) * 50
tab$mom3 <- normalize(mom[,4]) * 50
tab$mom4 <- normalize(mom[,5]) * 50
#Scaling Values 0 to 50...more comparable with age and moment
tab$pixelcount <- normalize(tab$pixelcount) * 50




#n = 542,  22 of the 565 patients did not have age information, 
#1 patient had no stage information
sapply(colnames(tab), function (x) sum(is.na(tab[[x]])))

#MV Cox model
res.cox <- coxph(Surv(surv, dead.alive) ~ mom1 + mom2 + mom3 + mom4 + 
                   age + pixelcount + stage + sex, data = tab)
summary(res.cox)

#Retaining tab as tab.supplement for the supplemental cox where we don't exclude 
#Stage 0
tab.supp <- tab

#Redoing the model without stage 0 patients, who cause an infinite fit for stage 0 patients
which(tab$stage == "0")
#Result is no different but prevents the erraneous 0 to inf HR for stage 0 patients
#n= 536
tab.diff <- tab[-which(tab$stage == "0"), ]

#Dropping the stage 0 level
tab.diff$stage <- droplevels(tab.diff$stage)



res.cox <- coxph(Surv(surv, dead.alive) ~ mom1 + mom2 + mom3 + mom4 + 
                   age + pixelcount + stage + sex, data = tab.diff)
summary(res.cox)

#Getting stats for manuscript
summary(coxph(Surv(surv, dead.alive) ~ pixelcount, data = tab.diff))

#SHoenfeld resideual graph
schoenfeld <- ggcoxzph(cox.zph(res.cox), font.main = 6,font.submain = 6,
                       font.caption = 6,font.x = 6,font.y = 6,
                       font.tickslab = 6,font.legend = 6)

schoenfeld

ggsave("./Figures/schoenfeld.png", 
       arrangeGrob(grobs = schoenfeld), scale = 1, 
       width = 8, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)

#Correlation plot between moment 1 and tumor image size, Supplement
cor <- round(cor(tab.diff$mom1, tab.diff$pixelcount)^2, 3)
corgplot <- ggplot(data=tab.diff, aes(x=mom1, y=pixelcount)) +
  geom_point(size=.5, shape=21) +
  annotate(geom="text", x=3, y=54, label=paste("R^2=", cor), color="black") +
  theme_bw() +   
  labs(title = "Image Size vs Moment 1", x = "Scaled Moment 1 of all Feature Curves", 
       y = "Scaled Tumor Image Sizes of All Feature Curves") + 
  theme(plot.title = element_text(hjust = 0.5))
corgplot
  
ggsave("./Figures/cormom1size.png", plot = corgplot,
       scale = 1, width = 8, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)




#Creating easy sumdmary data frame for forest plot
MV.Cox.Ob <- data.frame(cbind(exp(coefficients(res.cox)), exp(confint(res.cox))))

#Adding pvalues
MV.Cox.Ob <- cbind(MV.Cox.Ob, summary(res.cox)$coefficients[, 5])

#Removing stage 0 var, as that factor was not analyzed in tab.diff
#Stage 0 result shown in supplement
MV.Cox.Ob <- na.omit(MV.Cox.Ob)
colnames(MV.Cox.Ob) <- c("HR", "LL", "UL", "pval")
labs <- rownames(MV.Cox.Ob)
MV.Cox.Ob$label <- labs


ggplot(data=MV.Cox.Ob, aes(x=label, y=HR, ymin=LL, ymax=UL)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Label") + ylab("Mean (95% CI)") 
  theme_bw()  # use a white background

write.csv(MV.Cox.Ob, "./Results/MVCoxModel.csv")

####Univartiate Cox Analysis####
vars <- c("mom1", "mom2", "mom3", "mom4", "age", "pixelcount", "stage", "sex")


#Creating easy sumdmary data frame for forest plot
uv.list.holder <- list()
for (i in 1:length(vars)) {
  form <- as.formula(paste("Surv(surv, dead.alive) ~ ", vars[i], sep = ""))
  cox.uv <- coxph(form, data = tab.diff)
  UV.Cox.row <- data.frame(cbind(exp(coefficients(cox.uv)), 
                                   exp(confint(cox.uv)),
                                   summary(cox.uv)$coefficients[, 5]))
  uv.list.holder[[i]] <- UV.Cox.row
  
}

UV.Cox.Ob <- do.call(rbind, uv.list.holder)
#Removing erraneous stage 0 var
UV.Cox.Ob <- na.omit(UV.Cox.Ob)
colnames(UV.Cox.Ob) <- c("HR", "LL", "UL", "pval")
labs <- rownames(UV.Cox.Ob)
UV.Cox.Ob$label <- labs


ggplot(data=UV.Cox.Ob, aes(x=label, y=HR, ymin=LL, ymax=UL)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Label") + ylab("Mean (95% CI)") 
theme_bw()  # use a white background

write.csv(UV.Cox.Ob, "./Results/UVCoxModel.csv")


####KM Curves Analysis####
tab.KM <- tab.diff[, 1:3]
quart.mom <- quantile(tab.KM[,3])
quart.mom

#Adding A Moment 1 quartile Variable
tab.KM[,4] <- ifelse(tab.KM[, 3] <= quart.mom[2], "Mom25",
              ifelse(tab.KM[, 3] <= quart.mom[3], "Mom50",
                     ifelse(tab.KM[, 3] <= quart.mom[4], "Mom75",
                            ifelse(tab.KM[, 3] <= quart.mom[5], "Mom100", NA))))
colnames(tab.KM) <- c("surv", "dead.alive", "mom1", "quart")

tab.KM$quart <- factor(tab.KM$quart, levels = c("Mom25", "Mom50", "Mom75", "Mom100"))

fit <- survfit(Surv(surv, dead.alive) ~ quart, data = tab.KM)



ggsurvplot(fit, data = tab.KM, conf.int = F)


write.csv(tab.KM, "./Results/KMCurves.csv")



####Supplemental Cox Keeping Stage 0####
#Redoing full Cox but with tab.supp, which keeps the stage 0 patients
res.cox.sup <- coxph(Surv(surv, dead.alive) ~ mom1 + mom2 + mom3 + mom4 + 
                   age + pixelcount + stage + sex, data = tab.supp)
summary(res.cox.sup)

#Creating easy sumdmary data frame for forest plot
sup.MV.Cox.Ob <- data.frame(cbind(exp(coefficients(res.cox.sup)), 
                                  exp(confint(res.cox.sup))))

#Adding pvalues
sup.MV.Cox.Ob <- cbind(sup.MV.Cox.Ob, summary(res.cox.sup)$coefficients[, 5])

#Removing stage 0 var, as that factor was not analyzed in tab.diff
#Stage 0 result shown in supplement
colnames(sup.MV.Cox.Ob) <- c("HR", "LL", "UL", "pval")
labs <- rownames(sup.MV.Cox.Ob)
sup.MV.Cox.Ob$label <- labs

write.csv(sup.MV.Cox.Ob, "./Results/supMVCoxModel.csv")

#Univartiate Sup Cox Analysis
vars.sup <- c("mom1", "mom2", "mom3", "mom4", "age", "pixelcount", "stage", "sex")


#Creating easy sumdmary data frame for forest plot
uv.list.holder.sup <- list()
for (i in 1:length(vars.sup)) {
  form <- as.formula(paste("Surv(surv, dead.alive) ~ ", vars.sup[i], sep = ""))
  cox.uv <- coxph(form, data = tab.supp)
  UV.Cox.row <- data.frame(cbind(exp(coefficients(cox.uv)), 
                                 exp(confint(cox.uv)),
                                 summary(cox.uv)$coefficients[, 5]))
  uv.list.holder.sup[[i]] <- UV.Cox.row
  
}

sup.UV.Cox.Ob <- do.call(rbind, uv.list.holder.sup)

colnames(sup.UV.Cox.Ob) <- c("HR", "LL", "UL", "pval")
labs <- rownames(sup.UV.Cox.Ob)
sup.UV.Cox.Ob$label <- labs


write.csv(sup.UV.Cox.Ob, "./Results/supUVCoxModel.csv")





