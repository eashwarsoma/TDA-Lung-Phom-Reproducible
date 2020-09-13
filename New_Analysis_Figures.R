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
library(gt)
library(paletteer)
library(ggpubr)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(tableone)
library(TDAstats)
library(grid)
library(gridExtra)
library(ggplotify)
library(png)
library(rstatix)
library(ggfortify)

####Reading in Data####
#Reading in the fully formatted results
complete.lung.results <- readRDS("complete_lung_results.rds")

#Reading in the scans
segment.slices.rds <- c(readRDS("radiomics.segments.rds"), 
                        readRDS("radiogenomics.segments.rds"))
#Reading in original phom
phom.mat <- readRDS("formatted_phom_matrix.rds")


#####Set of Helpful Functions####
#Function to split groups based on clinical vairable keeping a certain
#Phom representation
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

#Function to linearly scale something to 0 to 1
normalize <- function(col) {
  return ((col - min(col)) / (max(col) - min(col)))
}

#Calculates raw moment
raw.moment <- function (vec, degree) {
  len <- length(vec)
  val <- sum(vec^degree)
  
  return(val/len)
}

#Calculates centralized moment (by subtracting mean)
cent.moment <- function (vec, degree) {
  avg <- mean(vec)
  len <- length(vec)
  val <- sum((vec-avg)^degree)
  
  return(val/len)
}

#Calculates standardized moment(by dividing sd^n)
stand.moment <- function (vec, degree) {
  sd <- sd(vec)
  len <- length(vec)
  val <- (sum((vec)^degree)/(sd^degree))
  
  return(val/len)
}

#Calculates standardized and centralized moment
standcent.moment <- function (vec, degree) {
  avg <- mean(vec)
  sd <- sd(vec)
  len <- length(vec)
  val <- sum(((vec-avg)/sd)^degree)
  return(val/len)
}

#Function to return median topological count at each HU
feat.medianer <- function(feat.list) {
  arr <- array(unlist(feat.list), dim = c(dim(feat.list[[1]]), length(feat.list)))
  feat.list.av <- apply(arr, c(1,2), median)
  return(feat.list.av)
}

#Function to get raw moments for each survival category
moment.puller <- function(which.feat, results) {
  
  len <- length(results)
  
  
  surv <- clinic.splitter.cust(results, which.feat, "surv.cat")
  
  surv25 <- surv[[1]]
  surv50 <- surv[[2]]
  surv75 <- surv[[3]]
  surv100 <- surv[[4]]
  
  surv25m1 <- sapply(surv25, raw.moment, 1)
  surv25m2 <- sapply(surv25, raw.moment, 2)
  surv25m3 <- sapply(surv25, raw.moment, 3)
  surv25m4 <- sapply(surv25, raw.moment, 4)
  surv25mom <- cbind(surv25m1, surv25m2, surv25m3, surv25m4)
  colnames(surv25mom) <- c("mom1", "mom2", "mom3", "mom4")
  
  surv50m1 <- sapply(surv50, raw.moment, 1)
  surv50m2 <- sapply(surv50, raw.moment, 2)
  surv50m3 <- sapply(surv50, raw.moment, 3)
  surv50m4 <- sapply(surv50, raw.moment, 4)
  surv50mom <- cbind(surv50m1, surv50m2, surv50m3, surv50m4)
  colnames(surv50mom) <- c("mom1", "mom2", "mom3", "mom4")
  
  surv75m1 <- sapply(surv75, raw.moment, 1)
  surv75m2 <- sapply(surv75, raw.moment, 2)
  surv75m3 <- sapply(surv75, raw.moment, 3)
  surv75m4 <- sapply(surv75, raw.moment, 4)
  surv75mom <- cbind(surv75m1, surv75m2, surv75m3, surv75m4)
  colnames(surv75mom) <- c("mom1", "mom2", "mom3", "mom4")
  
  surv100m1 <- sapply(surv100, raw.moment, 1)
  surv100m2 <- sapply(surv100, raw.moment, 2)
  surv100m3 <- sapply(surv100, raw.moment, 3)
  surv100m4 <- sapply(surv100, raw.moment, 4)
  surv100mom <- cbind(surv100m1, surv100m2, surv100m3, surv100m4)
  colnames(surv100mom) <- c("mom1", "mom2", "mom3", "mom4")
  
  list <- list(surv25mom, surv50mom, surv75mom, surv100mom)
  
  names(list) <- c("surv25", "surv50", "surv75", "surv100")
  
  meltlist <- melt(list)
  
  names(meltlist) <- c("scan_name", "moment", "value", "group")
  
  return(meltlist)
  
  
}

####Survival Plots####
#Getting the Zero Dim Curves based on survival group
surv <- clinic.splitter.cust(complete.lung.results, "Zero.Feat.Count", "surv.cat")
  
surv25 <- surv[["surv25"]]
surv50 <- surv[["surv50"]]
surv75 <- surv[["surv75"]]
surv100 <- surv[["surv100"]]

#Computing the mean zero feature curve of each survival group
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
levels(comb.data$survival_group) <- c("<303 days", "303 to 699 days", 
                                      "699 to 1634 days", ">1634 days")

surv.curve.plot <- ggplot(comb.data, aes(filtration, features)) +
  geom_point() + 
  facet_wrap(~survival_group) +
  theme_bw() + 
  labs(title = "Median Feature Curves For Survival Groups", x = "Filtration (Normalized HU)", 
       y = "Median Feature Count") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 5))
surv.curve.plot


ggsave("./Figures/surv.curve.plot.png", plot = surv.curve.plot,
       scale = 1, width = 8, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)


#####Box Plot: Comparing Moments for Zero Featuress####
#Extracting the zero dimensional feature curves
zero.feat.mom <- moment.puller("Zero.Feat.Count", complete.lung.results)

zero.feat.mom$group <- factor(zero.feat.mom$group, 
                              levels = c("surv25", "surv50", "surv75", "surv100"))
levels(zero.feat.mom$group) <- c("<303 days", "303 to 699 days", 
                                 "699 to 1634 days", ">1634 days")
levels(zero.feat.mom$moment) <- c("First Moment (AUC)", "Second Moment", 
                                      "Third Moment", "Fourth Moment")


options(scipen = 1)
box.whisk <- ggplot(data = zero.feat.mom, aes(x = group, y = value)) +
  geom_boxplot(outlier.shape = 18, notch = TRUE) + facet_wrap(~moment, scale = "free") + 
  scale_y_continuous(expand = expansion(mult = c(.03, .07), add = c(0, 0)), trans = "log10") +
  theme_bw() +   
  labs(title = "Box Plot of 0D Feature Curve Moments", x = "Survival Group", y = "Moment Value") +
  scale_x_discrete(labels = c("<303\ndays", "303 to\n699 days", 
                              "699 to\n1634 days", ">1634\ndays")) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 8))

zero.feat.mom %>%
  group_by(moment) %>%
  kruskal_test(value ~ group) %>%
  add_significance() 

#Note, bonferroni should multiple by 6, not 24 since we're making pairwise
#Comparisons within each moment which we consider separate variables
#Therewhere, we modify the padj from the forula
stat.test <- zero.feat.mom %>%
  group_by(moment) %>%
  dunn_test(value ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() 
stat.test$p.adj <- signif(stat.test$p*6, 2)
stat.test$p.adj.signif <- ifelse(stat.test$p.adj <= .05, "*", "ns")
print(stat.test)



box.whisk.stat <- box.whisk + stat_pvalue_manual(stat.test, label = "p.adj", 
                                            y.position = c(2.75, 3.25, 3,
                                                           6.9, 8, 7.50,
                                                           11.3, 13.3, 12.3,
                                                           15.3, 17.3, 16.3), 
                                            size = 2, hide.ns = TRUE) + 
  stat_compare_means(label.y.npc = .01, label.x = .8, size = 2)

box.whisk.stat


ggsave("./Figures/box.whisk.png", plot = box.whisk.stat,
       scale = 1, width = 8, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)

options(scipen = 10)


####Supplement: Stat Compare Group Moment####
kw.tab <- zero.feat.mom %>%
  group_by(moment) %>%
  kruskal_test(value ~ group) %>%
  add_significance() 
kw.tab$statistic <- round(kw.tab$statistic, 2)
kw.tab$p <- signif(kw.tab$p, 2)



relab.mom <- c(paste("Moment 1, KW H-statistic = ", kw.tab$statistic[1], 
               ", p = ", kw.tab$p[1], sep = ""),
               paste("Moment 2, KW H-statistic = ", kw.tab$statistic[2], 
                     ", p = ", kw.tab$p[2], sep = ""),
               paste("Moment 3, KW H-statistic = ", kw.tab$statistic[3], 
                     ", p = ", kw.tab$p[3], sep = ""),
               paste("Moment 4, KW H-statistic = ", kw.tab$statistic[4], 
                     ", p = ", kw.tab$p[4], sep = ""))

df.comp <- stat.test[, c("moment", "group1", "group2", "statistic", "p.adj")]
df.comp$statistic <- round(df.comp$statistic, 2)
df.comp$p.adj <- signif(df.comp$p.adj, 2)
df.comp$p.adj <- ifelse(df.comp$p.adj >= 1, 1, df.comp$p.adj)
levels(df.comp$moment) <- relab.mom
colnames(df.comp) <- c("moment", "Group 1", "Group 2", 
                       "Dunn's z-statistic", "adjusted p-value")


mom.tot.comp.tab <- df.comp %>%
  dplyr::group_by(moment) %>%
  gt() %>%
  tab_header(
    title = "Survival Group Feature Curve\nMoment Statistical Comparison",
  ) %>%
  cols_align("left")  %>%
  tab_style(
    style = list(
      cell_fill(color = "#F9E3D6"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = vars(`adjusted p-value`),
      rows = `adjusted p-value` < .05)
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#F9E3D6"),
      cell_text(weight = "bold")
    ),
    locations = cells_row_groups()
  ) 

mom.tot.comp.tab

gtsave(mom.tot.comp.tab, "./Figures/moments_surv_groups_statcomp.tex")



####Contineous ANalysis


####Supplement: Table of Moments For Each Survival Group####
melted.table.mom <- moment.puller("Zero.Feat.Count", complete.lung.results)
list.table.mom <- split(melted.table.mom[,-4], melted.table.mom$group, drop=TRUE)

surv25momrange <- list.table.mom[["surv25"]] %>% group_by(moment) %>%  
                   summarise(quantile = scales::percent(c(0, 0.25, 0.5, 0.75, 1)),
                   value = quantile(value, c(0, 0.25, 0.5, 0.75, 1)))
colnames(surv25momrange) <- c("moment", "Percentile Value", "<303 days")

surv50momrange <- list.table.mom[["surv50"]] %>% group_by(moment) %>%  
  summarise(quantile = scales::percent(c(0, 0.25, 0.5, 0.75, 1)),
            value = quantile(value, c(0, 0.25, 0.5, 0.75, 1)))
colnames(surv50momrange) <- c("moment", "Percentile Value", "303 to 699 days")

surv75momrange <- list.table.mom[["surv75"]] %>% group_by(moment) %>%  
  summarise(quantile = scales::percent(c(0, 0.25, 0.5, 0.75, 1)),
            value = quantile(value, c(0, 0.25, 0.5, 0.75, 1)))
colnames(surv75momrange) <- c("moment", "Percentile Value", "699 to 1634 days")

surv100momrange <- list.table.mom[["surv100"]] %>% group_by(moment) %>%  
  summarise(quantile = scales::percent(c(0, 0.25, 0.5, 0.75, 1)),
            value = quantile(value, c(0, 0.25, 0.5, 0.75, 1)))
colnames(surv100momrange) <- c("moment", "Percentile Value", ">1634 days")

#Merging Dataframes
allmomrange <- merge(surv25momrange, surv50momrange)
allmomrange <- merge(allmomrange, surv75momrange)
allmomrange <- merge(allmomrange, surv100momrange)

allmomrange$moment <- revalue(allmomrange$moment, c("mom1"="Moment 1", 
                                              "mom2"="Moment 2", 
                                              "mom3"="Moment 3", 
                                              "mom4"="Moment 4"))


allmomrange$`<303 days` <- signif(allmomrange$`<303 days`, 3)
allmomrange$`303 to 699 days` <- signif(allmomrange$`303 to 699 days`, 3)
allmomrange$`699 to 1634 days` <- signif(allmomrange$`699 to 1634 days`, 3)
allmomrange$`>1634 days` <- signif(allmomrange$`>1634 days`, 3)

allmomrange$`Percentile Value` <- factor(allmomrange$`Percentile Value`, 
                                         levels = c("0%", "25%", "50%", "75%", "100%"))
allmomrange <- allmomrange[order(allmomrange$`Percentile Value`),]

mom.group.table <- allmomrange %>%
  dplyr::group_by(moment) %>%
  gt(rowname_col = "Percentile Value") %>%
  tab_header(
    title = "Survival Group Moment Distributions",
  ) %>%
  fmt_scientific(columns = 3:6,
                 rows = `<303 days` >= 1000 |
                   `303 to 699 days` >= 1000 |
                   `699 to 1634 days` >= 1000 |
                   `>1634 days` >= 1000,
                 decimals = 1) %>%   
  data_color(
    columns = 3:6,
    colors = scales::col_bin(
      palette = paletteer::paletteer_d(
        palette = "ggsci::red_material") %>% as.character(),
      domain = NULL, bins = c(0,2^(0:58)), pretty = TRUE))
mom.group.table

gtsave(mom.group.table, "./Figures/moments_surv_groups.png")

####Cox Model Data Table####
#Survival Time
surv <- do.call(rbind, lapply(complete.lung.results, function (x) x[["Survival.time"]]))

#Dead vs Alive Status
dead.alive <- do.call(rbind, lapply(complete.lung.results, function (x) x[["deadstatus.event"]]))

#feat0 moments
f0.mom1 <- sapply(complete.lung.results, function (x) raw.moment(x[["Zero.Feat.Count"]][,2], 1))
f0.mom2 <- sapply(complete.lung.results, function (x) raw.moment(x[["Zero.Feat.Count"]][,2], 2))
f0.mom3 <- sapply(complete.lung.results, function (x) raw.moment(x[["Zero.Feat.Count"]][,2], 3))
f0.mom4 <- sapply(complete.lung.results, function (x) raw.moment(x[["Zero.Feat.Count"]][,2], 4))

#Pixel Count 
pixelcount <- do.call(rbind, lapply(segment.slices.rds, function (x) 
  length(x)*dim(x[[1]])[1]*dim(x[[1]])[2]))

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
tab <- as.data.frame(cbind(surv, dead.alive, 
                           f0.mom1, f0.mom2, f0.mom3, f0.mom4, 
                           pixelcount, age, stage, sex))

names(tab) <- c("surv", "dead.alive", 
                "f0.mom1", "f0.mom2", "f0.mom3", "f0.mom4", 
                "pixelcount", "age", "stage", "sex")

tab$stage <- factor(tab$stage, levels = c("I", "II", "IIIa", "IIIb", "IV", "0"))
tab$surv <- as.numeric(as.character(tab$surv))
tab$dead.alive <- as.numeric(as.character(tab$dead.alive))
tab$age <- as.numeric(as.character(tab$age))
tab$pixelcount <- as.numeric(as.character(tab$pixelcount))
tab$f0.mom1 <- normalize(as.numeric(as.character(tab$f0.mom1))) * 50
tab$f0.mom2 <- normalize(as.numeric(as.character(tab$f0.mom2))) * 50
tab$f0.mom3 <- normalize(as.numeric(as.character(tab$f0.mom3))) * 50
tab$f0.mom4 <- normalize(as.numeric(as.character(tab$f0.mom4))) * 50


#Scaling Values 0 to 50...more comparable with age and moment
tab$pixelcount <- normalize(tab$pixelcount) * 50

#n = 542,  22 of the 565 patients did not have age information, 
#1 patient had no stage information
#Keeping tab.supp for supplemental model showing no convergence
sapply(colnames(tab), function (x) sum(is.na(tab[[x]])))
tab.supp <- tab
tab.diff <- tab[-which(tab$stage == "0"), ]

#Dropping the stage 0 level
tab.diff$stage <- droplevels(tab.diff$stage)


####Cox MV Analysis####
vars <- c("f0.mom1", "f0.mom2", "f0.mom3", "f0.mom4", 
          "pixelcount", "age", "stage", "sex")

form <- as.formula(paste("Surv(surv, dead.alive) ~", paste(vars, collapse = "+")))

res.cox <- coxph(form, data = tab.diff)

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
cor <- round(cor(tab.diff$f0.mom1, tab.diff$pixelcount)^2, 3)
corgplot <- ggplot(data=tab.diff, aes(x=f0.mom1, y=pixelcount)) +
  geom_point(size=.5, shape=21) +
  annotate(geom="text", x=3, y=54, label=paste("R^2=", cor), color="black") +
  theme_bw() +   
  labs(title = "Image Size vs Moment 1", x = "Scaled Moment 1 of all Zero Feature Curves", 
       y = "Scaled Tumor Image Sizes") + 
  theme(plot.title = element_text(hjust = 0.5))
corgplot

ggsave("./Figures/cormom1size.png", plot = corgplot,
       scale = 1, width = 8, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)


#Creating easy sumdmary data frame for forest plot
MV.Cox.Ob <- data.frame(cbind(exp(coefficients(res.cox)), exp(confint(res.cox))))

#Adding pvalues
MV.Cox.Ob <- cbind(MV.Cox.Ob, summary(res.cox)$coefficients[, 5])

colnames(MV.Cox.Ob) <- c("HR", "LL", "UL", "pval")
labs <- rownames(MV.Cox.Ob)
MV.Cox.Ob$label <- labs


####Univariate analysis####
vars <- c("f0.mom1", "f0.mom2", "f0.mom3", "f0.mom4", 
          "pixelcount", "age", "stage", "sex")
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
colnames(UV.Cox.Ob) <- c("HR", "LL", "UL", "pval")
labs <- rownames(UV.Cox.Ob)
UV.Cox.Ob$label <- labs



####Cox Tables####
MV.Cox <- cbind(MV.Cox.Ob, "MV")
UV.Cox <- cbind(UV.Cox.Ob, "UV")
colnames(MV.Cox) <- c("HR", "LL", "UL", "pval", "label", "model")
colnames(UV.Cox) <- c("HR", "LL", "UL", "pval", "label", "model")
Cox.tot <- rbind(UV.Cox, MV.Cox)

colnames(Cox.tot) <- c("Hazard Ratio", "Lower Bound", "Upper Bound", 
                       "p-value", "label", "model")

Cox.tot$label <- Cox.tot$label %>%
  revalue(., c("f0.mom1"="Scaled Moment 1", 
               "f0.mom2"="Scaled Moment 2",
               "f0.mom3"="Scaled Moment 3",
               "f0.mom4"="Scaled Moment 4",
               "age" = "Age",
               "pixelcount" = "Scaled Tumor Image Size",
               "stageII" = "Stage II vs I",
               "stageIIIa" = "Stage IIIa vs I",
               "stageIIIb" = "Stage IIIb vs I",
               "stageIV" = "Stage IV vs I",
               "sexMale" = "Male vs Female"))

Cox.tot$model <- Cox.tot$model %>%
  revalue(., c("UV"="Univariate Model", 
               "MV"="Multivariate Model"))

Cox.tot$`Hazard Ratio` <- round(Cox.tot$`Hazard Ratio`, 3)
Cox.tot$`Lower Bound` <- round(Cox.tot$`Lower Bound`, 3)
Cox.tot$`Upper Bound` <- round(Cox.tot$`Upper Bound`, 3)

Cox.tot$`p-value` <- signif(Cox.tot$`p-value`, digits = 2)

modHR <- paste(Cox.tot$`Hazard Ratio`, " (", Cox.tot$`Lower Bound`, "-", 
               Cox.tot$`Upper Bound`, ")", sep = "")

UVmodtab <- cbind(modHR, Cox.tot$`p-value`, as.character(Cox.tot$`label`), 
                  as.character(Cox.tot$`model`)) %>% subset(.[,4] == "Univariate Model")
UVmodtab <- UVmodtab[,-4]
UVmodtab <- as.data.frame(UVmodtab)
colnames(UVmodtab) <- c("Univariate Hazard Ratio", "p-value", "label")
MVmodtab <- cbind(modHR, Cox.tot$`p-value`, as.character(Cox.tot$`label`), 
                  as.character(Cox.tot$`model`)) %>% subset(.[,4] == "Multivariate Model")
MVmodtab <- MVmodtab[,-4]
colnames(MVmodtab) <- c("Multivariate Hazard Ratio", "p-value", "label")
MVmodtab <- as.data.frame(MVmodtab)

comb.mod.HR.tab <- merge(UVmodtab, MVmodtab, by = "label")
colnames(comb.mod.HR.tab) <- c("label", "Univariate Model HR", "p-valueu", 
                               "Multivariate Model HR", "p-valuem")

options(scipen = 10)
comb.mod.HR.tab$`p-valueu` <- as.numeric(as.character(comb.mod.HR.tab$`p-valueu`))
comb.mod.HR.tab$`p-valuem` <- as.numeric(as.character(comb.mod.HR.tab$`p-valuem`))


Cox.table <- comb.mod.HR.tab %>%
  gt(rowname_col = "label") %>%
  tab_header(
    title = "Cox Proportional Hazard Model",
  ) %>%
  cols_align("center") %>%
  cols_label(
    `p-valueu` = "p-value",
    `p-valuem` = "p-value"
  ) %>%
  fmt_scientific(columns = 3,
                 rows = `p-valueu` < .0001,
                 decimals = 1) %>%
  fmt_scientific(columns = 5,
                 rows = `p-valuem` < .0001,
                 decimals = 1)

Cox.table

gtsave(Cox.table, "./Figures/coxtable.tex")


####Cox Forest Plot####
Cox.tot.mod <- subset(Cox.tot, `label` != "Stage IV vs I")

Cox.tot.mod$label <- as.factor(Cox.tot.mod$label)
levels(Cox.tot.mod$label)


coxforest <- ggplot(data=Cox.tot.mod, aes(x=label, y=`Hazard Ratio`, ymin=`Lower Bound`, ymax=`Upper Bound`)) +
  # geom_pointrange(color='black', shape=19, size = .4, fatten = .01) + 
  geom_errorbar(width=0.1, size=.5) + 
  geom_point(size=.5, shape=18) +
  facet_wrap(~model) +
  scale_y_continuous(trans='log10', limits = c(.2,4)) + 
  geom_hline(yintercept=1, lty=2, size = .2) +  
  coord_flip() +
  theme_bw() +   
  labs(title = "Survival Forest Plot", x = "Cox Variable", 
       y = "log(Hazard Ratio) (95% CI)") + 
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_x_discrete(labels = rev(c(expression(bold("Stage IIIb vs I")),
                                  expression(bold("Stage IIIa vs I")),
                                  "Stage II vs I",
                                  expression(italic("Scaled Tumor Image Size")),
                                  "Scaled Moment 4",
                                  "Scaled Moment 3",
                                  "Scaled Moment 2",
                                  expression(bold("Scaled Moment 1")),
                                  "Male vs Female",
                                  expression(bold("Age")))))
coxforest

ggsave("./Figures/coxforest.png", plot = coxforest,
       scale = 1, width = 8, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)


####KM Curve####
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

#Organizing the factor levels
tab.KM$quart <- revalue(tab.KM$quart, c("Mom100" = "First Moment 75-100 percentile", 
                                          "Mom25" = "First Moment 0-25 percentile", 
                                          "Mom50" = "First Moment 25-50 percentile", 
                                          "Mom75" = "First Moment 50-75 percentile"))
tab.KM$quart <- factor(tab.KM$quart, levels = c("First Moment 0-25 percentile", 
                                                  "First Moment 25-50 percentile", 
                                                  "First Moment 50-75 percentile",
                                                  "First Moment 75-100 percentile"))
#Getting the fit
fit <- survfit(Surv(surv, dead.alive) ~ quart, data = tab.KM)

#Getting the median survival
med.surv <- surv_median(fit)
med.surv$strata <- revalue(med.surv$strata, c("quart=First Moment 0-25 percentile" = "First Moment 0-25 percentile", 
                                              "quart=First Moment 25-50 percentile" = "First Moment 25-50 percentile", 
                                              "quart=First Moment 50-75 percentile" = "First Moment 50-75 percentile", 
                                              "quart=First Moment 75-100 percentile" = "First Moment 75-100 percentile"))
colnames(med.surv) <- c("strata", "median", "lowerlim", "upperlim")

#This log rank compares all groups
log.rank.all <- survdiff(Surv(surv, dead.alive) ~ quart, data = tab.KM)
pval.tot <- 1 - pchisq(log.rank.all$chisq, length(log.rank.all$n) - 1)

#Pairwise post hoc
tab.KM %>% subset(quart == "First Moment 0-25 percentile" | 
                     quart == "First Moment 75-100 percentile") %>% 
  survdiff(Surv(surv, dead.alive) ~ quart, data = .)

tab.KM %>% subset(quart == "First Moment 25-50 percentile" | 
                     quart == "First Moment 75-100 percentile") %>% 
  survdiff(Surv(surv, dead.alive) ~ quart, data = .)

#This log rank post hoc compares pairwise groups
statkm <- pairwise_survdiff(Surv(surv, dead.alive) ~ quart, 
                            data = tab.KM, p.adjust.method = "bonferroni")

#Creating a stat object DF for ggplot
stat.km.df <- data.frame(group1 = rep(NA, 6), group2 = rep(NA, 6), 
                         p.adj = rep(NA, 6), p.adj.signif = rep(NA, 6))

stat.km.df$group1 = c(rep("First Moment 0-25 percentile", 3),
                      rep("First Moment 25-50 percentile", 2),
                      rep("First Moment 50-75 percentile", 1))
stat.km.df$group2 = c("First Moment 25-50 percentile",
                      "First Moment 50-75 percentile",
                      "First Moment 75-100 percentile",
                      "First Moment 50-75 percentile",
                      "First Moment 75-100 percentile",
                      "First Moment 75-100 percentile")

stat.km.df$p.adj[1:3] <- statkm$p.value[, "First Moment 0-25 percentile"]
stat.km.df$p.adj[4:5] <- na.omit(statkm$p.value[, "First Moment 25-50 percentile"])
stat.km.df$p.adj[6] <- na.omit(statkm$p.value[, "First Moment 50-75 percentile"])
stat.km.df$p.adj.signif <- ifelse(stat.km.df$p.adj < .05, "*", "ns")
#Converting group names to x positions
stat.km.df$group1 <- case_when(stat.km.df$group1 == med.surv$strata[1] ~ med.surv$median[1],
                               stat.km.df$group1 == med.surv$strata[2] ~ med.surv$median[2],
                               stat.km.df$group1 == med.surv$strata[3] ~ med.surv$median[3],
                               stat.km.df$group1 == med.surv$strata[4] ~ med.surv$median[4])
stat.km.df$group2 <- case_when(stat.km.df$group2 == med.surv$strata[1] ~ med.surv$median[1],
                               stat.km.df$group2 == med.surv$strata[2] ~ med.surv$median[2],
                               stat.km.df$group2 == med.surv$strata[3] ~ med.surv$median[3],
                               stat.km.df$group2 == med.surv$strata[4] ~ med.surv$median[4])


stat.km.df$p.adj <- signif(stat.km.df$p.adj, 2)


kmcurve <- autoplot(fit, data = tab.KM, conf.int = F, censor.shape = "|", 
                    censor.size = 2) + 
  geom_vline(aes(xintercept = median, color = strata), alpha = .5, data = med.surv) +
  scale_y_continuous(expand = c(0, .03), breaks = c(0, .2, .4, .6, .8, 1)) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#9932CC"),
                     labels = c(paste("Q1 = 0 - 0.33, ", med.surv$median[1], " days (", med.surv$lowerlim[1], "-", 
                                      med.surv$upperlim[1], ")", sep = ""),
                                paste("Q2 = 0.33 - 1.46, ", med.surv$median[2], " days (", med.surv$lowerlim[2], "-", 
                                      med.surv$upperlim[2], ")", sep = ""),
                                paste("Q3 = 1.46 - 4.27, ", med.surv$median[3], " days (", med.surv$lowerlim[3], "-", 
                                      med.surv$upperlim[3], ")", sep = ""),
                                paste("Q4 = 4.27 - 50, ", med.surv$median[4], " days (", med.surv$lowerlim[4], "-", 
                                      med.surv$upperlim[4], ")", sep = ""))) + 
  theme_classic() + 
  theme(legend.justification=c(1,1), legend.position=c(1,.88),
        plot.title = element_text(hjust = 0.5)) +   
  labs(title = "KM Curves for Zero Feature Curve Scaled Moment 1 Groups", x = "Survival Time (days)", 
       y = "Survival Probability", color = "Zero Feature Curve\nScaled Moment 1 Quartiles") + 
  guides(colour = guide_legend(override.aes = list(shape = 15)))
kmcurve


kmcurve.stat <- kmcurve + stat_pvalue_manual(stat.km.df, label = "p.adj", 
                                        y.position = c(1, .85), 
                                        size = 3, hide.ns = TRUE, bracket.nudge.y = .1, tip.length = 0) +
  annotate("text", x = 85, y = .1, size = 3,
           label = paste("p = ", format(signif(pval.tot, 2), scientific = FALSE)))

kmcurve.stat

ggsave("./Figures/kmcurve.png", plot = kmcurve.stat,
       scale = 1, width = 8, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)





####Background Cubical Complex Figure####
#Selecting the first image
img <- segment.slices.rds[[1]]

#13th slice shows a good view of tumor and some calcifications
slice <- img[[13]]

#Viewing actual slice
tum <- image(slice, col=grey(0:2041/2041), axes=FALSE, xlab="", ylab="", mar = c(0,0,0,0))

#Viewing all the images
par(mfrow=c(5,5), mar = c(0, 0, 0, 0))
for (i in 1:length(img)) {
  image(img[[i]], col=grey(0:2041/2041), axes=FALSE, xlab="", ylab="", mar = c(0,0,0,0))
}
dev.off()

#Viewing all the images, but binary
img.bin <- lapply(img, function (x) ifelse(x >= -34, 1, 0))
par(mfrow=c(5,5), mar = c(0, 0, 0, 0))
for (i in 1:length(img.bin)) {
  image(img.bin[[i]], col=grey(0:2041/2041), axes=FALSE, xlab="", ylab="", mar = c(0,0,0,0))
}
dev.off()


#Choosing 3 HU values illustrate the binary images concept well
HUfilt1 <- -900
HUfilt2 <- -34
HUfilt3 <- 200

slice.filt1 <- ifelse(slice >= HUfilt1, 1, 0)

slice.filt2 <- ifelse(slice >= HUfilt2, 1, 0)

slice.filt3 <- ifelse(slice >= HUfilt3, 1, 0)

#Saving as png file
png(file = "./Figures/bitmaps.png", width = 3, height = 8, units = "in", res = 800)
par(mfrow=c(4,1), mar = c(1,1.8,0.4,0))
image(slice, col=grey(0:2041/2041), axes=FALSE, ylab="Slice 13", 
      xlab="", col.lab = "black", cex.lab=1.5, mgp=c(.2,1,0))
image(slice.filt1, col=grey(0:2041/2041), axes=FALSE, ylab="HU: -900", 
      xlab="", col.lab = "#08FFF0", cex.lab=1.5, mgp=c(.2,1,0))
image(slice.filt2, col=grey(0:2041/2041), axes=FALSE, ylab="HU: -34", 
      xlab="", col.lab = "#C39E05", cex.lab=1.5, mgp=c(.2,1,0))
image(slice.filt3, col=grey(0:2041/2041), axes=FALSE, ylab="HU: 200", 
      xlab="", col.lab = "#9500DA", cex.lab=1.5, mgp=c(.2,1,0))
while (!is.null(dev.list()))  dev.off()



#From Script 4, the absolute min and max of the HU units across whole sets
#Getting the corresponding normalized HU values
abs.max <- 3071
abs.min <- -1342

filt1norm <- (HUfilt1 - abs.min)/(abs.max - abs.min)
filt2norm <- (HUfilt2 - abs.min)/(abs.max - abs.min)
filt3norm <- (HUfilt3 - abs.min)/(abs.max - abs.min)


#Getting the barcode diagram of this phom
phom.barcode <- phom.mat[[1]]
barcode <- plot_barcode(as.matrix(phom.barcode)) + labs(x = "Filtration (Normalized HU)",y = " ", 
                                                        color = "Dimension") + 
  geom_vline(aes(xintercept = filt1norm), color = "#08FFF0", size = 1) + 
  geom_vline(aes(xintercept = filt2norm), color = "#C39E05", size = 1) + 
  geom_vline(aes(xintercept = filt3norm), color = "#9500DA", size = 1) + 
  theme(legend.position=c(.8, .5)) +
  theme(axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 9),
        axis.line.y = element_line(size = 3, colour = "white", linetype=2),
        plot.margin = unit(c(5.5, 5.5, 5.5, 20.5), "pt")) +
  scale_x_continuous(breaks = c(0, .2, .4, .6, .8, 1), limits = c(0,1), expand = c(0,0)) +
  scale_y_continuous(position = "left", expand = c(0,0)) 

barcode

#Getting all the topological feature curve representation
topfeatcurv <- complete.lung.results[[1]][["Tot.Feat.Count"]]
topfeatcurv <- cbind(topfeatcurv, 
                     complete.lung.results[[1]][["Zero.Feat.Count"]][,2],
                     complete.lung.results[[1]][["One.Feat.Count"]][,2],
                     complete.lung.results[[1]][["Two.Feat.Count"]][,2])
colnames(topfeatcurv) <- c("filtration", "total features", 
                           "zero features", "one features", "two features")
topfeatcurv <- as.data.frame(topfeatcurv)
topfeatcurv.melt <- melt(topfeatcurv, id.vars = c("filtration"), 
     measure.vars = c("total features", "zero features", "one features", "two features"))
colnames(topfeatcurv.melt) <- c("filtration", "feature.type", "feature.count")
levels(topfeatcurv.melt$feature.type) <- c("Total Features", "0D Features", "1D Features", "2D Features")

group.colors <- c("#000000", "#f15f36", "#25af35", "#5388fb")

curvplot <- ggplot(topfeatcurv.melt, aes(x = filtration, y = feature.count, 
                                         color = feature.type, alpha = feature.type)) +
   geom_path(size = 1) + theme_bw() + 
  scale_color_manual(values = group.colors) +
  scale_alpha_manual(values = c(0.2, 1, 0.2, 0.2)) +
  labs(title = NULL, x = "Filtration (Normalized HU)", 
       y = "Topological Feature Curves", color = "Feature Dimension") + 
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = c(0.85, 0.75),
        legend.background=element_blank()) + 
  guides(alpha = "none") +
  geom_vline(aes(xintercept = filt1norm), color = "#08FFF0", size = 1) + 
  geom_vline(aes(xintercept = filt2norm), color = "#C39E05", size = 1) + 
  geom_vline(aes(xintercept = filt3norm), color = "#9500DA", size = 1) +
  scale_x_continuous(breaks = c(0, .2, .4, .6, .8, 1), expand = c(0,0)) + 
  guides(color = guide_legend(override.aes = list(alpha = c(0.2, 1, 0.2, 0.2))))
curvplot

PNG <- readPNG("./Figures/bitmaps.png")

#creating layout
lay <- rbind(c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,3,3),
             c(1,1,3,3))

png(file = "./Figures/cubcomp.png", width = 10, height = 8, units = "in", res = 800)
grid.arrange(rasterGrob(PNG, interpolate=TRUE, just = "left"), 
             as.grob(barcode), as.grob(curvplot), layout_matrix = lay)
dev.off()
#I then manually edited this image after to crop white space and add astericks



#### Flow Diagram####
prism <- grViz("digraph flowchart {
      # node definitions with substituted label text
      node [fontname = Helvetica, shape = rectangle]        
      
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']
      tab6 [label = '@@6']
      tab7 [label = '@@7']
      

      # edge definitions with the node IDs
      tab1 -> tab2 -> tab3;
      tab4 -> tab5 -> tab3;
      tab3 -> tab6 -> tab7;
      }

      [1]: 'TCIA NSCLC-Radiomics Dataset\\n n = 422 CT scans'
      [2]: 'Scans with Segmentation Data\\n n=421'
      [3]: 'Total Scans Analyzed n=565\\n Performed Discrete Analysis\\n Performed UV Cox Models*'
      [4]: 'TCIA NSCLC Radiogenomics Dataset\\n n = 211 CT scans'
      [5]: 'Scans with Segmentation Data\\n n=144'
      [6]: 'Excluded Missing age (n=22)\\n Excluded Missing stage (n=1) \\n Excluded Stage 0 scans (n=6)**'
      [7]: 'Total Scans Analyzed n=536\\n Performed MV Cox Models'
      ")

export_svg(prism) %>% charToRaw %>% 
  rsvg_png(., file = "./Figures/prismatic.png", width = 2400, height = 2000)





####Table One####
tad.pats <- tab.diff
tad.pats$cohort <- NA

#First 421 patients are in rad, remaining are in radg
#The 6 stage 0 patients removed were all from radg
tad.pats$cohort[1:421] <- "NSCLC-Radiomics"
tad.pats$cohort[422:559] <- "NSCLC-Radiogenomics"

#Removing first column
tad.pats <- tad.pats[,-1]

colnames(tad.pats) <- c("Vital Status", "Moment 1", "Moment 2", 
                        "Moment 3", "Moment 4", "Tumor Image Size", "Age", 
                        "Stage", "Sex", "Cohort")

tad.pats$`Vital Status` <- case_when(tad.pats$`Vital Status` == "1" ~ "Dead",
                                     tad.pats$`Vital Status` == "0" ~ "Alive")

tad.pats$`Vital Status` <- as.factor(tad.pats$`Vital Status`)
tad.pats$`Cohort` <- as.factor(tad.pats$`Cohort`)

tad.pats$Sex <- factor(tad.pats$Sex, levels=rev(levels(tad.pats$Sex))) 




#Printing some stats for manuscript
chisq_test(table(tad.pats$Cohort, tad.pats$Stage))
tad.pats %>% t_test(Age ~ Cohort, var.equal	= TRUE)
chisq_test(table(tad.pats$Cohort, tad.pats$Sex))
tad.pats %>% wilcox_test(`Moment 1` ~ Cohort)
tad.pats %>% wilcox_test(`Tumor Image Size` ~ Cohort)
chisq_test(table(tad.pats$Cohort, tad.pats$`Vital Status`))





tabone <- CreateTableOne(vars = c("Vital Status", "Moment 1", "Moment 2", 
                                  "Moment 3", "Moment 4", "Tumor Image Size", "Age", 
                                  "Stage", "Sex"),
                         strata = "Cohort", 
                         factorVars = c("Stage", "Vital Status", "Sex"), 
                         addOverall = TRUE, data = tad.pats)



table.mat <- print(tabone, nonnormal = c("Moment 1", "Moment 2", 
                                         "Moment 3", "Moment 4", "Tumor Image Size"), 
                   missing = TRUE)

table.mat.df <- as.data.frame(table.mat)
table.mat.df$label <- c("Total Sample Size", "Vital Status (% Dead)", "Moment 1",
                        "Moment 2", "Moment 3", "Moment 4", "Tumor Image Size", "Age",
                        "Stage", "Stage I", "Stage II", "Stage IIIa", "Stage IIIb", 
                        "Stage IV", "Sex (% Female)")

#Removing the test column
table.mat.df <- table.mat.df[,-5]


colnames(table.mat.df) <- c("Whole Cohort", "NSCLC-Radiogenomics", "NSCLC-Radiomics", 
                            "p-value", "Proportion Missing", "label")



#Modifying the really big numbers to be in sci-not format 
col1_25 <- tabone$ContTable$Overall[2:5,"p25"] %>% formatC(., format = "e", digits = 2)
col2_25 <- tabone$ContTable$`NSCLC-Radiogenomics`[2:5,"p25"] %>% formatC(., format = "e", digits = 2)
col3_25 <- tabone$ContTable$`NSCLC-Radiomics`[2:5,"p25"] %>% formatC(., format = "e", digits = 2)

col1_50 <- tabone$ContTable$Overall[2:5,"median"] %>% formatC(., format = "e", digits = 2)
col2_50 <- tabone$ContTable$`NSCLC-Radiogenomics`[2:5,"median"] %>% formatC(., format = "e", digits = 2)
col3_50 <- tabone$ContTable$`NSCLC-Radiomics`[2:5,"median"] %>% formatC(., format = "e", digits = 2)

col1_75 <- tabone$ContTable$Overall[2:5,"p75"] %>% formatC(., format = "e", digits = 2)
col2_75 <- tabone$ContTable$`NSCLC-Radiogenomics`[2:5,"p75"] %>% formatC(., format = "e", digits = 2)
col3_75 <- tabone$ContTable$`NSCLC-Radiomics`[2:5,"p75"] %>% formatC(., format = "e", digits = 2)

table.mat.df$`Whole Cohort` <- as.character(table.mat.df$`Whole Cohort`)
table.mat.df[4:7, 'Whole Cohort'] <- paste(col1_50, " [", col1_25, ", ", col1_75, "]", sep = "")
table.mat.df$`NSCLC-Radiogenomics` <- as.character(table.mat.df$`NSCLC-Radiogenomics`)
table.mat.df[4:7, 'NSCLC-Radiogenomics'] <- paste(col2_50, " [", col2_25, ", ", col2_75, "]", sep = "")
table.mat.df$`NSCLC-Radiomics` <- as.character(table.mat.df$`NSCLC-Radiomics`)
table.mat.df[4:7, 'NSCLC-Radiomics'] <- paste(col3_50, " [", col3_25, ", ", col3_75, "]", sep = "")

#Actual Sci-not in overleaf
table.mat.df.gt <- table.mat.df %>% gt(rowname_col = "label") 

table.mat.df.gt

gtsave(table.mat.df.gt, "./Figures/tableone.tex")



####Supplemental Cox Model with Stage 0####

#Redoing full Cox but with tab.supp, which keeps the stage 0 patients
res.cox.sup <- coxph(Surv(surv, dead.alive) ~ f0.mom1 + f0.mom2 + f0.mom3 + f0.mom4 + 
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

#Univartiate Sup Cox Analysis
vars.sup <- c("f0.mom1", "f0.mom2", "f0.mom3", "f0.mom4", 
               "pixelcount", "age", "stage", "sex")

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



sup.MV.Cox <- cbind(sup.MV.Cox.Ob, "MV")
sup.UV.Cox <- cbind(sup.UV.Cox.Ob, "UV")
colnames(sup.MV.Cox) <- c(colnames(sup.MV.Cox.Ob), "model")
colnames(sup.UV.Cox) <- c(colnames(sup.UV.Cox.Ob), "model")

sup.Cox.tot <- rbind(sup.UV.Cox, sup.MV.Cox)
colnames(sup.Cox.tot) <- c("Hazard Ratio", "Lower Bound", "Upper Bound", 
                           "p-value", "label", "model")

sup.Cox.tot$label <- sup.Cox.tot$label %>%
  revalue(., c("f0.mom1"="Scaled Moment 1", 
               "f0.mom2"="Scaled Moment 2",
               "f0.mom3"="Scaled Moment 3",
               "f0.mom4"="Scaled Moment 4",
               "age" = "Age",
               "pixelcount" = "Scaled Tumor Image Size",
               "stageII" = "Stage II vs I",
               "stageIIIa" = "Stage IIIa vs I",
               "stageIIIb" = "Stage IIIb vs I",
               "stageIV" = "Stage IV vs I",
               "stage0" = "Stage 0 vs I",
               "sexMale" = "Male vs Female"))

sup.Cox.tot$model <- sup.Cox.tot$model %>%
  revalue(., c("UV"="Univariate Model", 
               "MV"="Multivariate Model"))

sup.Cox.tot$`Hazard Ratio` <- round(sup.Cox.tot$`Hazard Ratio`, 3)
sup.Cox.tot$`Lower Bound` <- round(sup.Cox.tot$`Lower Bound`, 3)
sup.Cox.tot$`Upper Bound` <- round(sup.Cox.tot$`Upper Bound`, 3)

sup.Cox.tot$`p-value` <- signif(sup.Cox.tot$`p-value`, digits = 2)

sup.modHR <- paste(sup.Cox.tot$`Hazard Ratio`, " (", sup.Cox.tot$`Lower Bound`, "-", 
                   sup.Cox.tot$`Upper Bound`, ")", sep = "")

sup.UVmodtab <- cbind(sup.modHR, sup.Cox.tot$`p-value`, as.character(sup.Cox.tot$`label`), 
                      as.character(sup.Cox.tot$`model`)) %>% subset(.[,4] == "Univariate Model")
sup.UVmodtab <- sup.UVmodtab[,-4]
sup.UVmodtab <- as.data.frame(sup.UVmodtab)
colnames(sup.UVmodtab) <- c("Univariate Hazard Ratio", "p-value", "label")
sup.MVmodtab <- cbind(sup.modHR, sup.Cox.tot$`p-value`, as.character(sup.Cox.tot$`label`), 
                      as.character(sup.Cox.tot$`model`)) %>% subset(.[,4] == "Multivariate Model")
sup.MVmodtab <- sup.MVmodtab[,-4]
colnames(sup.MVmodtab) <- c("Multivariate Hazard Ratio", "p-value", "label")
sup.MVmodtab <- as.data.frame(sup.MVmodtab)

sup.comb.mod.HR.tab <- merge(sup.UVmodtab, sup.MVmodtab, by = "label")
colnames(sup.comb.mod.HR.tab) <- c("label", "Univariate Model HR", "p-valueu", 
                                   "Multivariate Model HR", "p-valuem")

options(scipen = 10)
sup.comb.mod.HR.tab$`p-valueu` <- as.numeric(as.character(sup.comb.mod.HR.tab$`p-valueu`))
sup.comb.mod.HR.tab$`p-valuem` <- as.numeric(as.character(sup.comb.mod.HR.tab$`p-valuem`))


sup.Cox.table <- sup.comb.mod.HR.tab %>%
  gt(rowname_col = "label") %>%
  tab_header(
    title = "Cox Proportional Hazard Model",
  ) %>%
  cols_align("center") %>%
  cols_label(
    `p-valueu` = "p-value",
    `p-valuem` = "p-value"
  ) %>%
  fmt_scientific(columns = 3,
                 rows = `p-valueu` < .0001,
                 decimals = 1) %>%
  fmt_scientific(columns = 5,
                 rows = `p-valuem` < .0001,
                 decimals = 1)

sup.Cox.table

gtsave(sup.Cox.table, "./Figures/supcoxtable.tex")

















##If you have any questions, please don't hesistate to email evs27@case.edu



