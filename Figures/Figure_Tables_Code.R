library(gt)
library(reshape)
library(ggplot2)
library(paletteer)
library(plyr)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(moments)
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



####Background Cubical Complex Figure####
#Reading in the slices for radiomics
segment.slices <- readRDS("radiomics.segments.rds")
#Reading in the complete lung results
complete.lung.results <- readRDS("complete_lung_results.rds")
#Reading in original phom
phom.mat <- readRDS("formatted_phom_matrix.rds")


#Selecting the first image
img <- segment.slices[[1]]

#13th slice shows a good view of tumor
slice <- img[[13]]

#Viewing actual slice
tum <- image(slice, col=grey(0:2041/2041), axes=FALSE, xlab="", ylab="", mar = c(0,0,0,0))

#Choosing 4 random HU values to create the binary images 
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
      xlab="", col.lab = "#FFDF20", cex.lab=1.5, mgp=c(.2,1,0))
image(slice.filt3, col=grey(0:2041/2041), axes=FALSE, ylab="HU: 200", 
      xlab="", col.lab = "#FF0808", cex.lab=1.5, mgp=c(.2,1,0))
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
  geom_vline(aes(xintercept = filt1norm), color = "#08FFF0") + 
  geom_vline(aes(xintercept = filt2norm), color = "#FFDF20") + 
  geom_vline(aes(xintercept = filt3norm), color = "#FF0808") + 
  theme(legend.position=c(.8, .5)) +
  theme(axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 9),
        axis.line.y = element_line(size = 3, colour = "white", linetype=2),
        plot.margin = unit(c(5.5, 5.5, 5.5, 20.5), "pt")) +
  scale_x_continuous(breaks = c(0, .2, .4, .6, .8, 1), limits = c(0,1), expand = c(0,0)) +
  scale_y_continuous(position = "left", expand = c(0,0)) 

barcode

#Getting the topological feature curve representation
topfeatcurv <- complete.lung.results[[1]][["Tot.Feat.Count"]]
colnames(topfeatcurv) <- c("filtration", "features")
curvplot <- ggplot(topfeatcurv, aes(filtration, features)) +
  geom_point() + theme_bw() + 
  labs(title = NULL, x = "Filtration (Normalized HU)", 
       y = "Total Topological Feature Count") + 
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6)) + 
  geom_vline(aes(xintercept = filt1norm), color = "#08FFF0") + 
  geom_vline(aes(xintercept = filt2norm), color = "#FFDF20") + 
  geom_vline(aes(xintercept = filt3norm), color = "#FF0808") +
  scale_x_continuous(breaks = c(0, .2, .4, .6, .8, 1), expand = c(0,0))

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
tad.pats <- read.csv("./Results/tab1data.csv")
tad.pats$cohort <- NA

#First 421 patients are in rad, remaining are in radg
#The 6 stage 0 patients removed were all from radg
tad.pats$cohort[1:421] <- "NSCLC-Radiomics"
tad.pats$cohort[422:559] <- "NSCLC-Radiogenomics"

#Removing erraneous first column
tad.pats <- tad.pats[,-1]

colnames(tad.pats) <- c("Survival", "Vital Status", "Moment 1", "Moment 2", 
                        "Moment 3", "Moment 4", "Tumor Image Size", "Age", 
                        "Stage", "Sex", "Cohort")

tad.pats$`Vital Status` <- case_when(tad.pats$`Vital Status` == "1" ~ "Dead",
                                     tad.pats$`Vital Status` == "0" ~ "Alive")

tad.pats$`Vital Status` <- as.factor(tad.pats$`Vital Status`)
tad.pats$`Cohort` <- as.factor(tad.pats$`Cohort`)

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
                        "Stage IV", "Sex (% Male)")

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



####Discrete Analysis
####Survival Group Feature Curv####
feat.curves.data <- read.csv("./Results/discretesurvgraphs.csv")
feat.curves.data$survival_group <- factor(feat.curves.data$survival_group,
                                          levels = c("surv25", "surv50",
                                                     "surv75", "surv100"))

feat.curves.data$survival_group <- revalue(feat.curves.data$survival_group, 
                                          c("surv25"="Survival Group 1", 
                                            "surv50"="Survival Group 2", 
                                            "surv75"="Survival Group 3", 
                                           "surv100"="Survival Group 4"))


surv.curve.plot <- ggplot(feat.curves.data, aes(filtration, features)) +
  geom_point() + facet_wrap(~survival_group) + theme_bw() + 
  labs(title = "Median Feature Curves For Survival Groups", x = "Filtration (Normalized HU)", 
       y = "Median Feature Count") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 5))


surv.curve.plot


ggsave("./Figures/surv.curve.plot.png", plot = surv.curve.plot,
       scale = 1, width = 8, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)

####Survival Key####
key.table <- cbind(c("Survival Group 1", "Survival Group 2", "Survival Group 3", "Survival Group 4"), 
      c("Less than 303 Days", "303 to 699 Days", "699 Days to 1634 Days", "More than 1634 Days"))
key.table <- data.frame(key.table)
colnames(key.table) <- c("Survival Group", "Survival Length")


table.key <- key.table %>%
  gt() %>%
  tab_header(
    title = "Length of Survival in Each Group",
  ) %>%
  cols_align("left")
table.key

gtsave(table.key, "./Figures/tablekey.tex")

####Supplement: Table of Moments For Each Sur Group####
mom1.tab.data <- cbind(read.csv("./Results/mom1table.csv"), "mom1")
colnames(mom1.tab.data) <- c("Percentile Value", "Survival Group 1", "Survival Group 2", 
                             "Survival Group 3", "Survival Group 4", "moment")
mom2.tab.data <- cbind(read.csv("./Results/mom2table.csv"), "mom2")
colnames(mom2.tab.data) <- c("Percentile Value", "Survival Group 1", "Survival Group 2", 
                             "Survival Group 3", "Survival Group 4", "moment")
mom3.tab.data <- cbind(read.csv("./Results/mom3table.csv"), "mom3")
colnames(mom3.tab.data) <- c("Percentile Value", "Survival Group 1", "Survival Group 2", 
                             "Survival Group 3", "Survival Group 4", "moment")
mom4.tab.data <- cbind(read.csv("./Results/mom4table.csv"), "mom4")
colnames(mom4.tab.data) <- c("Percentile Value", "Survival Group 1", "Survival Group 2", 
                             "Survival Group 3", "Survival Group 4", "moment")

mom.tabs <- rbind(mom1.tab.data, mom2.tab.data, mom3.tab.data, mom4.tab.data)
mom.tabs$moment <- revalue(mom.tabs$moment, c("mom1"="Moment 1", 
                                              "mom2"="Moment 2", 
                                              "mom3"="Moment 3", 
                                              "mom4"="Moment 4"))


mom.tabs$`Survival Group 1` <- trunc(mom.tabs$`Survival Group 1`) 
mom.tabs$`Survival Group 2` <- trunc(mom.tabs$`Survival Group 2`) 
mom.tabs$`Survival Group 3` <- trunc(mom.tabs$`Survival Group 3`) 
mom.tabs$`Survival Group 4` <- trunc(mom.tabs$`Survival Group 4`) 

mom.group.table <- mom.tabs %>%
dplyr::group_by(moment) %>%
gt(rowname_col = "Percentile Value") %>%
  tab_header(
    title = "Survival Groups Moment Distribution",
  ) %>%
  fmt_scientific(columns = 2:5,
                 rows = `Survival Group 1` >= 1000 |
                   `Survival Group 2` >= 1000 |
                   `Survival Group 3` >= 1000 |
                   `Survival Group 4` >= 1000,
                 decimals = 1) %>%   
  data_color(
             columns = 2:5,
             colors = scales::col_bin(
             palette = paletteer::paletteer_d(
             palette = "ggsci::red_material") %>% as.character(),
                     domain = NULL, bins = c(0,2^(0:58)), pretty = TRUE))
mom.group.table

gtsave(mom.group.table, "./Figures/moments_surv_groups.png")


####Box and Whisker of Moments for Groups####
all.surv.moments <- read.csv("./Results/allsurvmoments.csv")
all.surv.moments.melt <- melt(all.surv.moments)
colnames(all.surv.moments.melt) <- c("ID", "surv.group", "moment", "feature.count")
all.surv.moments.melt.rem.zero <- subset(all.surv.moments.melt, moment != "mom0")
all.surv.moments.melt.rem.zero$surv.group <- factor(all.surv.moments.melt.rem.zero$surv.group, 
                                                         levels = c("surv25", "surv50", 
                                                                    "surv75", "surv100"))
all.surv.moments.melt.rem.zero$surv.group <- revalue(all.surv.moments.melt.rem.zero$surv.group, 
                                                     c("surv25"="Survival Group 1", 
                                                       "surv50"="Survival Group 2", 
                                                       "surv75"="Survival Group 3", 
                                                       "surv100"="Survival Group 4"))

all.surv.moments.melt.rem.zero$moment <- revalue(all.surv.moments.melt.rem.zero$moment, 
                                                     c("mom1"="Moment 1", 
                                                       "mom2"="Moment 2", 
                                                       "mom3"="Moment 3", 
                                                       "mom4"="Moment 4"))


#Adding in ymin and ymax for ggplot
all.surv.moments.melt.rem.zero <- all.surv.moments.melt.rem.zero %>% mutate(ymin = 0) %>% 
                                   mutate(ymax = case_when(
                                                 moment == "Moment 1" ~ 15,
                                                 moment == "Moment 2" ~ 30,
                                                 moment == "Moment 3" ~ 40,
                                                 moment == "Moment 4" ~ 90))



box.whisk <- ggplot(data = all.surv.moments.melt.rem.zero, aes(x = surv.group, y = feature.count)) +
  geom_boxplot(outlier.shape = 18) + facet_wrap(~moment, scale = "free") + 
  scale_y_continuous(expand = expand_scale(mult = c(.03, .15), add = c(0, 0)), trans = "log10") +
  theme_bw() +   
  labs(title = "Box Plot of Feature Curve Moments", x = "Survival Group", y = "Moment Value") + 
  scale_x_discrete(labels = c("Survival\nGroup 1", "Survival\nGroup 2", "Survival\nGroup 3", "Survival\nGroup 4")) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 5))



#Adding the non parametric stats
stat.test <- all.surv.moments.melt.rem.zero %>%
  group_by(moment) %>%
  dunn_test(feature.count ~ surv.group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() 
stat.test$p.adj <- signif(stat.test$p.adj, 2)


box.whisk <- box.whisk + stat_pvalue_manual(stat.test, label = "p.adj", 
                                            y.position = c(4, 8.5, 13.3, 18), 
                                            size = 2, hide.ns = TRUE) + 
            stat_compare_means(label.y.npc = .01, label.x = .8, size = 2)
box.whisk


ggsave("./Figures/box.whisk.png", plot = box.whisk,
       scale = 1, width = 8, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)











####Supplement: Stat Compare Group Moment####

mom1.comp.data <- cbind(read.csv("./Results/mom1comp.csv"), "mom1")
mom2.comp.data <- cbind(read.csv("./Results/mom2comp.csv"), "mom2")
mom3.comp.data <- cbind(read.csv("./Results/mom3comp.csv"), "mom3")
mom4.comp.data <- cbind(read.csv("./Results/mom4comp.csv"), "mom4")

df.comp <- data.frame(`Group 1` = rep(NA, 24), `Group 2` = rep(NA, 24), 
                      `Dunn Statistic` = rep(NA, 24), `p-value` = rep(NA, 24), 
                      `adjusted p-value` = rep(NA, 24), `moment` = rep(NA, 24))

df.comp[1:6, ] <- as.matrix(mom1.comp.data[3:8, 3:8])
df.comp[7:12, ] <- as.matrix(mom2.comp.data[3:8, 3:8])
df.comp[13:18, ] <- as.matrix(mom3.comp.data[3:8, 3:8])
df.comp[19:24, ] <- as.matrix(mom4.comp.data[3:8, 3:8])


df.comp$moment <- revalue(df.comp$moment, c("mom1"= paste("Moment 1, KW H-statistic = ", 
                                                          round(as.numeric(as.character((mom1.comp.data$statistic[1]))), 2), ", (p = ", 
                                                          signif(as.numeric(as.character(mom1.comp.data$p[1])), 2), ")", sep = ""),
                                            "mom2"= paste("Moment 2, KW H-statistic = ", 
                                                          round(as.numeric(as.character((mom2.comp.data$statistic[1]))), 2), ", (p = ", 
                                                          signif(as.numeric(as.character(mom2.comp.data$p[1])), 2), ")", sep = ""),
                                            "mom3"= paste("Moment 3, KW H-statistic = ", 
                                                          round(as.numeric(as.character((mom3.comp.data$statistic[1]))), 2), ", (p = ", 
                                                          signif(as.numeric(as.character(mom3.comp.data$p[1])), 2), ")", sep = ""),
                                            "mom4"= paste("Moment 4, KW H-statistic = ", 
                                                          round(as.numeric(as.character((mom4.comp.data$statistic[1]))), 2), ", (p = ", 
                                                          signif(as.numeric(as.character(mom4.comp.data$p[1])), 2), ")", sep = "")))


colnames(df.comp) <- c("Group 1", "Group 2", "Dunn's z test statistic", "p-value", "adjusted p-value", 
                       "moment")


df.comp$`Dunn's z test statistic` <- round(as.numeric(as.character(
  df.comp$`Dunn's z test statistic`)), 2)
df.comp$`p-value` <- signif(as.numeric(as.character(
  df.comp$`p-value`)), 2)
df.comp$`adjusted p-value` <- signif(as.numeric(as.character(
  df.comp$`adjusted p-value`)), 2)

df.comp$`Group 1` <- revalue(df.comp$`Group 1`, c("surv25"= "Survival Group 1",
                                                  "surv50"= "Survival Group 2",
                                                  "surv75"= "Survival Group 3",
                                                  "surv100"= "Survival Group 4"))
df.comp$`Group 2` <- revalue(df.comp$`Group 2`, c("surv25"= "Survival Group 1",
                                                  "surv50"= "Survival Group 2",
                                                  "surv75"= "Survival Group 3",
                                                  "surv100"= "Survival Group 4"))


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
####Table of UV and MV Hazard Ratios####
MV.Cox <- cbind(read.csv("./Results/MVCoxModel.csv")[,-1], "MV")
colnames(MV.Cox) <- c("HR", "LL", "UL", "pval", "label", "model")
UV.Cox <- cbind(read.csv("./Results/UVCoxModel.csv")[,-1], "UV")
colnames(UV.Cox) <- c("HR", "LL", "UL", "pval", "label", "model")

Cox.tot <- rbind(UV.Cox, MV.Cox)
colnames(Cox.tot) <- c("Hazard Ratio", "Lower Bound", "Upper Bound", 
                       "p-value", "label", "model")
Cox.tot$label <- Cox.tot$label %>%
                 revalue(., c("mom1"="Scaled Moment 1", 
                              "mom2"="Scaled Moment 2",
                              "mom3"="Scaled Moment 3",
                              "mom4"="Scaled Moment 4",
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



####MV and UV Forest Plot####
Cox.tot.mod <- subset(Cox.tot, `p-value` <= .05)
Cox.tot.mod <- subset(Cox.tot, `label` != "Stage IV vs I")


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
                              "Male vs Female",
                              expression(italic("Scaled Tumor Image Size")),
                              "Scaled Moment 4",
                              "Scaled Moment 3",
                              "Scaled Moment 2",
                              expression(bold("Scaled Moment 1")),
                              expression(bold("Age")))))
coxforest


ggsave("./Figures/coxforest.png", plot = coxforest,
       scale = 1, width = 8, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)
                              

####KM Curve####

KMcurve <- read.csv("./Results/KMCurves.csv")[,-1]

#Organizing the factor levels
KMcurve$quart <- revalue(KMcurve$quart, c("Mom100" = "First Moment 75-100 percentile", 
                                           "Mom25" = "First Moment 0-25 percentile", 
                                           "Mom50" = "First Moment 25-50 percentile", 
                                            "Mom75" = "First Moment 50-75 percentile"))
KMcurve$quart <- factor(KMcurve$quart, levels = c("First Moment 0-25 percentile", 
                                                  "First Moment 25-50 percentile", 
                                                  "First Moment 50-75 percentile",
                                                  "First Moment 75-100 percentile"))
#Getting the fit
fit <- survfit(Surv(surv, dead.alive) ~ quart, data = KMcurve)

#Getting the median survival
med.surv <- surv_median(fit)
med.surv$strata <- revalue(med.surv$strata, c("quart=First Moment 0-25 percentile" = "First Moment 0-25 percentile", 
                                              "quart=First Moment 25-50 percentile" = "First Moment 25-50 percentile", 
                                              "quart=First Moment 50-75 percentile" = "First Moment 50-75 percentile", 
                                              "quart=First Moment 75-100 percentile" = "First Moment 75-100 percentile"))
colnames(med.surv) <- c("strata", "median", "lowerlim", "upperlim")



#This log rank compares all groups
log.rank.all <- survdiff(Surv(surv, dead.alive) ~ quart, data = KMcurve)
pval.tot <- 1 - pchisq(log.rank.all$chisq, length(log.rank.all$n) - 1)

#Pairwise post hoc
KMcurve %>% subset(quart == "First Moment 0-25 percentile" | 
                     quart == "First Moment 75-100 percentile") %>% 
  survdiff(Surv(surv, dead.alive) ~ quart, data = .)

KMcurve %>% subset(quart == "First Moment 25-50 percentile" | 
                     quart == "First Moment 75-100 percentile") %>% 
  survdiff(Surv(surv, dead.alive) ~ quart, data = .)

#This log rank post hoc compares pairwise groups
statkm <- pairwise_survdiff(Surv(surv, dead.alive) ~ quart, 
                            data = KMcurve, p.adjust.method = "bonferroni")


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
                     labels = c(paste("0-25 percentile, ", med.surv$median[1], " days (", med.surv$lowerlim[1], "-", 
                                      med.surv$upperlim[1], ")", sep = ""),
                                paste("25-50 percentile, ", med.surv$median[2], " days (", med.surv$lowerlim[2], "-", 
                                      med.surv$upperlim[2], ")", sep = ""),
                                paste("50-75 percentile, ", med.surv$median[3], " days (", med.surv$lowerlim[3], "-", 
                                      med.surv$upperlim[3], ")", sep = ""),
                                paste("75-100 percentile, ", med.surv$median[4], " days (", med.surv$lowerlim[4], "-", 
                                      med.surv$upperlim[4], ")", sep = ""))) + 
  theme_classic() + 
  theme(legend.justification=c(1,1), legend.position=c(1,.88),
        plot.title = element_text(hjust = 0.5)) +   
  labs(title = "KM Curves for Feature Curve Moment 1 Groups", x = "Survival Time (days)", 
       y = "Survival Probability", color = "Feature Curve Moment 1 Quartiles") + 
  guides(colour = guide_legend(override.aes = list(shape = 15)))


  
kmcurve <- kmcurve + stat_pvalue_manual(stat.km.df, label = "p.adj", 
                             y.position = c(1, .85), 
                             size = 3, hide.ns = TRUE, bracket.nudge.y = .1, tip.length = 0) +
          annotate("text", x = 85, y = .1, size = 3,
                   label = paste("p = ", format(signif(pval.tot, 2), scientific = FALSE)))

kmcurve

ggsave("./Figures/kmcurve.png", plot = kmcurve,
       scale = 1, width = 8, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)








print.surv.diff(log.rank.all)






####Supplemental Cox Model with Stage 0####
sup.MV.Cox <- cbind(read.csv("./Results/supMVCoxModel.csv")[,-1], "MV")
colnames(sup.MV.Cox) <- c("HR", "LL", "UL", "pval", "label", "model")
sup.UV.Cox <- cbind(read.csv("./Results/supUVCoxModel.csv")[,-1], "UV")
colnames(sup.UV.Cox) <- c("HR", "LL", "UL", "pval", "label", "model")

sup.Cox.tot <- rbind(sup.UV.Cox, sup.MV.Cox)
colnames(sup.Cox.tot) <- c("Hazard Ratio", "Lower Bound", "Upper Bound", 
                       "p-value", "label", "model")
sup.Cox.tot$label <- sup.Cox.tot$label %>%
  revalue(., c("mom1"="Scaled Moment 1", 
               "mom2"="Scaled Moment 2",
               "mom3"="Scaled Moment 3",
               "mom4"="Scaled Moment 4",
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


