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
                        "Stage", "Cohort")

tad.pats$`Vital Status` <- case_when(tad.pats$`Vital Status` == "1" ~ "Dead",
                                     tad.pats$`Vital Status` == "0" ~ "Alive")

tad.pats$`Vital Status` <- as.factor(tad.pats$`Vital Status`)
tad.pats$`Cohort` <- as.factor(tad.pats$`Cohort`)


tabone <- CreateTableOne(vars = c("Vital Status", "Moment 1", "Moment 2", 
                        "Moment 3", "Moment 4", "Tumor Image Size", "Age", 
                        "Stage"),
               strata = "Cohort", 
               factorVars = c("Stage", "Vital Status"), 
               addOverall = TRUE, data = tad.pats)

table.mat <- print(tabone, nonnormal = c("Moment 1", "Moment 2", 
                        "Moment 3", "Moment 4", "Tumor Image Size"), missing = TRUE)

table.mat.df <- as.data.frame(table.mat)
table.mat.df$label <- c("Total Sample Size", "Vital Status (% Dead)", "Moment 1",
                        "Moment 2", "Moment 3", "Moment 4", "Tumor Image Size", "Age",
                        "Stage", "Stage I", "Stage II", "Stage IIIa", "Stage IIIb", 
                        "Stage IV")

#Removing the test column
table.mat.df <- table.mat.df[,-5]

colnames(table.mat.df) <- c("Whole Cohort", "NSCLC-Radiogenomics", "NSCLC-Radiomics", 
                            "p-value", "Proportion Missing", "label")

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
                     domain = NULL, bins = 2^(0:58), pretty = TRUE))

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



box.whisk <- ggplot(data = all.surv.moments.melt.rem.zero, aes(x = surv.group, y = log(feature.count))) +
  geom_boxplot(outlier.shape = 18) + facet_wrap(~moment, scale = "free") + 
  scale_y_continuous(expand = expand_scale(mult = c(.03, .15), add = c(0, 0))) +
  theme_bw() +   
  labs(title = "Box Plot of Feature Curve Moments", x = "Survival Group", y = "log(Moment Value)") + 
  scale_x_discrete(labels = c("Survival\nGroup 1", "Survival\nGroup 2", "Survival\nGroup 3", "Survival\nGroup 4")) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 5))

#Adding the non parametric stats
sig.comps <- list(c("Survival Group 1", "Survival Group 4"), 
                  c("Survival Group 1", "Survival Group 3"), 
                  c("Survival Group 2", "Survival Group 4"))
box.whisk <- box.whisk + stat_compare_means(comparisons = sig.comps, size = 2) + 
            stat_compare_means(label.y.npc = .01, label.x = .8, size = 2)

box.whisk

ggsave("./Figures/box.whisk.png", plot = box.whisk,
       scale = 1, width = 8, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)









####Supplement: Stat Compare Group Moment####
mom1.comp.data <- cbind(read.csv("./Results/mom1comp.csv"), "mom1")
colnames(mom1.comp.data) <- c("Survival Group Comparison", "Wilcoxon Statistic", "p-value", "Moment")
mom2.comp.data <- cbind(read.csv("./Results/mom2comp.csv"), "mom2")
colnames(mom2.comp.data) <- c("Survival Group Comparison", "Wilcoxon Statistic", "p-value", "Moment")
mom3.comp.data <- cbind(read.csv("./Results/mom3comp.csv"), "mom3")
colnames(mom3.comp.data) <- c("Survival Group Comparison", "Wilcoxon Statistic", "p-value", "Moment")
mom4.comp.data <- cbind(read.csv("./Results/mom4comp.csv"), "mom4")
colnames(mom4.comp.data) <- c("Survival Group Comparison", "Wilcoxon Statistic", "p-value", "Moment")
mom.tot.comp.data <- rbind(mom1.comp.data, mom2.comp.data, mom3.comp.data, mom4.comp.data)

mom.tot.comp.data$`Survival Group Comparison` <- revalue(mom.tot.comp.data$`Survival Group Comparison`, 
                                                     c("surv25 vs surv100"="Survival Group 1 vs Survival Group 4", 
                                                       "surv25 vs surv75"="Survival Group 1 vs Survival Group 3", 
                                                       "surv25 vs surv50"="Survival Group 1 vs Survival Group 2", 
                                                       "surv50 vs surv100"="Survival Group 2 vs Survival Group 4", 
                                                       "surv50 vs surv75"="Survival Group 2 vs Survival Group 3", 
                                                       "surv75 vs surv100"="Survival Group 3 vs Survival Group 4"))


mom.tot.comp.data$Moment <- revalue(mom.tot.comp.data$Moment, 
                                                 c("mom1"="Moment 1", 
                                                   "mom2"="Moment 2", 
                                                   "mom3"="Moment 3", 
                                                   "mom4"="Moment 4"))

mom.tot.comp.tab <- mom.tot.comp.data %>%
  dplyr::group_by(Moment) %>%
  gt(rowname_col = "Percentile Value") %>%
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
      columns = vars(`p-value`),
      rows = `p-value` < .05)
  )
 
mom.tot.comp.tab

gtsave(mom.tot.comp.tab, "./Figures/moments_surv_groups_statcomp.tex")



####Continuous Analysis
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
                 revalue(., c("mom1"="Moment 1", 
                              "mom2"="Moment 2",
                              "mom3"="Moment 3",
                              "mom4"="Moment 4",
                              "age" = "Age",
                              "pixelcount" = "Tumor Image Size",
                              "stageII" = "Stage II vs I",
                              "stageIIIa" = "Stage IIIa vs I",
                              "stageIIIb" = "Stage IIIb vs I",
                              "stageIV" = "Stage IV vs I"))

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
  scale_y_continuous(trans='log10', limits = c(.3,3)) + 
  geom_hline(yintercept=1, lty=2, size = .2) +  
  coord_flip() +
theme_bw() +   
  labs(title = "Survival Forest Plot", x = "Cox Variable", 
       y = "log(Hazard Ratio) (95% CI)") + 
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_x_discrete(labels = rev(c(expression(bold("Stage IIIb vs I")),
                              expression(bold("Stage IIIa vs I")),
                              "Stage II vs I",
                              expression(italic("Tumor Image Size")),
                              "Moment 4",
                              "Moment 3",
                              "Moment 2",
                              expression(bold("Moment 1")),
                              expression(bold("Age")))))

ggsave("./Figures/coxforest.png", plot = coxforest,
       scale = 1, width = 8, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)
                              

####KM Curve####
library(ggfortify)

KMcurve <- read.csv("./Results/KMCurves.csv")[,-1]

KMcurve$quart <- revalue(KMcurve$quart, c("Mom100" = "First Moment 75-100 percentile", 
                                           "Mom25" = "First Moment 0-25 percentile", 
                                           "Mom50" = "First Moment 25-50 percentile", 
                                            "Mom75" = "First Moment 50-75 percentile"))

KMcurve$quart <- factor(KMcurve$quart, levels = c("First Moment 0-25 percentile", 
                                                  "First Moment 25-50 percentile", 
                                                  "First Moment 50-75 percentile",
                                                  "First Moment 75-100 percentile"))

fit <- survfit(Surv(surv, dead.alive) ~ quart, data = KMcurve)

surv_pvalue(fit)

med.surv <- surv_median(fit)
med.surv

med.surv$strata <- revalue(med.surv$strata, c("quart=First Moment 0-25 percentile" = "First Moment 0-25 percentile", 
                                              "quart=First Moment 25-50 percentile" = "First Moment 25-50 percentile", 
                                              "quart=First Moment 50-75 percentile" = "First Moment 50-75 percentile", 
                                              "quart=First Moment 75-100 percentile" = "First Moment 75-100 percentile"))

colnames(med.surv) <- c("strata", "median", "upperlim", "lowerlim")


kmcurve <- autoplot(fit, data = tab.KM, conf.int = F, censor.shape = "|", censor.size = 2) + 
  geom_vline(aes(xintercept = median, color = strata), data = med.surv) +
  annotate("rect", xmin=med.surv$lowerlim[1], xmax=med.surv$upperlim[1], 
           ymin = 0, ymax = Inf, fill = "#00AFBB", alpha = 0.2) +
  annotate("rect", xmin=med.surv$lowerlim[2], xmax=med.surv$upperlim[2], 
           ymin = 0, ymax = Inf, fill = "#E7B800", alpha = 0.2) +
  annotate("rect", xmin=med.surv$lowerlim[3], xmax=med.surv$upperlim[3], 
         ymin = 0, ymax = Inf, fill = "#FC4E07", alpha = 0.2) +
  annotate("rect", xmin=med.surv$lowerlim[4], xmax=med.surv$upperlim[4], 
         ymin = 0, ymax = Inf, fill = "#9932CC", alpha = 0.2) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#9932CC")) + theme_classic() + 
  theme(legend.justification=c(1,1), legend.position=c(1,1),
        plot.title = element_text(hjust = 0.5)) +   
  labs(title = "KM Curves for Feature Curve Moment 1 Groups", x = "Survival Time (days)", 
       y = "Survival Probability", color = "Feature Curve Moment 1 Quartiles") + 
  guides(colour = guide_legend(override.aes = list(shape = 15)))
  


ggsave("./Figures/kmcurve.png", plot = kmcurve,
       scale = 1, width = 8, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)
















####Supplemental Cox Model with Stage 0####
sup.MV.Cox <- cbind(read.csv("./Results/supMVCoxModel.csv")[,-1], "MV")
colnames(sup.MV.Cox) <- c("HR", "LL", "UL", "pval", "label", "model")
sup.UV.Cox <- cbind(read.csv("./Results/supUVCoxModel.csv")[,-1], "UV")
colnames(sup.UV.Cox) <- c("HR", "LL", "UL", "pval", "label", "model")

sup.Cox.tot <- rbind(sup.UV.Cox, sup.MV.Cox)
colnames(sup.Cox.tot) <- c("Hazard Ratio", "Lower Bound", "Upper Bound", 
                       "p-value", "label", "model")
sup.Cox.tot$label <- sup.Cox.tot$label %>%
  revalue(., c("mom1"="Moment 1", 
               "mom2"="Moment 2",
               "mom3"="Moment 3",
               "mom4"="Moment 4",
               "age" = "Age",
               "pixelcount" = "Tumor Image Size",
               "stageII" = "Stage II vs I",
               "stageIIIa" = "Stage IIIa vs I",
               "stageIIIb" = "Stage IIIb vs I",
               "stageIV" = "Stage IV vs I",
               "stage0" = "Stage 0 vs I"))

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


