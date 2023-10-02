library("here")
library("ape")
library("phytools")
library("foreach")
library("doMC")
library("OUwie")
library("caper")
library("corHMM")
library("data.table")
library("expm")
library("doMC")
library("ggplot2")
library("MetBrewer")

i_am("R/houwie_results.R")

for(i in 1:length(list.files(here("output/houwie")))){
    assign(gsub("-", ".", gsub("-nohid.RDS", "", list.files(here("output/houwie"))[i])),
           readRDS(paste0("output/houwie/", list.files(here("output/houwie"))[i]))
           )
}

## Domatium Growth

### Corola Length

model.list.domgrow.corleng <- list(er.bm1 = domgrow.corleng.er[[1]],
     er.bms = domgrow.corleng.er[[2]],
     er.ou1 = domgrow.corleng.er[[3]],
     er.oum = domgrow.corleng.er[[4]],
     er.ouma = domgrow.corleng.er[[5]],
     er.oumv = domgrow.corleng.er[[6]],
     sym.bm1 = domgrow.corleng.sym[[1]],
     sym.bms = domgrow.corleng.sym[[2]],
     sym.ou1 = domgrow.corleng.sym[[3]],
     sym.oum = domgrow.corleng.sym[[4]],
     sym.ouma = domgrow.corleng.sym[[5]],
     sym.oumv = domgrow.corleng.sym[[6]],
     ard.bm1 = domgrow.corleng.ard[[1]],
     ard.bms = domgrow.corleng.ard[[2]],
     ard.ou1 = domgrow.corleng.ard[[3]],
     ard.oum = domgrow.corleng.ard[[4]],
     ard.ouma = domgrow.corleng.ard[[5]],
     ard.oumv = domgrow.corleng.ard[[6]])

pars.domgrow.corleng <- getModelAvgParams(model.list.domgrow.corleng, force = FALSE)

plot.data.domgrow.corleng <- melt(pars.domgrow.corleng)
plot.data.domgrow.corleng$discrete <- "Domatium Growth"
plot.data.domgrow.corleng$continuous <- "Corolla Length"
plot.data.domgrow.corleng$tip_state <- factor(c("Apical", "Diffuse")[as.integer(plot.data.domgrow.corleng$tip_state)], levels = c("Apical", "Diffuse"))

## ggplot(plot.data.domgrow.corleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Leaf Area

model.list.domgrow.leafarea <- list(er.bm1 = domgrow.leafarea.er[[1]],
     er.bms = domgrow.leafarea.er[[2]],
     er.ou1 = domgrow.leafarea.er[[3]],
     er.oum = domgrow.leafarea.er[[4]],
     er.ouma = domgrow.leafarea.er[[5]],
     er.oumv = domgrow.leafarea.er[[6]],
     sym.bm1 = domgrow.leafarea.sym[[1]],
     sym.bms = domgrow.leafarea.sym[[2]],
     sym.ou1 = domgrow.leafarea.sym[[3]],
     sym.oum = domgrow.leafarea.sym[[4]],
     sym.ouma = domgrow.leafarea.sym[[5]],
     sym.oumv = domgrow.leafarea.sym[[6]],
     ard.bm1 = domgrow.leafarea.ard[[1]],
     ard.bms = domgrow.leafarea.ard[[2]],
     ard.ou1 = domgrow.leafarea.ard[[3]],
     ard.oum = domgrow.leafarea.ard[[4]],
     ard.ouma = domgrow.leafarea.ard[[5]],
     ard.oumv = domgrow.leafarea.ard[[6]])

pars.domgrow.leafarea <- getModelAvgParams(model.list.domgrow.leafarea, force = FALSE)

plot.data.domgrow.leafarea <- melt(pars.domgrow.leafarea)
plot.data.domgrow.leafarea$discrete <- "Domatium Growth"
plot.data.domgrow.leafarea$continuous <- "Leaf Area"
plot.data.domgrow.leafarea$tip_state <- factor(c("Apical", "Diffuse")[as.integer(plot.data.domgrow.leafarea$tip_state)], levels = c("Apical", "Diffuse"))

## ggplot(plot.data.domgrow.leafarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Petiole Length

model.list.domgrow.petleng <- list(er.bm1 = domgrow.petleng.er[[1]],
     er.bms = domgrow.petleng.er[[2]],
     er.ou1 = domgrow.petleng.er[[3]],
     er.oum = domgrow.petleng.er[[4]],
     er.ouma = domgrow.petleng.er[[5]],
     er.oumv = domgrow.petleng.er[[6]],
     sym.bm1 = domgrow.petleng.sym[[1]],
     sym.bms = domgrow.petleng.sym[[2]],
     sym.ou1 = domgrow.petleng.sym[[3]],
     sym.oum = domgrow.petleng.sym[[4]],
     sym.ouma = domgrow.petleng.sym[[5]],
     sym.oumv = domgrow.petleng.sym[[6]],
     ard.bm1 = domgrow.petleng.ard[[1]],
     ard.bms = domgrow.petleng.ard[[2]],
     ard.ou1 = domgrow.petleng.ard[[3]],
     ard.oum = domgrow.petleng.ard[[4]],
     ard.ouma = domgrow.petleng.ard[[5]],
     ard.oumv = domgrow.petleng.ard[[6]])

pars.domgrow.petleng <- getModelAvgParams(model.list.domgrow.petleng, force = FALSE)

plot.data.domgrow.petleng <- melt(pars.domgrow.petleng)
plot.data.domgrow.petleng$discrete <- "Domatium Growth"
plot.data.domgrow.petleng$continuous <- "Petiole Length"
plot.data.domgrow.petleng$tip_state <- factor(c("Apical", "Diffuse")[as.integer(plot.data.domgrow.petleng$tip_state)], levels = c("Apical", "Diffuse"))

## ggplot(plot.data.domgrow.petleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")

### Stem Area

model.list.domgrow.stemarea <- list(er.bm1 = domgrow.stemarea.er[[1]],
     er.bms = domgrow.stemarea.er[[2]],
     er.ou1 = domgrow.stemarea.er[[3]],
     er.oum = domgrow.stemarea.er[[4]],
     er.ouma = domgrow.stemarea.er[[5]],
     er.oumv = domgrow.stemarea.er[[6]],
     sym.bm1 = domgrow.stemarea.sym[[1]],
     sym.bms = domgrow.stemarea.sym[[2]],
     sym.ou1 = domgrow.stemarea.sym[[3]],
     sym.oum = domgrow.stemarea.sym[[4]],
     sym.ouma = domgrow.stemarea.sym[[5]],
     sym.oumv = domgrow.stemarea.sym[[6]],
     ard.bm1 = domgrow.stemarea.ard[[1]],
     ard.bms = domgrow.stemarea.ard[[2]],
     ard.ou1 = domgrow.stemarea.ard[[3]],
     ard.oum = domgrow.stemarea.ard[[4]],
     ard.ouma = domgrow.stemarea.ard[[5]],
     ard.oumv = domgrow.stemarea.ard[[6]])

pars.domgrow.stemarea <- getModelAvgParams(model.list.domgrow.stemarea, force = FALSE)

plot.data.domgrow.stemarea <- melt(pars.domgrow.stemarea)
plot.data.domgrow.stemarea$discrete <- "Domatium Growth"
plot.data.domgrow.stemarea$continuous <- "Stem Area"
plot.data.domgrow.stemarea$tip_state <- factor(c("Apical", "Diffuse")[as.integer(plot.data.domgrow.stemarea$tip_state)], levels = c("Apical", "Diffuse"))

## ggplot(plot.data.domgrow.stemarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")






## Hole Diameter

### Corola Length

model.list.holediam.corleng <- list(er.bm1 = holediam.corleng.er[[1]],
     er.bms = holediam.corleng.er[[2]],
     er.ou1 = holediam.corleng.er[[3]],
     er.oum = holediam.corleng.er[[4]],
     er.ouma = holediam.corleng.er[[5]],
     er.oumv = holediam.corleng.er[[6]],
     sym.bm1 = holediam.corleng.sym[[1]],
     sym.bms = holediam.corleng.sym[[2]],
     sym.ou1 = holediam.corleng.sym[[3]],
     sym.oum = holediam.corleng.sym[[4]],
     sym.ouma = holediam.corleng.sym[[5]],
     sym.oumv = holediam.corleng.sym[[6]],
     ard.bm1 = holediam.corleng.ard[[1]],
     ard.bms = holediam.corleng.ard[[2]],
     ard.ou1 = holediam.corleng.ard[[3]],
     ard.oum = holediam.corleng.ard[[4]],
     ard.ouma = holediam.corleng.ard[[5]],
     ard.oumv = holediam.corleng.ard[[6]])

pars.holediam.corleng <- getModelAvgParams(model.list.holediam.corleng, force = FALSE)

plot.data.holediam.corleng <- melt(pars.holediam.corleng)
plot.data.holediam.corleng$discrete <- "Hole Diameter"
plot.data.holediam.corleng$continuous <- "Corolla Length"
plot.data.holediam.corleng$tip_state <- factor(c("Several Large\nat Base", "One Large\nat Base", "All Large")[as.integer(plot.data.holediam.corleng$tip_state)], levels = c("Several Large\nat Base", "One Large\nat Base", "All Large"))

## ggplot(plot.data.holediam.corleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Leaf Area

model.list.holediam.leafarea <- list(er.bm1 = holediam.leafarea.er[[1]],
     er.bms = holediam.leafarea.er[[2]],
     er.ou1 = holediam.leafarea.er[[3]],
     er.oum = holediam.leafarea.er[[4]],
     er.ouma = holediam.leafarea.er[[5]],
     er.oumv = holediam.leafarea.er[[6]],
     sym.bm1 = holediam.leafarea.sym[[1]],
     sym.bms = holediam.leafarea.sym[[2]],
     sym.ou1 = holediam.leafarea.sym[[3]],
     sym.oum = holediam.leafarea.sym[[4]],
     sym.ouma = holediam.leafarea.sym[[5]],
     sym.oumv = holediam.leafarea.sym[[6]],
     ard.bm1 = holediam.leafarea.ard[[1]],
     ard.bms = holediam.leafarea.ard[[2]],
     ard.ou1 = holediam.leafarea.ard[[3]],
     ard.oum = holediam.leafarea.ard[[4]],
     ard.ouma = holediam.leafarea.ard[[5]],
     ard.oumv = holediam.leafarea.ard[[6]])

pars.holediam.leafarea <- getModelAvgParams(model.list.holediam.leafarea, force = FALSE)

plot.data.holediam.leafarea <- melt(pars.holediam.leafarea)
plot.data.holediam.leafarea$discrete <- "Hole Diameter"
plot.data.holediam.leafarea$continuous <- "Leaf Area"
plot.data.holediam.leafarea$tip_state <- factor(c("Several Large\nat Base", "One Large\nat Base", "All Large")[as.integer(plot.data.holediam.leafarea$tip_state)], levels = c("Several Large\nat Base", "One Large\nat Base", "All Large"))

## ggplot(plot.data.holediam.leafarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Petiole Length

model.list.holediam.petleng <- list(er.bm1 = holediam.petleng.er[[1]],
     er.bms = holediam.petleng.er[[2]],
     er.ou1 = holediam.petleng.er[[3]],
     er.oum = holediam.petleng.er[[4]],
     er.ouma = holediam.petleng.er[[5]],
     er.oumv = holediam.petleng.er[[6]],
     sym.bm1 = holediam.petleng.sym[[1]],
     sym.bms = holediam.petleng.sym[[2]],
     sym.ou1 = holediam.petleng.sym[[3]],
     sym.oum = holediam.petleng.sym[[4]],
     sym.ouma = holediam.petleng.sym[[5]],
     sym.oumv = holediam.petleng.sym[[6]],
     ard.bm1 = holediam.petleng.ard[[1]],
     ard.bms = holediam.petleng.ard[[2]],
     ard.ou1 = holediam.petleng.ard[[3]],
     ard.oum = holediam.petleng.ard[[4]],
     ard.ouma = holediam.petleng.ard[[5]],
     ard.oumv = holediam.petleng.ard[[6]])

pars.holediam.petleng <- getModelAvgParams(model.list.holediam.petleng, force = FALSE)

plot.data.holediam.petleng <- melt(pars.holediam.petleng)
plot.data.holediam.petleng$discrete <- "Hole Diameter"
plot.data.holediam.petleng$continuous <- "Petiole Length"
plot.data.holediam.petleng$tip_state <- factor(c("Several Large\nat Base", "One Large\nat Base", "All Large")[as.integer(plot.data.holediam.petleng$tip_state)], levels = c("Several Large\nat Base", "One Large\nat Base", "All Large"))

## ggplot(plot.data.holediam.petleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")

### Stem Area

model.list.holediam.stemarea <- list(er.bm1 = holediam.stemarea.er[[1]],
     er.bms = holediam.stemarea.er[[2]],
     er.ou1 = holediam.stemarea.er[[3]],
     er.oum = holediam.stemarea.er[[4]],
     er.ouma = holediam.stemarea.er[[5]],
     er.oumv = holediam.stemarea.er[[6]],
     sym.bm1 = holediam.stemarea.sym[[1]],
     sym.bms = holediam.stemarea.sym[[2]],
     sym.ou1 = holediam.stemarea.sym[[3]],
     sym.oum = holediam.stemarea.sym[[4]],
     sym.ouma = holediam.stemarea.sym[[5]],
     sym.oumv = holediam.stemarea.sym[[6]],
     ard.bm1 = holediam.stemarea.ard[[1]],
     ard.bms = holediam.stemarea.ard[[2]],
     ard.ou1 = holediam.stemarea.ard[[3]],
     ard.oum = holediam.stemarea.ard[[4]],
     ard.ouma = holediam.stemarea.ard[[5]],
     ard.oumv = holediam.stemarea.ard[[6]])

pars.holediam.stemarea <- getModelAvgParams(model.list.holediam.stemarea, force = FALSE)

plot.data.holediam.stemarea <- melt(pars.holediam.stemarea)
plot.data.holediam.stemarea$discrete <- "Hole Diameter"
plot.data.holediam.stemarea$continuous <- "Stem Area"
plot.data.holediam.stemarea$tip_state <- factor(c("Several Large\nat Base", "One Large\nat Base", "All Large")[as.integer(plot.data.holediam.stemarea$tip_state)], levels = c("Several Large\nat Base", "One Large\nat Base", "All Large"))

## ggplot(plot.data.holediam.stemarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")





## Reward

### Corola Length

model.list.reward.corleng <- list(er.bm1 = reward.corleng.er[[1]],
     er.bms = reward.corleng.er[[2]],
     er.ou1 = reward.corleng.er[[3]],
     er.oum = reward.corleng.er[[4]],
     er.ouma = reward.corleng.er[[5]],
     er.oumv = reward.corleng.er[[6]],
     sym.bm1 = reward.corleng.sym[[1]],
     sym.bms = reward.corleng.sym[[2]],
     sym.ou1 = reward.corleng.sym[[3]],
     sym.oum = reward.corleng.sym[[4]],
     sym.ouma = reward.corleng.sym[[5]],
     #sym.oumv = reward.corleng.sym[[6]],
     ard.bm1 = reward.corleng.ard[[1]],
     ard.bms = reward.corleng.ard[[2]],
     ard.ou1 = reward.corleng.ard[[3]],
     ard.oum = reward.corleng.ard[[4]],
     ard.ouma = reward.corleng.ard[[5]],
     ard.oumv = reward.corleng.ard[[6]])

pars.reward.corleng <- getModelAvgParams(model.list.reward.corleng, force = FALSE)

plot.data.reward.corleng <- melt(pars.reward.corleng)
plot.data.reward.corleng$discrete <- "Reward"
plot.data.reward.corleng$continuous <- "Corolla Length"
plot.data.reward.corleng$tip_state <- factor(c("Absent", "Present")[as.integer(plot.data.reward.corleng$tip_state) + 1], levels = c("Absent", "Present"))

## ggplot(plot.data.reward.corleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Leaf Area

model.list.reward.leafarea <- list(er.bm1 = reward.leafarea.er[[1]],
     er.bms = reward.leafarea.er[[2]],
     er.ou1 = reward.leafarea.er[[3]],
     er.oum = reward.leafarea.er[[4]],
     er.ouma = reward.leafarea.er[[5]],
     er.oumv = reward.leafarea.er[[6]],
     sym.bm1 = reward.leafarea.sym[[1]],
     sym.bms = reward.leafarea.sym[[2]],
     sym.ou1 = reward.leafarea.sym[[3]],
     sym.oum = reward.leafarea.sym[[4]],
     sym.ouma = reward.leafarea.sym[[5]],
     sym.oumv = reward.leafarea.sym[[6]],
     ard.bm1 = reward.leafarea.ard[[1]],
     ard.bms = reward.leafarea.ard[[2]],
     ard.ou1 = reward.leafarea.ard[[3]],
     ard.oum = reward.leafarea.ard[[4]],
     ard.ouma = reward.leafarea.ard[[5]],
     ard.oumv = reward.leafarea.ard[[6]])

pars.reward.leafarea <- getModelAvgParams(model.list.reward.leafarea, force = FALSE)

plot.data.reward.leafarea <- melt(pars.reward.leafarea)
plot.data.reward.leafarea$discrete <- "Reward"
plot.data.reward.leafarea$continuous <- "Leaf Area"
plot.data.reward.leafarea$tip_state <- factor(c("Absent", "Present")[as.integer(plot.data.reward.leafarea$tip_state) + 1], levels = c("Absent", "Present"))

## ggplot(plot.data.reward.leafarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Petiole Length

model.list.reward.petleng <- list(er.bm1 = reward.petleng.er[[1]],
     er.bms = reward.petleng.er[[2]],
     er.ou1 = reward.petleng.er[[3]],
     er.oum = reward.petleng.er[[4]],
     er.ouma = reward.petleng.er[[5]],
     er.oumv = reward.petleng.er[[6]],
     sym.bm1 = reward.petleng.sym[[1]],
     sym.bms = reward.petleng.sym[[2]],
     sym.ou1 = reward.petleng.sym[[3]],
     sym.oum = reward.petleng.sym[[4]],
     sym.ouma = reward.petleng.sym[[5]],
     sym.oumv = reward.petleng.sym[[6]],
     ard.bm1 = reward.petleng.ard[[1]],
     ard.bms = reward.petleng.ard[[2]],
     ard.ou1 = reward.petleng.ard[[3]],
     ard.oum = reward.petleng.ard[[4]],
     ard.ouma = reward.petleng.ard[[5]],
     ard.oumv = reward.petleng.ard[[6]])

pars.reward.petleng <- getModelAvgParams(model.list.reward.petleng, force = FALSE)

plot.data.reward.petleng <- melt(pars.reward.petleng)
plot.data.reward.petleng$discrete <- "Reward"
plot.data.reward.petleng$continuous <- "Petiole Length"
plot.data.reward.petleng$tip_state <- factor(c("Absent", "Present")[as.integer(plot.data.reward.petleng$tip_state) + 1], levels = c("Absent", "Present"))

## ggplot(plot.data.reward.petleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")

### Stem Area

model.list.reward.stemarea <- list(er.bm1 = reward.stemarea.er[[1]],
     er.bms = reward.stemarea.er[[2]],
     er.ou1 = reward.stemarea.er[[3]],
     er.oum = reward.stemarea.er[[4]],
     er.ouma = reward.stemarea.er[[5]],
     er.oumv = reward.stemarea.er[[6]],
     sym.bm1 = reward.stemarea.sym[[1]],
     sym.bms = reward.stemarea.sym[[2]],
     sym.ou1 = reward.stemarea.sym[[3]],
     sym.oum = reward.stemarea.sym[[4]],
     sym.ouma = reward.stemarea.sym[[5]],
     sym.oumv = reward.stemarea.sym[[6]],
     ard.bm1 = reward.stemarea.ard[[1]],
     ard.bms = reward.stemarea.ard[[2]],
     ard.ou1 = reward.stemarea.ard[[3]],
     ard.oum = reward.stemarea.ard[[4]],
     ard.ouma = reward.stemarea.ard[[5]],
     ard.oumv = reward.stemarea.ard[[6]])

pars.reward.stemarea <- getModelAvgParams(model.list.reward.stemarea, force = FALSE)

plot.data.reward.stemarea <- melt(pars.reward.stemarea)
plot.data.reward.stemarea$discrete <- "Reward"
plot.data.reward.stemarea$continuous <- "Stem Area"
plot.data.reward.stemarea$tip_state <- factor(c("Absent", "Present")[as.integer(plot.data.reward.stemarea$tip_state) + 1], levels = c("Absent", "Present"))

## ggplot(plot.data.reward.stemarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")




## Warts

### Corola Length

model.list.warts.corleng <- list(er.bm1 = warts.corleng.er[[1]],
     er.bms = warts.corleng.er[[2]],
     er.ou1 = warts.corleng.er[[3]],
     er.oum = warts.corleng.er[[4]],
     er.ouma = warts.corleng.er[[5]],
     er.oumv = warts.corleng.er[[6]],
     sym.bm1 = warts.corleng.sym[[1]],
     sym.bms = warts.corleng.sym[[2]],
     sym.ou1 = warts.corleng.sym[[3]],
     sym.oum = warts.corleng.sym[[4]],
     sym.ouma = warts.corleng.sym[[5]],
     sym.oumv = warts.corleng.sym[[6]],
     ard.bm1 = warts.corleng.ard[[1]],
     ard.bms = warts.corleng.ard[[2]],
     ard.ou1 = warts.corleng.ard[[3]],
     ard.oum = warts.corleng.ard[[4]],
     ard.ouma = warts.corleng.ard[[5]],
     ard.oumv = warts.corleng.ard[[6]])

pars.warts.corleng <- getModelAvgParams(model.list.warts.corleng, force = FALSE)

plot.data.warts.corleng <- melt(pars.warts.corleng)
plot.data.warts.corleng$discrete <- "Warts"
plot.data.warts.corleng$continuous <- "Corolla Length"
plot.data.warts.corleng$tip_state <- factor(c("Lost", "Differentiated", "Variable")[as.integer(plot.data.warts.corleng$tip_state)], levels = c("Lost", "Differentiated", "Variable"))

## ggplot(plot.data.warts.corleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Leaf Area

model.list.warts.leafarea <- list(er.bm1 = warts.leafarea.er[[1]],
     er.bms = warts.leafarea.er[[2]],
     er.ou1 = warts.leafarea.er[[3]],
     er.oum = warts.leafarea.er[[4]],
     er.ouma = warts.leafarea.er[[5]],
     er.oumv = warts.leafarea.er[[6]],
     sym.bm1 = warts.leafarea.sym[[1]],
     sym.bms = warts.leafarea.sym[[2]],
     sym.ou1 = warts.leafarea.sym[[3]],
     sym.oum = warts.leafarea.sym[[4]],
     sym.ouma = warts.leafarea.sym[[5]],
     sym.oumv = warts.leafarea.sym[[6]],
     ard.bm1 = warts.leafarea.ard[[1]],
     ard.bms = warts.leafarea.ard[[2]],
     ard.ou1 = warts.leafarea.ard[[3]],
     ard.oum = warts.leafarea.ard[[4]],
     ard.ouma = warts.leafarea.ard[[5]],
     ard.oumv = warts.leafarea.ard[[6]])

pars.warts.leafarea <- getModelAvgParams(model.list.warts.leafarea, force = FALSE)

plot.data.warts.leafarea <- melt(pars.warts.leafarea)
plot.data.warts.leafarea$discrete <- "Warts"
plot.data.warts.leafarea$continuous <- "Leaf Area"
plot.data.warts.leafarea$tip_state <- factor(c("Lost", "Differentiated", "Variable")[as.integer(plot.data.warts.leafarea$tip_state)], levels = c("Lost", "Differentiated", "Variable"))

## ggplot(plot.data.warts.leafarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Petiole Length

model.list.warts.petleng <- list(er.bm1 = warts.petleng.er[[1]],
     er.bms = warts.petleng.er[[2]],
     er.ou1 = warts.petleng.er[[3]],
     er.oum = warts.petleng.er[[4]],
     er.ouma = warts.petleng.er[[5]],
     er.oumv = warts.petleng.er[[6]],
     sym.bm1 = warts.petleng.sym[[1]],
     sym.bms = warts.petleng.sym[[2]],
     sym.ou1 = warts.petleng.sym[[3]],
     sym.oum = warts.petleng.sym[[4]],
     sym.ouma = warts.petleng.sym[[5]],
     sym.oumv = warts.petleng.sym[[6]],
     ard.bm1 = warts.petleng.ard[[1]],
     ard.bms = warts.petleng.ard[[2]],
     ard.ou1 = warts.petleng.ard[[3]],
     ard.oum = warts.petleng.ard[[4]],
     ard.ouma = warts.petleng.ard[[5]],
     ard.oumv = warts.petleng.ard[[6]])

pars.warts.petleng <- getModelAvgParams(model.list.warts.petleng, force = FALSE)

plot.data.warts.petleng <- melt(pars.warts.petleng)
plot.data.warts.petleng$discrete <- "Warts"
plot.data.warts.petleng$continuous <- "Petiole Length"
plot.data.warts.petleng$tip_state <- factor(c("Lost", "Differentiated", "Variable")[as.integer(plot.data.warts.petleng$tip_state)], levels = c("Lost", "Differentiated", "Variable"))

## ggplot(plot.data.warts.petleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")

### Stem Area

model.list.warts.stemarea <- list(er.bm1 = warts.stemarea.er[[1]],
     er.bms = warts.stemarea.er[[2]],
     er.ou1 = warts.stemarea.er[[3]],
     er.oum = warts.stemarea.er[[4]],
     er.ouma = warts.stemarea.er[[5]],
     er.oumv = warts.stemarea.er[[6]],
     sym.bm1 = warts.stemarea.sym[[1]],
     sym.bms = warts.stemarea.sym[[2]],
     sym.ou1 = warts.stemarea.sym[[3]],
     sym.oum = warts.stemarea.sym[[4]],
     sym.ouma = warts.stemarea.sym[[5]],
     #sym.oumv = warts.stemarea.sym[[6]],
     ard.bm1 = warts.stemarea.ard[[1]],
     ard.bms = warts.stemarea.ard[[2]],
     ard.ou1 = warts.stemarea.ard[[3]],
     ard.oum = warts.stemarea.ard[[4]],
     ard.ouma = warts.stemarea.ard[[5]],
     ard.oumv = warts.stemarea.ard[[6]])

pars.warts.stemarea <- getModelAvgParams(model.list.warts.stemarea, force = FALSE)

plot.data.warts.stemarea <- melt(pars.warts.stemarea)
plot.data.warts.stemarea$discrete <- "Warts"
plot.data.warts.stemarea$continuous <- "Stem Area"
plot.data.warts.stemarea$tip_state <- factor(c("Lost", "Differentiated", "Variable")[as.integer(plot.data.warts.stemarea$tip_state)], levels = c("Lost", "Differentiated", "Variable"))

## ggplot(plot.data.warts.stemarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")




theta.plot <- rbind(plot.data.domgrow.corleng[grep("theta", plot.data.domgrow.corleng$variable),],
                    plot.data.domgrow.leafarea[grep("theta", plot.data.domgrow.leafarea$variable),],
                    plot.data.domgrow.petleng[grep("theta", plot.data.domgrow.petleng$variable),],
                    plot.data.domgrow.stemarea[grep("theta", plot.data.domgrow.stemarea$variable),],
                    plot.data.holediam.corleng[grep("theta", plot.data.holediam.corleng$variable),],
                    plot.data.holediam.leafarea[grep("theta", plot.data.holediam.leafarea$variable),],
                    plot.data.holediam.petleng[grep("theta", plot.data.holediam.petleng$variable),],
                    plot.data.holediam.stemarea[grep("theta", plot.data.holediam.stemarea$variable),],
                    plot.data.reward.corleng[grep("theta", plot.data.reward.corleng$variable),],
                    plot.data.reward.leafarea[grep("theta", plot.data.reward.leafarea$variable),],
                    plot.data.reward.petleng[grep("theta", plot.data.reward.petleng$variable),],
                    plot.data.reward.stemarea[grep("theta", plot.data.reward.stemarea$variable),],
                    plot.data.warts.corleng[grep("theta", plot.data.warts.corleng$variable),],
                    plot.data.warts.leafarea[grep("theta", plot.data.warts.leafarea$variable),],
                    plot.data.warts.petleng[grep("theta", plot.data.warts.petleng$variable),],
                    plot.data.warts.stemarea[grep("theta", plot.data.warts.stemarea$variable),]
                    )


alpha.plot <- rbind(plot.data.domgrow.corleng[grep("alpha", plot.data.domgrow.corleng$variable),],
                    plot.data.domgrow.leafarea[grep("alpha", plot.data.domgrow.leafarea$variable),],
                    plot.data.domgrow.petleng[grep("alpha", plot.data.domgrow.petleng$variable),],
                    plot.data.domgrow.stemarea[grep("alpha", plot.data.domgrow.stemarea$variable),],
                    plot.data.holediam.corleng[grep("alpha", plot.data.holediam.corleng$variable),],
                    plot.data.holediam.leafarea[grep("alpha", plot.data.holediam.leafarea$variable),],
                    plot.data.holediam.petleng[grep("alpha", plot.data.holediam.petleng$variable),],
                    plot.data.holediam.stemarea[grep("alpha", plot.data.holediam.stemarea$variable),],
                    plot.data.reward.corleng[grep("alpha", plot.data.reward.corleng$variable),],
                    plot.data.reward.leafarea[grep("alpha", plot.data.reward.leafarea$variable),],
                    plot.data.reward.petleng[grep("alpha", plot.data.reward.petleng$variable),],
                    plot.data.reward.stemarea[grep("alpha", plot.data.reward.stemarea$variable),],
                    plot.data.warts.corleng[grep("alpha", plot.data.warts.corleng$variable),],
                    plot.data.warts.leafarea[grep("alpha", plot.data.warts.leafarea$variable),],
                    plot.data.warts.petleng[grep("alpha", plot.data.warts.petleng$variable),],
                    plot.data.warts.stemarea[grep("alpha", plot.data.warts.stemarea$variable),]
                    )


sigma.plot <- rbind(plot.data.domgrow.corleng[grep("sigma", plot.data.domgrow.corleng$variable),],
                    plot.data.domgrow.leafarea[grep("sigma", plot.data.domgrow.leafarea$variable),],
                    plot.data.domgrow.petleng[grep("sigma", plot.data.domgrow.petleng$variable),],
                    plot.data.domgrow.stemarea[grep("sigma", plot.data.domgrow.stemarea$variable),],
                    plot.data.holediam.corleng[grep("sigma", plot.data.holediam.corleng$variable),],
                    plot.data.holediam.leafarea[grep("sigma", plot.data.holediam.leafarea$variable),],
                    plot.data.holediam.petleng[grep("sigma", plot.data.holediam.petleng$variable),],
                    plot.data.holediam.stemarea[grep("sigma", plot.data.holediam.stemarea$variable),],
                    plot.data.reward.corleng[grep("sigma", plot.data.reward.corleng$variable),],
                    plot.data.reward.leafarea[grep("sigma", plot.data.reward.leafarea$variable),],
                    plot.data.reward.petleng[grep("sigma", plot.data.reward.petleng$variable),],
                    plot.data.reward.stemarea[grep("sigma", plot.data.reward.stemarea$variable),],
                    plot.data.warts.corleng[grep("sigma", plot.data.warts.corleng$variable),],
                    plot.data.warts.leafarea[grep("sigma", plot.data.warts.leafarea$variable),],
                    plot.data.warts.petleng[grep("sigma", plot.data.warts.petleng$variable),],
                    plot.data.warts.stemarea[grep("sigma", plot.data.warts.stemarea$variable),]
                    )


ggplot(theta.plot, aes(x = tip_state, y = value, colour = tip_state)) +
    geom_point(size = 5, shape = 21) +
    stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    scale_colour_manual(values = c(met.brewer("Greek", n = 2, type = "continuous"),
                                   met.brewer("Veronese", n = 3, type = "continuous"),
                                   met.brewer("Hiroshige", n = 2, type = "continuous"),
                                   met.brewer("Isfahan2", n = 3, type = "continuous"))) +
    scale_alpha_manual(values = 0.5) +
    theme_classic() +
    facet_wrap(discrete ~ continuous, scales = "free")

## FIX LEGEND ORDER (COLOURS, LABELS, ETC.)
