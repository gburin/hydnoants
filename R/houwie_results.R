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
library("ggh4x")

i_am("R/houwie_results.R")

for(i in 1:length(list.files(here("output/houwie/main_analysis"), pattern = ".RDS"))){
    assign(gsub("-", ".", gsub("-nohid.RDS", "", list.files(here("output/houwie/main_analysis"), pattern = ".RDS")[i])),
           readRDS(paste0("output/houwie/main_analysis/", list.files(here("output/houwie/main_analysis"), pattern = ".RDS")[i]))
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

bic.index.dgcl <- sapply(model.list.domgrow.corleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
domgrow.corleng.table <- getModelTable(model.list.domgrow.corleng[-which(abs(bic.index.dgcl) > 1e04)])

if(any(abs(bic.index.dgcl) > 1e04, na.rm = TRUE)){
    pars.domgrow.corleng <- getModelAvgParams(model.list.domgrow.corleng[-which(abs(bic.index.dgcl) > 1e04)], force = FALSE)
} else {                                                                              pars.domgrow.corleng <- getModelAvgParams(model.list.domgrow.corleng, force = FALSE)}

plot.data.domgrow.corleng <- reshape2::melt(pars.domgrow.corleng)
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

bic.index.dgla <- sapply(model.list.domgrow.leafarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
domgrow.leafarea.table <- getModelTable(model.list.domgrow.leafarea[-which(abs(bic.index.dgla) > 1e04)])

if(any(abs(bic.index.dgla) > 1e04, na.rm = TRUE)){
    pars.domgrow.leafarea <- getModelAvgParams(model.list.domgrow.leafarea[-which(abs(bic.index.dgla) > 1e04)], force = FALSE)
} else {                                                                              pars.domgrow.leafarea <- getModelAvgParams(model.list.domgrow.leafarea, force = FALSE)}

plot.data.domgrow.leafarea <- reshape2::melt(pars.domgrow.leafarea)
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

bic.index.dgpl <- sapply(model.list.domgrow.petleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
domgrow.petleng.table <- getModelTable(model.list.domgrow.petleng[-which(abs(bic.index.dgpl) > 1e04)])

if(any(abs(bic.index.dgpl) > 1e04, na.rm = TRUE)){
    pars.domgrow.petleng <- getModelAvgParams(model.list.domgrow.petleng[-which(abs(bic.index.dgpl) > 1e04)], force = FALSE)
} else {                                                                              pars.domgrow.petleng <- getModelAvgParams(model.list.domgrow.petleng, force = FALSE)}

plot.data.domgrow.petleng <- reshape2::melt(pars.domgrow.petleng)
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

bic.index.dgsa <- sapply(model.list.domgrow.stemarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
domgrow.stemarea.table <- getModelTable(model.list.domgrow.stemarea[-which(abs(bic.index.dgsa) > 1e04)])

if(any(abs(bic.index.dgsa) > 1e04, na.rm = TRUE)){
    pars.domgrow.stemarea <- getModelAvgParams(model.list.domgrow.stemarea[-which(abs(bic.index.dgsa) > 1e04)], force = FALSE)
} else {                                                                              pars.domgrow.stemarea <- getModelAvgParams(model.list.domgrow.stemarea, force = FALSE)}

plot.data.domgrow.stemarea <- reshape2::melt(pars.domgrow.stemarea)
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

bic.index.hdcl <- sapply(model.list.holediam.corleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
holediam.corleng.table <- getModelTable(model.list.holediam.corleng[-which(abs(bic.index.hdcl) > 1e04)])

if(any(abs(bic.index.hdcl) > 1e04, na.rm = TRUE)){
    pars.holediam.corleng <- getModelAvgParams(model.list.holediam.corleng[-which(abs(bic.index.hdcl) > 1e04)], force = FALSE)
} else {                                                                              pars.holediam.corleng <- getModelAvgParams(model.list.holediam.corleng, force = FALSE)}

plot.data.holediam.corleng <- reshape2::melt(pars.holediam.corleng)
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

bic.index.hdla <- sapply(model.list.holediam.leafarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
holediam.leafarea.table <- getModelTable(model.list.holediam.leafarea[-which(abs(bic.index.hdla) > 1e04)])

if(any(abs(bic.index.hdla) > 1e04, na.rm = TRUE)){
    pars.holediam.leafarea <- getModelAvgParams(model.list.holediam.leafarea[-which(abs(bic.index.hdla) > 1e04)], force = FALSE)
} else {                                                                              pars.holediam.leafarea <- getModelAvgParams(model.list.holediam.leafarea, force = FALSE)}

plot.data.holediam.leafarea <- reshape2::melt(pars.holediam.leafarea)
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

bic.index.hdpl <- sapply(model.list.holediam.petleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
holediam.petleng.table <- getModelTable(model.list.holediam.petleng[-which(abs(bic.index.hdpl) > 1e04)])

if(any(abs(bic.index.hdpl) > 1e04, na.rm = TRUE)){
    pars.holediam.petleng <- getModelAvgParams(model.list.holediam.petleng[-which(abs(bic.index.hdpl) > 1e04)], force = FALSE)
} else {                                                                              pars.holediam.petleng <- getModelAvgParams(model.list.holediam.petleng, force = FALSE)}

plot.data.holediam.petleng <- reshape2::melt(pars.holediam.petleng)
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

bic.index.hdsa <- sapply(model.list.holediam.stemarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
holediam.stemarea.table <- getModelTable(model.list.holediam.stemarea[-which(abs(bic.index.hdsa) > 1e04)])

if(any(abs(bic.index.hdsa) > 1e04, na.rm = TRUE)){
    pars.holediam.stemarea <- getModelAvgParams(model.list.holediam.stemarea[-which(abs(bic.index.hdsa) > 1e04)], force = FALSE)
} else {                                                                              pars.holediam.stemarea <- getModelAvgParams(model.list.holediam.stemarea, force = FALSE)}

plot.data.holediam.stemarea <- reshape2::melt(pars.holediam.stemarea)
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

bic.index.recl <- sapply(model.list.reward.corleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
reward.corleng.table <- getModelTable(model.list.reward.corleng[-which(abs(bic.index.recl) > 1e04)])

if(any(abs(bic.index.recl) > 1e04, na.rm = TRUE)){
    pars.reward.corleng <- getModelAvgParams(model.list.reward.corleng[-which(abs(bic.index.recl) > 1e04)], force = FALSE)
} else {                                                                              pars.reward.corleng <- getModelAvgParams(model.list.reward.corleng, force = FALSE)}

plot.data.reward.corleng <- reshape2::melt(pars.reward.corleng)
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

bic.index.rela <- sapply(model.list.reward.leafarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
reward.leafarea.table <- getModelTable(model.list.reward.leafarea[-which(abs(bic.index.rela) > 1e04)])

if(any(abs(bic.index.rela) > 1e04, na.rm = TRUE)){
    pars.reward.leafarea <- getModelAvgParams(model.list.reward.leafarea[-which(abs(bic.index.rela) > 1e04)], force = FALSE)
} else {                                                                              pars.reward.leafarea <- getModelAvgParams(model.list.reward.leafarea, force = FALSE)}

plot.data.reward.leafarea <- reshape2::melt(pars.reward.leafarea)
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

bic.index.repl <- sapply(model.list.reward.petleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
reward.petleng.table <- getModelTable(model.list.reward.petleng[-which(abs(bic.index.repl) > 1e04)])

if(any(abs(bic.index.repl) > 1e04, na.rm = TRUE)){
    pars.reward.petleng <- getModelAvgParams(model.list.reward.petleng[-which(abs(bic.index.repl) > 1e04)], force = FALSE)
} else {                                                                              pars.reward.petleng <- getModelAvgParams(model.list.reward.petleng, force = FALSE)}

plot.data.reward.petleng <- reshape2::melt(pars.reward.petleng)
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

bic.index.resa <- sapply(model.list.reward.stemarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
reward.stemarea.table <- getModelTable(model.list.reward.stemarea[-which(abs(bic.index.resa) > 1e04)])

if(any(abs(bic.index.resa) > 1e04, na.rm = TRUE)){
    pars.reward.stemarea <- getModelAvgParams(model.list.reward.stemarea[-which(abs(bic.index.resa) > 1e04)], force = FALSE)
} else {                                                                              pars.reward.stemarea <- getModelAvgParams(model.list.reward.stemarea, force = FALSE)}

plot.data.reward.stemarea <- reshape2::melt(pars.reward.stemarea)
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

bic.index.wacl <- sapply(model.list.warts.corleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
warts.corleng.table <- getModelTable(model.list.warts.corleng[-which(abs(bic.index.wacl) > 1e04)])

if(any(abs(bic.index.wacl) > 1e04, na.rm = TRUE)){
    pars.warts.corleng <- getModelAvgParams(model.list.warts.corleng[-which(abs(bic.index.wacl) > 1e04)], force = FALSE)
} else {                                                                              pars.warts.corleng <- getModelAvgParams(model.list.warts.corleng, force = FALSE)}

plot.data.warts.corleng <- reshape2::melt(pars.warts.corleng)
plot.data.warts.corleng$discrete <- "Warts"
plot.data.warts.corleng$continuous <- "Corolla Length"
plot.data.warts.corleng$tip_state <- factor(c("Variable", "Differentiated", "Lost")[as.integer(plot.data.warts.corleng$tip_state)], levels = c("Variable", "Differentiated", "Lost"))

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

bic.index.wala <- sapply(model.list.warts.leafarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
warts.leafarea.table <- getModelTable(model.list.warts.leafarea) ## no model has a BIC of over 1e04 in absolute

if(any(abs(bic.index.wala) > 1e04, na.rm = TRUE)){
    pars.warts.leafarea <- getModelAvgParams(model.list.warts.leafarea[-which(abs(bic.index.wala) > 1e04)], force = FALSE)
} else {                                                                              pars.warts.leafarea <- getModelAvgParams(model.list.warts.leafarea, force = FALSE)}

plot.data.warts.leafarea <- reshape2::melt(pars.warts.leafarea)
plot.data.warts.leafarea$discrete <- "Warts"
plot.data.warts.leafarea$continuous <- "Leaf Area"
plot.data.warts.leafarea$tip_state <- factor(c("Variable", "Differentiated", "Lost")[as.integer(plot.data.warts.leafarea$tip_state)], levels = c("Variable", "Differentiated", "Lost"))

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

bic.index.wapl <- sapply(model.list.warts.petleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
warts.petleng.table <- getModelTable(model.list.warts.petleng[-which(abs(bic.index.wapl) > 1e04)])

if(any(abs(bic.index.wapl) > 1e04, na.rm = TRUE)){
    pars.warts.petleng <- getModelAvgParams(model.list.warts.petleng[-which(abs(bic.index.wapl) > 1e04)], force = FALSE)
} else {                                                                              pars.warts.petleng <- getModelAvgParams(model.list.warts.petleng, force = FALSE)}

plot.data.warts.petleng <- reshape2::melt(pars.warts.petleng)
plot.data.warts.petleng$discrete <- "Warts"
plot.data.warts.petleng$continuous <- "Petiole Length"
plot.data.warts.petleng$tip_state <- factor(c("Variable", "Differentiated", "Lost")[as.integer(plot.data.warts.petleng$tip_state)], levels = c("Variable", "Differentiated", "Lost"))

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
     sym.oumv = warts.stemarea.sym[[6]],
     ard.bm1 = warts.stemarea.ard[[1]],
     ard.bms = warts.stemarea.ard[[2]],
     ard.ou1 = warts.stemarea.ard[[3]],
     ard.oum = warts.stemarea.ard[[4]],
     ard.ouma = warts.stemarea.ard[[5]],
     ard.oumv = warts.stemarea.ard[[6]])

bic.index.wasa <- sapply(model.list.warts.stemarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
warts.stemarea.table <- getModelTable(model.list.warts.stemarea[-which(abs(bic.index.wasa) > 1e04)])

if(any(abs(bic.index.wasa) > 1e04, na.rm = TRUE)){
    pars.warts.stemarea <- getModelAvgParams(model.list.warts.stemarea[-which(abs(bic.index.wasa) > 1e04)], force = FALSE)
} else {                                                                              pars.warts.stemarea <- getModelAvgParams(model.list.warts.stemarea, force = FALSE)}

plot.data.warts.stemarea <- reshape2::melt(pars.warts.stemarea)
plot.data.warts.stemarea$discrete <- "Warts"
plot.data.warts.stemarea$continuous <- "Stem Area"
plot.data.warts.stemarea$tip_state <- factor(c("Variable", "Differentiated", "Lost")[as.integer(plot.data.warts.stemarea$tip_state)], levels = c("Variable", "Differentiated", "Lost"))

## ggplot(plot.data.warts.stemarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


load(here("output/houwie/main_analysis/param_values_main.RData"))

theta.plot <- rbind(plot.data.domgrow.corleng[grep("expected_mean", plot.data.domgrow.corleng$variable),],
                    plot.data.domgrow.leafarea[grep("expected_mean", plot.data.domgrow.leafarea$variable),],
                    plot.data.domgrow.petleng[grep("expected_mean", plot.data.domgrow.petleng$variable),],
                    plot.data.domgrow.stemarea[grep("expected_mean", plot.data.domgrow.stemarea$variable),],
                    plot.data.holediam.corleng[grep("expected_mean", plot.data.holediam.corleng$variable),],
                    plot.data.holediam.leafarea[grep("expected_mean", plot.data.holediam.leafarea$variable),],
                    plot.data.holediam.petleng[grep("expected_mean", plot.data.holediam.petleng$variable),],
                    plot.data.holediam.stemarea[grep("expected_mean", plot.data.holediam.stemarea$variable),],
                    plot.data.reward.corleng[grep("expected_mean", plot.data.reward.corleng$variable),],
                    plot.data.reward.leafarea[grep("expected_mean", plot.data.reward.leafarea$variable),],
                    plot.data.reward.petleng[grep("expected_mean", plot.data.reward.petleng$variable),],
                    plot.data.reward.stemarea[grep("expected_mean", plot.data.reward.stemarea$variable),],
                    plot.data.warts.corleng[grep("expected_mean", plot.data.warts.corleng$variable),],
                    plot.data.warts.leafarea[grep("expected_mean", plot.data.warts.leafarea$variable),],
                    plot.data.warts.petleng[grep("expected_mean", plot.data.warts.petleng$variable),],
                    plot.data.warts.stemarea[grep("expected_mean", plot.data.warts.stemarea$variable),]
                    )
theta.plot$discrete <- factor(theta.plot$discrete, levels = c("Domatium Growth", "Reward", "Warts", "Hole Diameter"))

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
alpha.plot$discrete <- factor(alpha.plot$discrete, levels = c("Domatium Growth", "Reward", "Warts", "Hole Diameter"))

sigma.plot <- rbind(plot.data.domgrow.corleng[grep("expected_var", plot.data.domgrow.corleng$variable),],
                    plot.data.domgrow.leafarea[grep("expected_var", plot.data.domgrow.leafarea$variable),],
                    plot.data.domgrow.petleng[grep("expected_var", plot.data.domgrow.petleng$variable),],
                    plot.data.domgrow.stemarea[grep("expected_var", plot.data.domgrow.stemarea$variable),],
                    plot.data.holediam.corleng[grep("expected_var", plot.data.holediam.corleng$variable),],
                    plot.data.holediam.leafarea[grep("expected_var", plot.data.holediam.leafarea$variable),],
                    plot.data.holediam.petleng[grep("expected_var", plot.data.holediam.petleng$variable),],
                    plot.data.holediam.stemarea[grep("expected_var", plot.data.holediam.stemarea$variable),],
                    plot.data.reward.corleng[grep("expected_var", plot.data.reward.corleng$variable),],
                    plot.data.reward.leafarea[grep("expected_var", plot.data.reward.leafarea$variable),],
                    plot.data.reward.petleng[grep("expected_var", plot.data.reward.petleng$variable),],
                    plot.data.reward.stemarea[grep("expected_var", plot.data.reward.stemarea$variable),],
                    plot.data.warts.corleng[grep("expected_var", plot.data.warts.corleng$variable),],
                    plot.data.warts.leafarea[grep("expected_var", plot.data.warts.leafarea$variable),],
                    plot.data.warts.petleng[grep("expected_var", plot.data.warts.petleng$variable),],
                    plot.data.warts.stemarea[grep("expected_var", plot.data.warts.stemarea$variable),]
                    )
sigma.plot$discrete <- factor(sigma.plot$discrete, levels = c("Domatium Growth", "Reward", "Warts", "Hole Diameter"))


ggplot(theta.plot, aes(x = tip_state, y = value, colour = tip_state, alpha = tip_state)) +
    geom_jitter(size = 2, shape = 21, width = 0.1, alpha = 0.3) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 5) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1, colour = tip_state), width = 0.3) +
    scale_colour_manual(values = c(met.brewer("Greek", n = 2, type = "continuous"),
                                   met.brewer("Veronese", n = 3, type = "continuous"),
                                   met.brewer("Hiroshige", n = 2, type = "continuous"),
                                   met.brewer("Isfahan2", n = 3, type = "continuous"))) +
    scale_alpha_manual(values = rep(1, 10)) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    ggh4x::facet_grid2(c("discrete", "continuous"), scales = "free", independent = "all") +
    cowplot::theme_cowplot() +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "white", colour = "darkgrey", linewidth = 1.5),
          strip.text = element_text(size = 14))

ggsave(here("manuscript/figs/houwie/theta_houwie.png"), bg = "white", dpi = 300, height = 10, width = 18, units = "in")
ggsave(here("manuscript/figs/houwie/theta_houwie.pdf"), bg = "white", height = 10, width = 18, units = "in")

ggplot(alpha.plot, aes(x = tip_state, y = value, colour = tip_state, alpha = tip_state)) +
    geom_point(size = 5, shape = 21) +
    stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1, colour = tip_state), width = 0.15) +
    scale_colour_manual(values = c(met.brewer("Greek", n = 2, type = "continuous"),
                                   met.brewer("Veronese", n = 3, type = "continuous"),
                                   met.brewer("Hiroshige", n = 2, type = "continuous"),
                                   met.brewer("Isfahan2", n = 3, type = "continuous"))) +
    scale_alpha_manual(values = rep(1, 10)) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    ggh4x::facet_grid2(c("discrete", "continuous"), scales = "free", independent = "all") +
    cowplot::theme_cowplot() +
    theme(legend.position = "none")

ggsave(here("manuscript/figs/houwie/alpha_houwie.png"), bg = "white", dpi = 300, height = 10, width = 18, units = "in")
ggsave(here("manuscript/figs/houwie/alpha_houwie.pdf"), bg = "white", height = 10, width = 18, units = "in")


ggplot(sigma.plot, aes(x = tip_state, y = value, colour = tip_state, alpha = tip_state)) +
    geom_jitter(size = 2, shape = 21, width = 0.1, alpha = 0.3) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 5) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1, colour = tip_state), width = 0.3) +
    scale_colour_manual(values = c(met.brewer("Greek", n = 2, type = "continuous"),
                                   met.brewer("Veronese", n = 3, type = "continuous"),
                                   met.brewer("Hiroshige", n = 2, type = "continuous"),
                                   met.brewer("Isfahan2", n = 3, type = "continuous"))) +
    scale_alpha_manual(values = rep(1, 10)) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    ggh4x::facet_grid2(c("discrete", "continuous"), scales = "free", independent = "all") +
    cowplot::theme_cowplot() +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "white", colour = "darkgrey", linewidth = 1.5),
          strip.text = element_text(size = 14))

ggsave(here("manuscript/figs/houwie/sigma_houwie.png"), bg = "white", dpi = 300, height = 10, width = 18, units = "in")
ggsave(here("manuscript/figs/houwie/sigma_houwie.pdf"), bg = "white", height = 10, width = 18, units = "in")

## PGLS to test for differences in parameter values
load(here("output/houwie/main_analysis/pgls_data_main.RData"))
mcc.tree <- read.tree(here("data/MCCtree135taxa.tre"))

## Expected Mean

### Domatium Growth

dgcl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.domgrow.corleng))]),
                                   data.frame(pars.domgrow.corleng,
                                              species = rownames(pars.domgrow.corleng)),
                                   names.col = "species")
summary(pgls(expected_mean ~ tip_state, data = dgcl.comp.data, lambda = "ML"))

dgla.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.domgrow.leafarea))]),
                                   data.frame(pars.domgrow.leafarea,
                                              species = rownames(pars.domgrow.leafarea)),
                                   names.col = "species")
summary(pgls(expected_mean ~ tip_state, data = dgla.comp.data, lambda = "ML"))

dgpl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.domgrow.petleng))]),
                                   data.frame(pars.domgrow.petleng,
                                              species = rownames(pars.domgrow.petleng)),
                                   names.col = "species")
summary(pgls(expected_mean ~ tip_state, data = dgpl.comp.data, lambda = "ML"))

dgsa.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.domgrow.stemarea))]),
                                   data.frame(pars.domgrow.stemarea,
                                              species = rownames(pars.domgrow.stemarea)),
                                   names.col = "species")
summary(pgls(expected_mean ~ tip_state, data = dgsa.comp.data, lambda = "ML"))


### Hole Diameter

hdcl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.holediam.corleng))]),
                                   data.frame(pars.holediam.corleng,
                                              species = rownames(pars.holediam.corleng)),
                                   names.col = "species")
phylANOVA(hdcl.comp.data$phy, x = setNames(hdcl.comp.data$data$tip_state, rownames(hdcl.comp.data$data)), y = setNames(hdcl.comp.data$data$expected_mean, rownames(hdcl.comp.data$data)))

hdla.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.holediam.leafarea))]),
                                   data.frame(pars.holediam.leafarea,
                                              species = rownames(pars.holediam.leafarea)),
                                   names.col = "species")
phylANOVA(hdla.comp.data$phy, x = setNames(hdla.comp.data$data$tip_state, rownames(hdla.comp.data$data)), y = setNames(hdla.comp.data$data$expected_mean, rownames(hdla.comp.data$data)))

hdpl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.holediam.petleng))]),
                                   data.frame(pars.holediam.petleng,
                                              species = rownames(pars.holediam.petleng)),
                                   names.col = "species")
phylANOVA(hdpl.comp.data$phy, x = setNames(hdpl.comp.data$data$tip_state, rownames(hdpl.comp.data$data)), y = setNames(hdpl.comp.data$data$expected_mean, rownames(hdpl.comp.data$data)))

hdsa.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.holediam.stemarea))]),
                                   data.frame(pars.holediam.stemarea,
                                              species = rownames(pars.holediam.stemarea)),
                                   names.col = "species")
phylANOVA(hdsa.comp.data$phy, x = setNames(hdsa.comp.data$data$tip_state, rownames(hdsa.comp.data$data)), y = setNames(hdsa.comp.data$data$expected_mean, rownames(hdsa.comp.data$data)))


### Reward

rewcl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.reward.corleng))]),
                                   data.frame(pars.reward.corleng,
                                              species = rownames(pars.reward.corleng)),
                                   names.col = "species")
summary(pgls(expected_mean ~ tip_state, data = rewcl.comp.data, lambda = "ML"))

rewla.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.reward.leafarea))]),
                                   data.frame(pars.reward.leafarea,
                                              species = rownames(pars.reward.leafarea)),
                                   names.col = "species")
summary(pgls(expected_mean ~ tip_state, data = rewla.comp.data, lambda = "ML"))

rewpl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.reward.petleng))]),
                                   data.frame(pars.reward.petleng,
                                              species = rownames(pars.reward.petleng)),
                                   names.col = "species")
summary(pgls(expected_mean ~ tip_state, data = rewpl.comp.data, lambda = "ML"))

rewsa.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.reward.stemarea))]),
                                   data.frame(pars.reward.stemarea,
                                              species = rownames(pars.reward.stemarea)),
                                   names.col = "species")
summary(pgls(expected_mean ~ tip_state, data = rewsa.comp.data, lambda = "ML"))


### Warts

warcl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.warts.corleng))]),
                                   data.frame(pars.warts.corleng,
                                              species = rownames(pars.warts.corleng)),
                                   names.col = "species")
phylANOVA(warcl.comp.data$phy, x = setNames(warcl.comp.data$data$tip_state, rownames(warcl.comp.data$data)), y = setNames(warcl.comp.data$data$expected_mean, rownames(warcl.comp.data$data)))

warla.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.warts.leafarea))]),
                                   data.frame(pars.warts.leafarea,
                                              species = rownames(pars.warts.leafarea)),
                                   names.col = "species")
phylANOVA(warla.comp.data$phy, x = setNames(warla.comp.data$data$tip_state, rownames(warla.comp.data$data)), y = setNames(warla.comp.data$data$expected_mean, rownames(warla.comp.data$data)))

warpl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.warts.petleng))]),
                                   data.frame(pars.warts.petleng,
                                              species = rownames(pars.warts.petleng)),
                                   names.col = "species")
phylANOVA(warpl.comp.data$phy, x = setNames(warpl.comp.data$data$tip_state, rownames(warpl.comp.data$data)), y = setNames(warpl.comp.data$data$expected_mean, rownames(warpl.comp.data$data)))

warsa.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.warts.stemarea))]),
                                   data.frame(pars.warts.stemarea,
                                              species = rownames(pars.warts.stemarea)),
                                   names.col = "species")
phylANOVA(warsa.comp.data$phy, x = setNames(warsa.comp.data$data$tip_state, rownames(warsa.comp.data$data)), y = setNames(warsa.comp.data$data$expected_mean, rownames(warsa.comp.data$data)))


## Expected Var

### Domatium Growth

dgcl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.domgrow.corleng))]),
                                   data.frame(pars.domgrow.corleng,
                                              species = rownames(pars.domgrow.corleng)),
                                   names.col = "species")
summary(pgls(expected_var ~ tip_state, data = dgcl.comp.data, lambda = "ML"))

dgla.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.domgrow.leafarea))]),
                                   data.frame(pars.domgrow.leafarea,
                                              species = rownames(pars.domgrow.leafarea)),
                                   names.col = "species")
summary(pgls(expected_var ~ tip_state, data = dgla.comp.data, lambda = "ML"))

dgpl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.domgrow.petleng))]),
                                   data.frame(pars.domgrow.petleng,
                                              species = rownames(pars.domgrow.petleng)),
                                   names.col = "species")
summary(pgls(expected_var ~ tip_state, data = dgpl.comp.data, lambda = "ML"))

dgsa.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.domgrow.stemarea))]),
                                   data.frame(pars.domgrow.stemarea,
                                              species = rownames(pars.domgrow.stemarea)),
                                   names.col = "species")
summary(pgls(expected_var ~ tip_state, data = dgsa.comp.data, lambda = "ML"))


### Hole Diameter

hdcl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.holediam.corleng))]),
                                   data.frame(pars.holediam.corleng,
                                              species = rownames(pars.holediam.corleng)),
                                   names.col = "species")
phylANOVA(hdcl.comp.data$phy, x = setNames(hdcl.comp.data$data$tip_state, rownames(hdcl.comp.data$data)), y = setNames(hdcl.comp.data$data$expected_var, rownames(hdcl.comp.data$data)))

hdla.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.holediam.leafarea))]),
                                   data.frame(pars.holediam.leafarea,
                                              species = rownames(pars.holediam.leafarea)),
                                   names.col = "species")
phylANOVA(hdla.comp.data$phy, x = setNames(hdla.comp.data$data$tip_state, rownames(hdla.comp.data$data)), y = setNames(hdla.comp.data$data$expected_var, rownames(hdla.comp.data$data)))

hdpl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.holediam.petleng))]),
                                   data.frame(pars.holediam.petleng,
                                              species = rownames(pars.holediam.petleng)),
                                   names.col = "species")
phylANOVA(hdpl.comp.data$phy, x = setNames(hdpl.comp.data$data$tip_state, rownames(hdpl.comp.data$data)), y = setNames(hdpl.comp.data$data$expected_var, rownames(hdpl.comp.data$data)))

hdsa.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.holediam.stemarea))]),
                                   data.frame(pars.holediam.stemarea,
                                              species = rownames(pars.holediam.stemarea)),
                                   names.col = "species")
phylANOVA(hdsa.comp.data$phy, x = setNames(hdsa.comp.data$data$tip_state, rownames(hdsa.comp.data$data)), y = setNames(hdsa.comp.data$data$expected_var, rownames(hdsa.comp.data$data)))


### Reward

rewcl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.reward.corleng))]),
                                   data.frame(pars.reward.corleng,
                                              species = rownames(pars.reward.corleng)),
                                   names.col = "species")
summary(pgls(expected_var ~ tip_state, data = rewcl.comp.data, lambda = "ML"))

rewla.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.reward.leafarea))]),
                                   data.frame(pars.reward.leafarea,
                                              species = rownames(pars.reward.leafarea)),
                                   names.col = "species")
summary(pgls(expected_var ~ tip_state, data = rewla.comp.data, lambda = "ML"))

rewpl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.reward.petleng))]),
                                   data.frame(pars.reward.petleng,
                                              species = rownames(pars.reward.petleng)),
                                   names.col = "species")
summary(pgls(expected_var ~ tip_state, data = rewpl.comp.data, lambda = "ML"))

rewsa.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.reward.stemarea))]),
                                   data.frame(pars.reward.stemarea,
                                              species = rownames(pars.reward.stemarea)),
                                   names.col = "species")
summary(pgls(expected_var ~ tip_state, data = rewsa.comp.data, lambda = "ML"))


### Warts

warcl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.warts.corleng))]),
                                   data.frame(pars.warts.corleng,
                                              species = rownames(pars.warts.corleng)),
                                   names.col = "species")
phylANOVA(warcl.comp.data$phy, x = setNames(warcl.comp.data$data$tip_state, rownames(warcl.comp.data$data)), y = setNames(warcl.comp.data$data$expected_var, rownames(warcl.comp.data$data)))

warla.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.warts.leafarea))]),
                                   data.frame(pars.warts.leafarea,
                                              species = rownames(pars.warts.leafarea)),
                                   names.col = "species")
phylANOVA(warla.comp.data$phy, x = setNames(warla.comp.data$data$tip_state, rownames(warla.comp.data$data)), y = setNames(warla.comp.data$data$expected_var, rownames(warla.comp.data$data)))

warpl.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.warts.petleng))]),
                                   data.frame(pars.warts.petleng,
                                              species = rownames(pars.warts.petleng)),
                                   names.col = "species")
phylANOVA(warpl.comp.data$phy, x = setNames(warpl.comp.data$data$tip_state, rownames(warpl.comp.data$data)), y = setNames(warpl.comp.data$data$expected_var, rownames(warpl.comp.data$data)))

warsa.comp.data <- comparative.data(drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% rownames(pars.warts.stemarea))]),
                                   data.frame(pars.warts.stemarea,
                                              species = rownames(pars.warts.stemarea)),
                                   names.col = "species")
phylANOVA(warsa.comp.data$phy, x = setNames(warsa.comp.data$data$tip_state, rownames(warsa.comp.data$data)), y = setNames(warsa.comp.data$data$expected_var, rownames(warsa.comp.data$data)))
