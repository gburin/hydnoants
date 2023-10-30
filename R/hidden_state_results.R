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

i_am("R/hidden_state_results.R")

for(i in 1:length(list.files(here("output/houwie/hidden_state"), pattern = ".RDS"))){
    assign(gsub("-", ".", gsub("-hid.RDS", ".hid", list.files(here("output/houwie/hidden_state"), pattern = ".RDS")[i])),
           readRDS(paste0("output/houwie/hidden_state/", list.files(here("output/houwie/hidden_state"), pattern = ".RDS"))[i])
           )
}

for(i in 1:length(list.files(here("output/houwie/main_analysis/"), pattern = ".RDS"))){
    assign(gsub("-", ".", gsub("-nohid.RDS", ".nohid", list.files(here("output/houwie/main_analysis/"), pattern = ".RDS")[i])),
           readRDS(paste0("output/houwie/main_analysis/", list.files(here("output/houwie/main_analysis/"), pattern = ".RDS"))[i])
           )
}

## Domatium Growth

### Corola Length

model.list.domgrow.corleng <- list(er.bm1.hid = domgrow.corleng.er.hid[[1]],
     er.bms.hid = domgrow.corleng.er.hid[[2]],
     er.ou1.hid = domgrow.corleng.er.hid[[3]],
     er.oum.hid = domgrow.corleng.er.hid[[4]],
     er.ouma.hid = domgrow.corleng.er.hid[[5]],
     er.oumv.hid = domgrow.corleng.er.hid[[6]],
     sym.bm1.hid = domgrow.corleng.sym.hid[[1]],
     sym.bms.hid = domgrow.corleng.sym.hid[[2]],
     sym.ou1.hid = domgrow.corleng.sym.hid[[3]],
     sym.oum.hid = domgrow.corleng.sym.hid[[4]],
     sym.ouma.hid = domgrow.corleng.sym.hid[[5]],
     sym.oumv.hid = domgrow.corleng.sym.hid[[6]],
     ard.bm1.hid = domgrow.corleng.ard.hid[[1]],
     ard.bms.hid = domgrow.corleng.ard.hid[[2]],
     ard.ou1.hid = domgrow.corleng.ard.hid[[3]],
     ard.oum.hid = domgrow.corleng.ard.hid[[4]],
     ard.ouma.hid = domgrow.corleng.ard.hid[[5]],
     ard.oumv.hid = domgrow.corleng.ard.hid[[6]],
     er.bm1 = domgrow.corleng.er.nohid[[1]],
     er.bms = domgrow.corleng.er.nohid[[2]],
     er.ou1 = domgrow.corleng.er.nohid[[3]],
     er.oum = domgrow.corleng.er.nohid[[4]],
     er.ouma = domgrow.corleng.er.nohid[[5]],
     er.oumv = domgrow.corleng.er.nohid[[6]],
     sym.bm1 = domgrow.corleng.sym.nohid[[1]],
     sym.bms = domgrow.corleng.sym.nohid[[2]],
     sym.ou1 = domgrow.corleng.sym.nohid[[3]],
     sym.oum = domgrow.corleng.sym.nohid[[4]],
     sym.ouma = domgrow.corleng.sym.nohid[[5]],
     sym.oumv = domgrow.corleng.sym.nohid[[6]],
     ard.bm1 = domgrow.corleng.ard.nohid[[1]],
     ard.bms = domgrow.corleng.ard.nohid[[2]],
     ard.ou1 = domgrow.corleng.ard.nohid[[3]],
     ard.oum = domgrow.corleng.ard.nohid[[4]],
     ard.ouma = domgrow.corleng.ard.nohid[[5]],
     ard.oumv = domgrow.corleng.ard.nohid[[6]])

bic.index.dgcl <- sapply(model.list.domgrow.corleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
domgrow.corleng.table.hid <- getModelTable(model.list.domgrow.corleng[-which(abs(bic.index.dgcl) > 1e04)])

if(any(abs(bic.index.dgcl) > 1e04, na.rm = TRUE)){
    pars.domgrow.corleng <- getModelAvgParams(model.list.domgrow.corleng[-which(abs(bic.index.dgcl) > 1e04)], force = FALSE)
} else {                                                                              pars.domgrow.corleng <- getModelAvgParams(model.list.domgrow.corleng, force = FALSE)}

plot.data.hid.domgrow.corleng <- reshape2::melt(pars.domgrow.corleng)
plot.data.hid.domgrow.corleng$discrete <- "Domatium Growth"
plot.data.hid.domgrow.corleng$continuous <- "Corolla Length"
plot.data.hid.domgrow.corleng$tip_state <- factor(c("Apical", "Diffuse")[as.integer(plot.data.hid.domgrow.corleng$tip_state)], levels = c("Apical", "Diffuse"))

## ggplot(plot.data.hid.domgrow.corleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Leaf Area

model.list.domgrow.leafarea <- list(er.bm1.hid = domgrow.leafarea.er.hid[[1]],
     er.bms.hid = domgrow.leafarea.er.hid[[2]],
     er.ou1.hid = domgrow.leafarea.er.hid[[3]],
     er.oum.hid = domgrow.leafarea.er.hid[[4]],
     er.ouma.hid = domgrow.leafarea.er.hid[[5]],
     er.oumv.hid = domgrow.leafarea.er.hid[[6]],
     sym.bm1.hid = domgrow.leafarea.sym.hid[[1]],
     sym.bms.hid = domgrow.leafarea.sym.hid[[2]],
     sym.ou1.hid = domgrow.leafarea.sym.hid[[3]],
     sym.oum.hid = domgrow.leafarea.sym.hid[[4]],
     sym.ouma.hid = domgrow.leafarea.sym.hid[[5]],
     sym.oumv.hid = domgrow.leafarea.sym.hid[[6]],
     ard.bm1.hid = domgrow.leafarea.ard.hid[[1]],
     ard.bms.hid = domgrow.leafarea.ard.hid[[2]],
     ard.ou1.hid = domgrow.leafarea.ard.hid[[3]],
     ard.oum.hid = domgrow.leafarea.ard.hid[[4]],
     ard.ouma.hid = domgrow.leafarea.ard.hid[[5]],
     ard.oumv.hid = domgrow.leafarea.ard.hid[[6]],
     er.bm1 = domgrow.leafarea.er.nohid[[1]],
     er.bms = domgrow.leafarea.er.nohid[[2]],
     er.ou1 = domgrow.leafarea.er.nohid[[3]],
     er.oum = domgrow.leafarea.er.nohid[[4]],
     er.ouma = domgrow.leafarea.er.nohid[[5]],
     er.oumv = domgrow.leafarea.er.nohid[[6]],
     sym.bm1 = domgrow.leafarea.sym.nohid[[1]],
     sym.bms = domgrow.leafarea.sym.nohid[[2]],
     sym.ou1 = domgrow.leafarea.sym.nohid[[3]],
     sym.oum = domgrow.leafarea.sym.nohid[[4]],
     sym.ouma = domgrow.leafarea.sym.nohid[[5]],
     sym.oumv = domgrow.leafarea.sym.nohid[[6]],
     ard.bm1 = domgrow.leafarea.ard.nohid[[1]],
     ard.bms = domgrow.leafarea.ard.nohid[[2]],
     ard.ou1 = domgrow.leafarea.ard.nohid[[3]],
     ard.oum = domgrow.leafarea.ard.nohid[[4]],
     ard.ouma = domgrow.leafarea.ard.nohid[[5]],
     ard.oumv = domgrow.leafarea.ard.nohid[[6]])

bic.index.dgla <- sapply(model.list.domgrow.leafarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
domgrow.leafarea.table.hid <- getModelTable(model.list.domgrow.leafarea[-which(abs(bic.index.dgla) > 1e04)])

if(any(abs(bic.index.dgla) > 1e04, na.rm = TRUE)){
    pars.domgrow.leafarea <- getModelAvgParams(model.list.domgrow.leafarea[-which(abs(bic.index.dgla) > 1e04)], force = FALSE)
} else {                                                                              pars.domgrow.leafarea <- getModelAvgParams(model.list.domgrow.leafarea, force = FALSE)}

plot.data.hid.domgrow.leafarea <- reshape2::melt(pars.domgrow.leafarea)
plot.data.hid.domgrow.leafarea$discrete <- "Domatium Growth"
plot.data.hid.domgrow.leafarea$continuous <- "Leaf Area"
plot.data.hid.domgrow.leafarea$tip_state <- factor(c("Apical", "Diffuse")[as.integer(plot.data.hid.domgrow.leafarea$tip_state)], levels = c("Apical", "Diffuse"))

## ggplot(plot.data.hid.domgrow.leafarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Petiole Length

model.list.domgrow.petleng <- list(er.bm1.hid = domgrow.petleng.er.hid[[1]],
     er.bms.hid = domgrow.petleng.er.hid[[2]],
     er.ou1.hid = domgrow.petleng.er.hid[[3]],
     er.oum.hid = domgrow.petleng.er.hid[[4]],
     er.ouma.hid = domgrow.petleng.er.hid[[5]],
     er.oumv.hid = domgrow.petleng.er.hid[[6]],
     sym.bm1.hid = domgrow.petleng.sym.hid[[1]],
     sym.bms.hid = domgrow.petleng.sym.hid[[2]],
     sym.ou1.hid = domgrow.petleng.sym.hid[[3]],
     sym.oum.hid = domgrow.petleng.sym.hid[[4]],
     sym.ouma.hid = domgrow.petleng.sym.hid[[5]],
     sym.oumv.hid = domgrow.petleng.sym.hid[[6]],
     ard.bm1.hid = domgrow.petleng.ard.hid[[1]],
     ard.bms.hid = domgrow.petleng.ard.hid[[2]],
     ard.ou1.hid = domgrow.petleng.ard.hid[[3]],
     ard.oum.hid = domgrow.petleng.ard.hid[[4]],
     ard.ouma.hid = domgrow.petleng.ard.hid[[5]],
     ard.oumv.hid = domgrow.petleng.ard.hid[[6]],
     er.bm1 = domgrow.petleng.er.nohid[[1]],
     er.bms = domgrow.petleng.er.nohid[[2]],
     er.ou1 = domgrow.petleng.er.nohid[[3]],
     er.oum = domgrow.petleng.er.nohid[[4]],
     er.ouma = domgrow.petleng.er.nohid[[5]],
     er.oumv = domgrow.petleng.er.nohid[[6]],
     sym.bm1 = domgrow.petleng.sym.nohid[[1]],
     sym.bms = domgrow.petleng.sym.nohid[[2]],
     sym.ou1 = domgrow.petleng.sym.nohid[[3]],
     sym.oum = domgrow.petleng.sym.nohid[[4]],
     sym.ouma = domgrow.petleng.sym.nohid[[5]],
     sym.oumv = domgrow.petleng.sym.nohid[[6]],
     ard.bm1 = domgrow.petleng.ard.nohid[[1]],
     ard.bms = domgrow.petleng.ard.nohid[[2]],
     ard.ou1 = domgrow.petleng.ard.nohid[[3]],
     ard.oum = domgrow.petleng.ard.nohid[[4]],
     ard.ouma = domgrow.petleng.ard.nohid[[5]],
     ard.oumv = domgrow.petleng.ard.nohid[[6]])

bic.index.dgpl <- sapply(model.list.domgrow.petleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
domgrow.petleng.table.hid <- getModelTable(model.list.domgrow.petleng[-which(abs(bic.index.dgpl) > 1e04)])

if(any(abs(bic.index.dgpl) > 1e04, na.rm = TRUE)){
    pars.domgrow.petleng <- getModelAvgParams(model.list.domgrow.petleng[-which(abs(bic.index.dgpl) > 1e04)], force = FALSE)
} else {                                                                              pars.domgrow.petleng <- getModelAvgParams(model.list.domgrow.petleng, force = FALSE)}

plot.data.hid.domgrow.petleng <- reshape2::melt(pars.domgrow.petleng)
plot.data.hid.domgrow.petleng$discrete <- "Domatium Growth"
plot.data.hid.domgrow.petleng$continuous <- "Petiole Length"
plot.data.hid.domgrow.petleng$tip_state <- factor(c("Apical", "Diffuse")[as.integer(plot.data.hid.domgrow.petleng$tip_state)], levels = c("Apical", "Diffuse"))

## ggplot(plot.data.hid.domgrow.petleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")

### Stem Area

model.list.domgrow.stemarea <- list(er.bm1.hid = domgrow.stemarea.er.hid[[1]],
     er.bms.hid = domgrow.stemarea.er.hid[[2]],
     er.ou1.hid = domgrow.stemarea.er.hid[[3]],
     er.oum.hid = domgrow.stemarea.er.hid[[4]],
     er.ouma.hid = domgrow.stemarea.er.hid[[5]],
     er.oumv.hid = domgrow.stemarea.er.hid[[6]],
     sym.bm1.hid = domgrow.stemarea.sym.hid[[1]],
     sym.bms.hid = domgrow.stemarea.sym.hid[[2]],
     sym.ou1.hid = domgrow.stemarea.sym.hid[[3]],
     sym.oum.hid = domgrow.stemarea.sym.hid[[4]],
     sym.ouma.hid = domgrow.stemarea.sym.hid[[5]],
     sym.oumv.hid = domgrow.stemarea.sym.hid[[6]],
     ard.bm1.hid = domgrow.stemarea.ard.hid[[1]],
     ard.bms.hid = domgrow.stemarea.ard.hid[[2]],
     ard.ou1.hid = domgrow.stemarea.ard.hid[[3]],
     ard.oum.hid = domgrow.stemarea.ard.hid[[4]],
     ard.ouma.hid = domgrow.stemarea.ard.hid[[5]],
     ard.oumv.hid = domgrow.stemarea.ard.hid[[6]],
     er.bm1 = domgrow.stemarea.er.nohid[[1]],
     er.bms = domgrow.stemarea.er.nohid[[2]],
     er.ou1 = domgrow.stemarea.er.nohid[[3]],
     er.oum = domgrow.stemarea.er.nohid[[4]],
     er.ouma = domgrow.stemarea.er.nohid[[5]],
     er.oumv = domgrow.stemarea.er.nohid[[6]],
     sym.bm1 = domgrow.stemarea.sym.nohid[[1]],
     sym.bms = domgrow.stemarea.sym.nohid[[2]],
     sym.ou1 = domgrow.stemarea.sym.nohid[[3]],
     sym.oum = domgrow.stemarea.sym.nohid[[4]],
     sym.ouma = domgrow.stemarea.sym.nohid[[5]],
     sym.oumv = domgrow.stemarea.sym.nohid[[6]],
     ard.bm1 = domgrow.stemarea.ard.nohid[[1]],
     ard.bms = domgrow.stemarea.ard.nohid[[2]],
     ard.ou1 = domgrow.stemarea.ard.nohid[[3]],
     ard.oum = domgrow.stemarea.ard.nohid[[4]],
     ard.ouma = domgrow.stemarea.ard.nohid[[5]],
     ard.oumv = domgrow.stemarea.ard.nohid[[6]])

bic.index.dgsa <- sapply(model.list.domgrow.stemarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
domgrow.stemarea.table.hid <- getModelTable(model.list.domgrow.stemarea[-which(abs(bic.index.dgsa) > 1e04)])

if(any(abs(bic.index.dgsa) > 1e04, na.rm = TRUE)){
    pars.domgrow.stemarea <- getModelAvgParams(model.list.domgrow.stemarea[-which(abs(bic.index.dgsa) > 1e04)], force = FALSE)
} else {                                                                              pars.domgrow.stemarea <- getModelAvgParams(model.list.domgrow.stemarea, force = FALSE)}

plot.data.hid.domgrow.stemarea <- reshape2::melt(pars.domgrow.stemarea)
plot.data.hid.domgrow.stemarea$discrete <- "Domatium Growth"
plot.data.hid.domgrow.stemarea$continuous <- "Stem Area"
plot.data.hid.domgrow.stemarea$tip_state <- factor(c("Apical", "Diffuse")[as.integer(plot.data.hid.domgrow.stemarea$tip_state)], levels = c("Apical", "Diffuse"))

## ggplot(plot.data.hid.domgrow.stemarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")






## Hole Diameter

### Corola Length

model.list.holediam.corleng <- list(er.bm1.hid = holediam.corleng.er.hid[[1]],
     er.bms.hid = holediam.corleng.er.hid[[2]],
     er.ou1.hid = holediam.corleng.er.hid[[3]],
     er.oum.hid = holediam.corleng.er.hid[[4]],
     er.ouma.hid = holediam.corleng.er.hid[[5]],
     er.oumv.hid = holediam.corleng.er.hid[[6]],
     sym.bm1.hid = holediam.corleng.sym.hid[[1]],
     sym.bms.hid = holediam.corleng.sym.hid[[2]],
     sym.ou1.hid = holediam.corleng.sym.hid[[3]],
     sym.oum.hid = holediam.corleng.sym.hid[[4]],
     sym.ouma.hid = holediam.corleng.sym.hid[[5]],
     sym.oumv.hid = holediam.corleng.sym.hid[[6]],
     ard.bm1.hid = holediam.corleng.ard.hid[[1]],
     ard.bms.hid = holediam.corleng.ard.hid[[2]],
     ard.ou1.hid = holediam.corleng.ard.hid[[3]],
     ard.oum.hid = holediam.corleng.ard.hid[[4]],
     ard.ouma.hid = holediam.corleng.ard.hid[[5]],
     ard.oumv.hid = holediam.corleng.ard.hid[[6]],
     er.bm1 = holediam.corleng.er.nohid[[1]],
     er.bms = holediam.corleng.er.nohid[[2]],
     er.ou1 = holediam.corleng.er.nohid[[3]],
     er.oum = holediam.corleng.er.nohid[[4]],
     er.ouma = holediam.corleng.er.nohid[[5]],
     er.oumv = holediam.corleng.er.nohid[[6]],
     sym.bm1 = holediam.corleng.sym.nohid[[1]],
     sym.bms = holediam.corleng.sym.nohid[[2]],
     sym.ou1 = holediam.corleng.sym.nohid[[3]],
     sym.oum = holediam.corleng.sym.nohid[[4]],
     sym.ouma = holediam.corleng.sym.nohid[[5]],
     sym.oumv = holediam.corleng.sym.nohid[[6]],
     ard.bm1 = holediam.corleng.ard.nohid[[1]],
     ard.bms = holediam.corleng.ard.nohid[[2]],
     ard.ou1 = holediam.corleng.ard.nohid[[3]],
     ard.oum = holediam.corleng.ard.nohid[[4]],
     ard.ouma = holediam.corleng.ard.nohid[[5]],
     ard.oumv = holediam.corleng.ard.nohid[[6]])

bic.index.hdcl <- sapply(model.list.holediam.corleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
holediam.corleng.table.hid <- getModelTable(model.list.holediam.corleng[-which(abs(bic.index.hdcl) > 1e04)])

if(any(abs(bic.index.hdcl) > 1e04, na.rm = TRUE)){
    pars.holediam.corleng <- getModelAvgParams(model.list.holediam.corleng[-which(abs(bic.index.hdcl) > 1e04)], force = FALSE)
} else {                                                                              pars.holediam.corleng <- getModelAvgParams(model.list.holediam.corleng, force = FALSE)}

plot.data.hid.holediam.corleng <- reshape2::melt(pars.holediam.corleng)
plot.data.hid.holediam.corleng$discrete <- "Hole Diameter"
plot.data.hid.holediam.corleng$continuous <- "Corolla Length"
plot.data.hid.holediam.corleng$tip_state <- factor(c("Several Large\nat Base", "One Large\nat Base", "All Large")[as.integer(plot.data.hid.holediam.corleng$tip_state)], levels = c("Several Large\nat Base", "One Large\nat Base", "All Large"))

## ggplot(plot.data.hid.holediam.corleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Leaf Area

model.list.holediam.leafarea <- list(er.bm1.hid = holediam.leafarea.er.hid[[1]],
     er.bms.hid = holediam.leafarea.er.hid[[2]],
     er.ou1.hid = holediam.leafarea.er.hid[[3]],
     er.oum.hid = holediam.leafarea.er.hid[[4]],
     er.ouma.hid = holediam.leafarea.er.hid[[5]],
     er.oumv.hid = holediam.leafarea.er.hid[[6]],
     sym.bm1.hid = holediam.leafarea.sym.hid[[1]],
     sym.bms.hid = holediam.leafarea.sym.hid[[2]],
     sym.ou1.hid = holediam.leafarea.sym.hid[[3]],
     sym.oum.hid = holediam.leafarea.sym.hid[[4]],
     sym.ouma.hid = holediam.leafarea.sym.hid[[5]],
     sym.oumv.hid = holediam.leafarea.sym.hid[[6]],
     ard.bm1.hid = holediam.leafarea.ard.hid[[1]],
     ard.bms.hid = holediam.leafarea.ard.hid[[2]],
     ard.ou1.hid = holediam.leafarea.ard.hid[[3]],
     ard.oum.hid = holediam.leafarea.ard.hid[[4]],
     ard.ouma.hid = holediam.leafarea.ard.hid[[5]],
     ard.oumv.hid = holediam.leafarea.ard.hid[[6]],
     er.bm1 = holediam.leafarea.er.nohid[[1]],
     er.bms = holediam.leafarea.er.nohid[[2]],
     er.ou1 = holediam.leafarea.er.nohid[[3]],
     er.oum = holediam.leafarea.er.nohid[[4]],
     er.ouma = holediam.leafarea.er.nohid[[5]],
     er.oumv = holediam.leafarea.er.nohid[[6]],
     sym.bm1 = holediam.leafarea.sym.nohid[[1]],
     sym.bms = holediam.leafarea.sym.nohid[[2]],
     sym.ou1 = holediam.leafarea.sym.nohid[[3]],
     sym.oum = holediam.leafarea.sym.nohid[[4]],
     sym.ouma = holediam.leafarea.sym.nohid[[5]],
     sym.oumv = holediam.leafarea.sym.nohid[[6]],
     ard.bm1 = holediam.leafarea.ard.nohid[[1]],
     ard.bms = holediam.leafarea.ard.nohid[[2]],
     ard.ou1 = holediam.leafarea.ard.nohid[[3]],
     ard.oum = holediam.leafarea.ard.nohid[[4]],
     ard.ouma = holediam.leafarea.ard.nohid[[5]],
     ard.oumv = holediam.leafarea.ard.nohid[[6]])

bic.index.hdla <- sapply(model.list.holediam.leafarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
holediam.leafarea.table.hid <- getModelTable(model.list.holediam.leafarea[-which(abs(bic.index.hdla) > 1e04)])

if(any(abs(bic.index.hdla) > 1e04, na.rm = TRUE)){
    pars.holediam.leafarea <- getModelAvgParams(model.list.holediam.leafarea[-which(abs(bic.index.hdla) > 1e04)], force = FALSE)
} else {                                                                              pars.holediam.leafarea <- getModelAvgParams(model.list.holediam.leafarea, force = FALSE)}

plot.data.hid.holediam.leafarea <- reshape2::melt(pars.holediam.leafarea)
plot.data.hid.holediam.leafarea$discrete <- "Hole Diameter"
plot.data.hid.holediam.leafarea$continuous <- "Leaf Area"
plot.data.hid.holediam.leafarea$tip_state <- factor(c("Several Large\nat Base", "One Large\nat Base", "All Large")[as.integer(plot.data.hid.holediam.leafarea$tip_state)], levels = c("Several Large\nat Base", "One Large\nat Base", "All Large"))

## ggplot(plot.data.hid.holediam.leafarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Petiole Length

model.list.holediam.petleng <- list(er.bm1.hid = holediam.petleng.er.hid[[1]],
     er.bms.hid = holediam.petleng.er.hid[[2]],
     er.ou1.hid = holediam.petleng.er.hid[[3]],
     er.oum.hid = holediam.petleng.er.hid[[4]],
     er.ouma.hid = holediam.petleng.er.hid[[5]],
     er.oumv.hid = holediam.petleng.er.hid[[6]],
     sym.bm1.hid = holediam.petleng.sym.hid[[1]],
     sym.bms.hid = holediam.petleng.sym.hid[[2]],
     sym.ou1.hid = holediam.petleng.sym.hid[[3]],
     sym.oum.hid = holediam.petleng.sym.hid[[4]],
     sym.ouma.hid = holediam.petleng.sym.hid[[5]],
     sym.oumv.hid = holediam.petleng.sym.hid[[6]],
     ard.bm1.hid = holediam.petleng.ard.hid[[1]],
     ard.bms.hid = holediam.petleng.ard.hid[[2]],
     ard.ou1.hid = holediam.petleng.ard.hid[[3]],
     ard.oum.hid = holediam.petleng.ard.hid[[4]],
     ard.ouma.hid = holediam.petleng.ard.hid[[5]],
     ard.oumv.hid = holediam.petleng.ard.hid[[6]],
     er.bm1 = holediam.petleng.er.nohid[[1]],
     er.bms = holediam.petleng.er.nohid[[2]],
     er.ou1 = holediam.petleng.er.nohid[[3]],
     er.oum = holediam.petleng.er.nohid[[4]],
     er.ouma = holediam.petleng.er.nohid[[5]],
     er.oumv = holediam.petleng.er.nohid[[6]],
     sym.bm1 = holediam.petleng.sym.nohid[[1]],
     sym.bms = holediam.petleng.sym.nohid[[2]],
     sym.ou1 = holediam.petleng.sym.nohid[[3]],
     sym.oum = holediam.petleng.sym.nohid[[4]],
     sym.ouma = holediam.petleng.sym.nohid[[5]],
     sym.oumv = holediam.petleng.sym.nohid[[6]],
     ard.bm1 = holediam.petleng.ard.nohid[[1]],
     ard.bms = holediam.petleng.ard.nohid[[2]],
     ard.ou1 = holediam.petleng.ard.nohid[[3]],
     ard.oum = holediam.petleng.ard.nohid[[4]],
     ard.ouma = holediam.petleng.ard.nohid[[5]],
     ard.oumv = holediam.petleng.ard.nohid[[6]])

bic.index.hdpl <- sapply(model.list.holediam.petleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
holediam.petleng.table.hid <- getModelTable(model.list.holediam.petleng[-which(abs(bic.index.hdpl) > 1e04)])

if(any(abs(bic.index.hdpl) > 1e04, na.rm = TRUE)){
    pars.holediam.petleng <- getModelAvgParams(model.list.holediam.petleng[-which(abs(bic.index.hdpl) > 1e04)], force = FALSE)
} else {                                                                              pars.holediam.petleng <- getModelAvgParams(model.list.holediam.petleng, force = FALSE)}

plot.data.hid.holediam.petleng <- reshape2::melt(pars.holediam.petleng)
plot.data.hid.holediam.petleng$discrete <- "Hole Diameter"
plot.data.hid.holediam.petleng$continuous <- "Petiole Length"
plot.data.hid.holediam.petleng$tip_state <- factor(c("Several Large\nat Base", "One Large\nat Base", "All Large")[as.integer(plot.data.hid.holediam.petleng$tip_state)], levels = c("Several Large\nat Base", "One Large\nat Base", "All Large"))

## ggplot(plot.data.hid.holediam.petleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")

### Stem Area

model.list.holediam.stemarea <- list(er.bm1.hid = holediam.stemarea.er.hid[[1]],
     er.bms.hid = holediam.stemarea.er.hid[[2]],
     er.ou1.hid = holediam.stemarea.er.hid[[3]],
     er.oum.hid = holediam.stemarea.er.hid[[4]],
     er.ouma.hid = holediam.stemarea.er.hid[[5]],
     er.oumv.hid = holediam.stemarea.er.hid[[6]],
     sym.bm1.hid = holediam.stemarea.sym.hid[[1]],
     sym.bms.hid = holediam.stemarea.sym.hid[[2]],
     sym.ou1.hid = holediam.stemarea.sym.hid[[3]],
     sym.oum.hid = holediam.stemarea.sym.hid[[4]],
     sym.ouma.hid = holediam.stemarea.sym.hid[[5]],
     sym.oumv.hid = holediam.stemarea.sym.hid[[6]],
     ard.bm1.hid = holediam.stemarea.ard.hid[[1]],
     ard.bms.hid = holediam.stemarea.ard.hid[[2]],
     ard.ou1.hid = holediam.stemarea.ard.hid[[3]],
     ard.oum.hid = holediam.stemarea.ard.hid[[4]],
     ard.ouma.hid = holediam.stemarea.ard.hid[[5]],
     ard.oumv.hid = holediam.stemarea.ard.hid[[6]],
     er.bm1 = holediam.stemarea.er.nohid[[1]],
     er.bms = holediam.stemarea.er.nohid[[2]],
     er.ou1 = holediam.stemarea.er.nohid[[3]],
     er.oum = holediam.stemarea.er.nohid[[4]],
     er.ouma = holediam.stemarea.er.nohid[[5]],
     er.oumv = holediam.stemarea.er.nohid[[6]],
     sym.bm1 = holediam.stemarea.sym.nohid[[1]],
     sym.bms = holediam.stemarea.sym.nohid[[2]],
     sym.ou1 = holediam.stemarea.sym.nohid[[3]],
     sym.oum = holediam.stemarea.sym.nohid[[4]],
     sym.ouma = holediam.stemarea.sym.nohid[[5]],
     sym.oumv = holediam.stemarea.sym.nohid[[6]],
     ard.bm1 = holediam.stemarea.ard.nohid[[1]],
     ard.bms = holediam.stemarea.ard.nohid[[2]],
     ard.ou1 = holediam.stemarea.ard.nohid[[3]],
     ard.oum = holediam.stemarea.ard.nohid[[4]],
     ard.ouma = holediam.stemarea.ard.nohid[[5]],
     ard.oumv = holediam.stemarea.ard.nohid[[6]])

bic.index.hdsa <- sapply(model.list.holediam.stemarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
holediam.stemarea.table.hid <- getModelTable(model.list.holediam.stemarea[-which(abs(bic.index.hdsa) > 1e04)])

if(any(abs(bic.index.hdsa) > 1e04, na.rm = TRUE)){
    pars.holediam.stemarea <- getModelAvgParams(model.list.holediam.stemarea[-which(abs(bic.index.hdsa) > 1e04)], force = FALSE)
} else {                                                                              pars.holediam.stemarea <- getModelAvgParams(model.list.holediam.stemarea, force = FALSE)}

plot.data.hid.holediam.stemarea <- reshape2::melt(pars.holediam.stemarea)
plot.data.hid.holediam.stemarea$discrete <- "Hole Diameter"
plot.data.hid.holediam.stemarea$continuous <- "Stem Area"
plot.data.hid.holediam.stemarea$tip_state <- factor(c("Several Large\nat Base", "One Large\nat Base", "All Large")[as.integer(plot.data.hid.holediam.stemarea$tip_state)], levels = c("Several Large\nat Base", "One Large\nat Base", "All Large"))

## ggplot(plot.data.hid.holediam.stemarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")





## Reward

### Corola Length

model.list.reward.corleng <- list(er.bm1.hid = reward.corleng.er.hid[[1]],
     er.bms.hid = reward.corleng.er.hid[[2]],
     er.ou1.hid = reward.corleng.er.hid[[3]],
     er.oum.hid = reward.corleng.er.hid[[4]],
     er.ouma.hid = reward.corleng.er.hid[[5]],
     er.oumv.hid = reward.corleng.er.hid[[6]],
     sym.bm1.hid = reward.corleng.sym.hid[[1]],
     sym.bms.hid = reward.corleng.sym.hid[[2]],
     sym.ou1.hid = reward.corleng.sym.hid[[3]],
     sym.oum.hid = reward.corleng.sym.hid[[4]],
     sym.ouma.hid = reward.corleng.sym.hid[[5]],
     sym.oumv.hid = reward.corleng.sym.hid[[6]],
     ard.bm1.hid = reward.corleng.ard.hid[[1]],
     ard.bms.hid = reward.corleng.ard.hid[[2]],
     ard.ou1.hid = reward.corleng.ard.hid[[3]],
     ard.oum.hid = reward.corleng.ard.hid[[4]],
     ard.ouma.hid = reward.corleng.ard.hid[[5]],
     ard.oumv.hid = reward.corleng.ard.hid[[6]],
     er.bm1 = reward.corleng.er.nohid[[1]],
     er.bms = reward.corleng.er.nohid[[2]],
     er.ou1 = reward.corleng.er.nohid[[3]],
     er.oum = reward.corleng.er.nohid[[4]],
     er.ouma = reward.corleng.er.nohid[[5]],
     er.oumv = reward.corleng.er.nohid[[6]],
     sym.bm1 = reward.corleng.sym.nohid[[1]],
     sym.bms = reward.corleng.sym.nohid[[2]],
     sym.ou1 = reward.corleng.sym.nohid[[3]],
     sym.oum = reward.corleng.sym.nohid[[4]],
     sym.ouma = reward.corleng.sym.nohid[[5]],
     ## sym.oumv = reward.corleng.sym.nohid[[6]],
     ard.bm1 = reward.corleng.ard.nohid[[1]],
     ard.bms = reward.corleng.ard.nohid[[2]],
     ard.ou1 = reward.corleng.ard.nohid[[3]],
     ard.oum = reward.corleng.ard.nohid[[4]],
     ard.ouma = reward.corleng.ard.nohid[[5]],
     ard.oumv = reward.corleng.ard.nohid[[6]])

bic.index.recl <- sapply(model.list.reward.corleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
reward.corleng.table.hid <- getModelTable(model.list.reward.corleng[-which(abs(bic.index.recl) > 1e04)])

if(any(abs(bic.index.recl) > 1e04, na.rm = TRUE)){
    pars.reward.corleng <- getModelAvgParams(model.list.reward.corleng[-which(abs(bic.index.recl) > 1e04)], force = FALSE)
} else {                                                                              pars.reward.corleng <- getModelAvgParams(model.list.reward.corleng, force = FALSE)}

plot.data.hid.reward.corleng <- reshape2::melt(pars.reward.corleng)
plot.data.hid.reward.corleng$discrete <- "Reward"
plot.data.hid.reward.corleng$continuous <- "Corolla Length"
plot.data.hid.reward.corleng$tip_state <- factor(c("Absent", "Present")[as.integer(plot.data.hid.reward.corleng$tip_state) + 1], levels = c("Absent", "Present"))

## ggplot(plot.data.hid.reward.corleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Leaf Area

model.list.reward.leafarea <- list(er.bm1.hid = reward.leafarea.er.hid[[1]],
     er.bms.hid = reward.leafarea.er.hid[[2]],
     er.ou1.hid = reward.leafarea.er.hid[[3]],
     er.oum.hid = reward.leafarea.er.hid[[4]],
     er.ouma.hid = reward.leafarea.er.hid[[5]],
     er.oumv.hid = reward.leafarea.er.hid[[6]],
     sym.bm1.hid = reward.leafarea.sym.hid[[1]],
     sym.bms.hid = reward.leafarea.sym.hid[[2]],
     sym.ou1.hid = reward.leafarea.sym.hid[[3]],
     sym.oum.hid = reward.leafarea.sym.hid[[4]],
     sym.ouma.hid = reward.leafarea.sym.hid[[5]],
     sym.oumv.hid = reward.leafarea.sym.hid[[6]],
     ard.bm1.hid = reward.leafarea.ard.hid[[1]],
     ard.bms.hid = reward.leafarea.ard.hid[[2]],
     ard.ou1.hid = reward.leafarea.ard.hid[[3]],
     ard.oum.hid = reward.leafarea.ard.hid[[4]],
     ard.ouma.hid = reward.leafarea.ard.hid[[5]],
     ard.oumv.hid = reward.leafarea.ard.hid[[6]],
     er.bm1 = reward.leafarea.er.nohid[[1]],
     er.bms = reward.leafarea.er.nohid[[2]],
     er.ou1 = reward.leafarea.er.nohid[[3]],
     er.oum = reward.leafarea.er.nohid[[4]],
     er.ouma = reward.leafarea.er.nohid[[5]],
     er.oumv = reward.leafarea.er.nohid[[6]],
     sym.bm1 = reward.leafarea.sym.nohid[[1]],
     sym.bms = reward.leafarea.sym.nohid[[2]],
     sym.ou1 = reward.leafarea.sym.nohid[[3]],
     sym.oum = reward.leafarea.sym.nohid[[4]],
     sym.ouma = reward.leafarea.sym.nohid[[5]],
     sym.oumv = reward.leafarea.sym.nohid[[6]],
     ard.bm1 = reward.leafarea.ard.nohid[[1]],
     ard.bms = reward.leafarea.ard.nohid[[2]],
     ard.ou1 = reward.leafarea.ard.nohid[[3]],
     ard.oum = reward.leafarea.ard.nohid[[4]],
     ard.ouma = reward.leafarea.ard.nohid[[5]],
     ard.oumv = reward.leafarea.ard.nohid[[6]])

bic.index.rela <- sapply(model.list.reward.leafarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
reward.leafarea.table.hid <- getModelTable(model.list.reward.leafarea[-which(abs(bic.index.rela) > 1e04)])

if(any(abs(bic.index.rela) > 1e04, na.rm = TRUE)){
    pars.reward.leafarea <- getModelAvgParams(model.list.reward.leafarea[-which(abs(bic.index.rela) > 1e04)], force = FALSE)
} else {                                                                              pars.reward.leafarea <- getModelAvgParams(model.list.reward.leafarea, force = FALSE)}

plot.data.hid.reward.leafarea <- reshape2::melt(pars.reward.leafarea)
plot.data.hid.reward.leafarea$discrete <- "Reward"
plot.data.hid.reward.leafarea$continuous <- "Leaf Area"
plot.data.hid.reward.leafarea$tip_state <- factor(c("Absent", "Present")[as.integer(plot.data.hid.reward.leafarea$tip_state) + 1], levels = c("Absent", "Present"))

## ggplot(plot.data.hid.reward.leafarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Petiole Length

model.list.reward.petleng <- list(er.bm1.hid = reward.petleng.er.hid[[1]],
     er.bms.hid = reward.petleng.er.hid[[2]],
     er.ou1.hid = reward.petleng.er.hid[[3]],
     er.oum.hid = reward.petleng.er.hid[[4]],
     er.ouma.hid = reward.petleng.er.hid[[5]],
     er.oumv.hid = reward.petleng.er.hid[[6]],
     sym.bm1.hid = reward.petleng.sym.hid[[1]],
     sym.bms.hid = reward.petleng.sym.hid[[2]],
     sym.ou1.hid = reward.petleng.sym.hid[[3]],
     sym.oum.hid = reward.petleng.sym.hid[[4]],
     sym.ouma.hid = reward.petleng.sym.hid[[5]],
     #sym.oumv.hid = reward.petleng.sym.hid[[6]],
     ard.bm1.hid = reward.petleng.ard.hid[[1]],
     ard.bms.hid = reward.petleng.ard.hid[[2]],
     ard.ou1.hid = reward.petleng.ard.hid[[3]],
     ard.oum.hid = reward.petleng.ard.hid[[4]],
     ard.ouma.hid = reward.petleng.ard.hid[[5]],
     ard.oumv.hid = reward.petleng.ard.hid[[6]],
     er.bm1 = reward.petleng.er.nohid[[1]],
     er.bms = reward.petleng.er.nohid[[2]],
     er.ou1 = reward.petleng.er.nohid[[3]],
     er.oum = reward.petleng.er.nohid[[4]],
     er.ouma = reward.petleng.er.nohid[[5]],
     er.oumv = reward.petleng.er.nohid[[6]],
     sym.bm1 = reward.petleng.sym.nohid[[1]],
     sym.bms = reward.petleng.sym.nohid[[2]],
     sym.ou1 = reward.petleng.sym.nohid[[3]],
     sym.oum = reward.petleng.sym.nohid[[4]],
     sym.ouma = reward.petleng.sym.nohid[[5]],
     sym.oumv = reward.petleng.sym.nohid[[6]],
     ard.bm1 = reward.petleng.ard.nohid[[1]],
     ard.bms = reward.petleng.ard.nohid[[2]],
     ard.ou1 = reward.petleng.ard.nohid[[3]],
     ard.oum = reward.petleng.ard.nohid[[4]],
     ard.ouma = reward.petleng.ard.nohid[[5]],
     ard.oumv = reward.petleng.ard.nohid[[6]])

bic.index.repl <- sapply(model.list.reward.petleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
reward.petleng.table.hid <- getModelTable(model.list.reward.petleng[-which(abs(bic.index.repl) > 1e04)])

if(any(abs(bic.index.repl) > 1e04, na.rm = TRUE)){
    pars.reward.petleng <- getModelAvgParams(model.list.reward.petleng[-which(abs(bic.index.repl) > 1e04)], force = FALSE)
} else {                                                                              pars.reward.petleng <- getModelAvgParams(model.list.reward.petleng, force = FALSE)}

plot.data.hid.reward.petleng <- reshape2::melt(pars.reward.petleng)
plot.data.hid.reward.petleng$discrete <- "Reward"
plot.data.hid.reward.petleng$continuous <- "Petiole Length"
plot.data.hid.reward.petleng$tip_state <- factor(c("Absent", "Present")[as.integer(plot.data.hid.reward.petleng$tip_state) + 1], levels = c("Absent", "Present"))

## ggplot(plot.data.hid.reward.petleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")

### Stem Area

model.list.reward.stemarea <- list(er.bm1.hid = reward.stemarea.er.hid[[1]],
     er.bms.hid = reward.stemarea.er.hid[[2]],
     er.ou1.hid = reward.stemarea.er.hid[[3]],
     er.oum.hid = reward.stemarea.er.hid[[4]],
     er.ouma.hid = reward.stemarea.er.hid[[5]],
     er.oumv.hid = reward.stemarea.er.hid[[6]],
     sym.bm1.hid = reward.stemarea.sym.hid[[1]],
     sym.bms.hid = reward.stemarea.sym.hid[[2]],
     sym.ou1.hid = reward.stemarea.sym.hid[[3]],
     sym.oum.hid = reward.stemarea.sym.hid[[4]],
     sym.ouma.hid = reward.stemarea.sym.hid[[5]],
     sym.oumv.hid = reward.stemarea.sym.hid[[6]],
     ard.bm1.hid = reward.stemarea.ard.hid[[1]],
     ard.bms.hid = reward.stemarea.ard.hid[[2]],
     ard.ou1.hid = reward.stemarea.ard.hid[[3]],
     ard.oum.hid = reward.stemarea.ard.hid[[4]],
     ard.ouma.hid = reward.stemarea.ard.hid[[5]],
     ard.oumv.hid = reward.stemarea.ard.hid[[6]],
     er.bm1 = reward.stemarea.er.nohid[[1]],
     er.bms = reward.stemarea.er.nohid[[2]],
     er.ou1 = reward.stemarea.er.nohid[[3]],
     er.oum = reward.stemarea.er.nohid[[4]],
     er.ouma = reward.stemarea.er.nohid[[5]],
     er.oumv = reward.stemarea.er.nohid[[6]],
     sym.bm1 = reward.stemarea.sym.nohid[[1]],
     sym.bms = reward.stemarea.sym.nohid[[2]],
     sym.ou1 = reward.stemarea.sym.nohid[[3]],
     sym.oum = reward.stemarea.sym.nohid[[4]],
     sym.ouma = reward.stemarea.sym.nohid[[5]],
     sym.oumv = reward.stemarea.sym.nohid[[6]],
     ard.bm1 = reward.stemarea.ard.nohid[[1]],
     ard.bms = reward.stemarea.ard.nohid[[2]],
     ard.ou1 = reward.stemarea.ard.nohid[[3]],
     ard.oum = reward.stemarea.ard.nohid[[4]],
     ard.ouma = reward.stemarea.ard.nohid[[5]],
     ard.oumv = reward.stemarea.ard.nohid[[6]])

bic.index.resa <- sapply(model.list.reward.stemarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
reward.stemarea.table.hid <- getModelTable(model.list.reward.stemarea[-which(abs(bic.index.resa) > 1e04)])

if(any(abs(bic.index.resa) > 1e04, na.rm = TRUE)){
    pars.reward.stemarea <- getModelAvgParams(model.list.reward.stemarea[-which(abs(bic.index.resa) > 1e04)], force = FALSE)
} else {                                                                              pars.reward.stemarea <- getModelAvgParams(model.list.reward.stemarea, force = FALSE)}

plot.data.hid.reward.stemarea <- reshape2::melt(pars.reward.stemarea)
plot.data.hid.reward.stemarea$discrete <- "Reward"
plot.data.hid.reward.stemarea$continuous <- "Stem Area"
plot.data.hid.reward.stemarea$tip_state <- factor(c("Absent", "Present")[as.integer(plot.data.hid.reward.stemarea$tip_state) + 1], levels = c("Absent", "Present"))

## ggplot(plot.data.hid.reward.stemarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")




## Warts

### Corola Length

model.list.warts.corleng <- list(er.bm1.hid = warts.corleng.er.hid[[1]],
     er.bms.hid = warts.corleng.er.hid[[2]],
     er.ou1.hid = warts.corleng.er.hid[[3]],
     er.oum.hid = warts.corleng.er.hid[[4]],
     er.ouma.hid = warts.corleng.er.hid[[5]],
     er.oumv.hid = warts.corleng.er.hid[[6]],
     sym.bm1.hid = warts.corleng.sym.hid[[1]],
     sym.bms.hid = warts.corleng.sym.hid[[2]],
     sym.ou1.hid = warts.corleng.sym.hid[[3]],
     sym.oum.hid = warts.corleng.sym.hid[[4]],
     sym.ouma.hid = warts.corleng.sym.hid[[5]],
     sym.oumv.hid = warts.corleng.sym.hid[[6]],
     ard.bm1.hid = warts.corleng.ard.hid[[1]],
     ard.bms.hid = warts.corleng.ard.hid[[2]],
     ard.ou1.hid = warts.corleng.ard.hid[[3]],
     ard.oum.hid = warts.corleng.ard.hid[[4]],
     ard.ouma.hid = warts.corleng.ard.hid[[5]],
     ard.oumv.hid = warts.corleng.ard.hid[[6]],
     er.bm1 = warts.corleng.er.nohid[[1]],
     er.bms = warts.corleng.er.nohid[[2]],
     er.ou1 = warts.corleng.er.nohid[[3]],
     er.oum = warts.corleng.er.nohid[[4]],
     er.ouma = warts.corleng.er.nohid[[5]],
     er.oumv = warts.corleng.er.nohid[[6]],
     sym.bm1 = warts.corleng.sym.nohid[[1]],
     sym.bms = warts.corleng.sym.nohid[[2]],
     sym.ou1 = warts.corleng.sym.nohid[[3]],
     sym.oum = warts.corleng.sym.nohid[[4]],
     sym.ouma = warts.corleng.sym.nohid[[5]],
     sym.oumv = warts.corleng.sym.nohid[[6]],
     ard.bm1 = warts.corleng.ard.nohid[[1]],
     ard.bms = warts.corleng.ard.nohid[[2]],
     ard.ou1 = warts.corleng.ard.nohid[[3]],
     ard.oum = warts.corleng.ard.nohid[[4]],
     ard.ouma = warts.corleng.ard.nohid[[5]],
     ard.oumv = warts.corleng.ard.nohid[[6]])

bic.index.wacl <- sapply(model.list.warts.corleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
warts.corleng.table.hid <- getModelTable(model.list.warts.corleng[-which(abs(bic.index.wacl) > 1e04)])

if(any(abs(bic.index.wacl) > 1e04, na.rm = TRUE)){
    pars.warts.corleng <- getModelAvgParams(model.list.warts.corleng[-which(abs(bic.index.wacl) > 1e04)], force = FALSE)
} else {                                                                              pars.warts.corleng <- getModelAvgParams(model.list.warts.corleng, force = FALSE)}

plot.data.hid.warts.corleng <- reshape2::melt(pars.warts.corleng)
plot.data.hid.warts.corleng$discrete <- "Warts"
plot.data.hid.warts.corleng$continuous <- "Corolla Length"
plot.data.hid.warts.corleng$tip_state <- factor(c("Variable", "Differentiated", "Lost")[as.integer(plot.data.hid.warts.corleng$tip_state)], levels = c("Variable", "Differentiated", "Lost"))

## ggplot(plot.data.hid.warts.corleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Leaf Area

model.list.warts.leafarea <- list(er.bm1.hid = warts.leafarea.er.hid[[1]],
     er.bms.hid = warts.leafarea.er.hid[[2]],
     er.ou1.hid = warts.leafarea.er.hid[[3]],
     er.oum.hid = warts.leafarea.er.hid[[4]],
     er.ouma.hid = warts.leafarea.er.hid[[5]],
     er.oumv.hid = warts.leafarea.er.hid[[6]],
     sym.bm1.hid = warts.leafarea.sym.hid[[1]],
     sym.bms.hid = warts.leafarea.sym.hid[[2]],
     sym.ou1.hid = warts.leafarea.sym.hid[[3]],
     sym.oum.hid = warts.leafarea.sym.hid[[4]],
     sym.ouma.hid = warts.leafarea.sym.hid[[5]],
     sym.oumv.hid = warts.leafarea.sym.hid[[6]],
     ard.bm1.hid = warts.leafarea.ard.hid[[1]],
     ard.bms.hid = warts.leafarea.ard.hid[[2]],
     ard.ou1.hid = warts.leafarea.ard.hid[[3]],
     ard.oum.hid = warts.leafarea.ard.hid[[4]],
     ard.ouma.hid = warts.leafarea.ard.hid[[5]],
     ard.oumv.hid = warts.leafarea.ard.hid[[6]],
     er.bm1 = warts.leafarea.er.nohid[[1]],
     er.bms = warts.leafarea.er.nohid[[2]],
     er.ou1 = warts.leafarea.er.nohid[[3]],
     er.oum = warts.leafarea.er.nohid[[4]],
     er.ouma = warts.leafarea.er.nohid[[5]],
     er.oumv = warts.leafarea.er.nohid[[6]],
     sym.bm1 = warts.leafarea.sym.nohid[[1]],
     sym.bms = warts.leafarea.sym.nohid[[2]],
     sym.ou1 = warts.leafarea.sym.nohid[[3]],
     sym.oum = warts.leafarea.sym.nohid[[4]],
     sym.ouma = warts.leafarea.sym.nohid[[5]],
     sym.oumv = warts.leafarea.sym.nohid[[6]],
     ard.bm1 = warts.leafarea.ard.nohid[[1]],
     ard.bms = warts.leafarea.ard.nohid[[2]],
     ard.ou1 = warts.leafarea.ard.nohid[[3]],
     ard.oum = warts.leafarea.ard.nohid[[4]],
     ard.ouma = warts.leafarea.ard.nohid[[5]],
     ard.oumv = warts.leafarea.ard.nohid[[6]])

bic.index.wala <- sapply(model.list.warts.leafarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
warts.leafarea.table.hid <- getModelTable(model.list.warts.leafarea[-which(abs(bic.index.wala) > 1e04)])

if(any(abs(bic.index.wala) > 1e04, na.rm = TRUE)){
    pars.warts.leafarea <- getModelAvgParams(model.list.warts.leafarea[-which(abs(bic.index.wala) > 1e04)], force = FALSE)
} else {                                                                              pars.warts.leafarea <- getModelAvgParams(model.list.warts.leafarea, force = FALSE)}

plot.data.hid.warts.leafarea <- reshape2::melt(pars.warts.leafarea)
plot.data.hid.warts.leafarea$discrete <- "Warts"
plot.data.hid.warts.leafarea$continuous <- "Leaf Area"
plot.data.hid.warts.leafarea$tip_state <- factor(c("Variable", "Differentiated", "Lost")[as.integer(plot.data.hid.warts.leafarea$tip_state)], levels = c("Variable", "Differentiated", "Lost"))

## ggplot(plot.data.hid.warts.leafarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")


### Petiole Length

model.list.warts.petleng <- list(er.bm1.hid = warts.petleng.er.hid[[1]],
     er.bms.hid = warts.petleng.er.hid[[2]],
     er.ou1.hid = warts.petleng.er.hid[[3]],
     er.oum.hid = warts.petleng.er.hid[[4]],
     er.ouma.hid = warts.petleng.er.hid[[5]],
     er.oumv.hid = warts.petleng.er.hid[[6]],
     sym.bm1.hid = warts.petleng.sym.hid[[1]],
     sym.bms.hid = warts.petleng.sym.hid[[2]],
     sym.ou1.hid = warts.petleng.sym.hid[[3]],
     sym.oum.hid = warts.petleng.sym.hid[[4]],
     sym.ouma.hid = warts.petleng.sym.hid[[5]],
     sym.oumv.hid = warts.petleng.sym.hid[[6]],
     ard.bm1.hid = warts.petleng.ard.hid[[1]],
     ard.bms.hid = warts.petleng.ard.hid[[2]],
     ard.ou1.hid = warts.petleng.ard.hid[[3]],
     ard.oum.hid = warts.petleng.ard.hid[[4]],
     ard.ouma.hid = warts.petleng.ard.hid[[5]],
     #ard.oumv.hid = warts.petleng.ard.hid[[6]],
     er.bm1 = warts.petleng.er.nohid[[1]],
     er.bms = warts.petleng.er.nohid[[2]],
     er.ou1 = warts.petleng.er.nohid[[3]],
     er.oum = warts.petleng.er.nohid[[4]],
     er.ouma = warts.petleng.er.nohid[[5]],
     er.oumv = warts.petleng.er.nohid[[6]],
     sym.bm1 = warts.petleng.sym.nohid[[1]],
     sym.bms = warts.petleng.sym.nohid[[2]],
     sym.ou1 = warts.petleng.sym.nohid[[3]],
     sym.oum = warts.petleng.sym.nohid[[4]],
     sym.ouma = warts.petleng.sym.nohid[[5]],
     sym.oumv = warts.petleng.sym.nohid[[6]],
     ard.bm1 = warts.petleng.ard.nohid[[1]],
     ard.bms = warts.petleng.ard.nohid[[2]],
     ard.ou1 = warts.petleng.ard.nohid[[3]],
     ard.oum = warts.petleng.ard.nohid[[4]],
     ard.ouma = warts.petleng.ard.nohid[[5]],
     ard.oumv = warts.petleng.ard.nohid[[6]])

bic.index.wapl <- sapply(model.list.warts.petleng, function(x){tryCatch(x$BIC, error = function(y){NA})})
warts.petleng.table.hid <- getModelTable(model.list.warts.petleng[-which(abs(bic.index.wapl) > 1e04)])

if(any(abs(bic.index.wapl) > 1e04, na.rm = TRUE)){
    pars.warts.petleng <- getModelAvgParams(model.list.warts.petleng[-which(abs(bic.index.wapl) > 1e04)], force = FALSE)
} else {                                                                              pars.warts.petleng <- getModelAvgParams(model.list.warts.petleng, force = FALSE)}

plot.data.hid.warts.petleng <- reshape2::melt(pars.warts.petleng)
plot.data.hid.warts.petleng$discrete <- "Warts"
plot.data.hid.warts.petleng$continuous <- "Petiole Length"
plot.data.hid.warts.petleng$tip_state <- factor(c("Variable", "Differentiated", "Lost")[as.integer(plot.data.hid.warts.petleng$tip_state)], levels = c("Variable", "Differentiated", "Lost"))

## ggplot(plot.data.hid.warts.petleng, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")

### Stem Area

model.list.warts.stemarea <- list(er.bm1.hid = warts.stemarea.er.hid[[1]],
     er.bms.hid = warts.stemarea.er.hid[[2]],
     er.ou1.hid = warts.stemarea.er.hid[[3]],
     er.oum.hid = warts.stemarea.er.hid[[4]],
     er.ouma.hid = warts.stemarea.er.hid[[5]],
     er.oumv.hid = warts.stemarea.er.hid[[6]],
     sym.bm1.hid = warts.stemarea.sym.hid[[1]],
     sym.bms.hid = warts.stemarea.sym.hid[[2]],
     sym.ou1.hid = warts.stemarea.sym.hid[[3]],
     sym.oum.hid = warts.stemarea.sym.hid[[4]],
     sym.ouma.hid = warts.stemarea.sym.hid[[5]],
     #sym.oumv.hid = warts.stemarea.sym.hid[[6]],
     ard.bm1.hid = warts.stemarea.ard.hid[[1]],
     ard.bms.hid = warts.stemarea.ard.hid[[2]],
     ard.ou1.hid = warts.stemarea.ard.hid[[3]],
     ard.oum.hid = warts.stemarea.ard.hid[[4]],
     ard.ouma.hid = warts.stemarea.ard.hid[[5]],
     ard.oumv.hid = warts.stemarea.ard.hid[[6]],
     er.bm1 = warts.stemarea.er.nohid[[1]],
     er.bms = warts.stemarea.er.nohid[[2]],
     er.ou1 = warts.stemarea.er.nohid[[3]],
     er.oum = warts.stemarea.er.nohid[[4]],
     er.ouma = warts.stemarea.er.nohid[[5]],
     er.oumv = warts.stemarea.er.nohid[[6]],
     sym.bm1 = warts.stemarea.sym.nohid[[1]],
     sym.bms = warts.stemarea.sym.nohid[[2]],
     sym.ou1 = warts.stemarea.sym.nohid[[3]],
     sym.oum = warts.stemarea.sym.nohid[[4]],
     sym.ouma = warts.stemarea.sym.nohid[[5]],
     #sym.oumv = warts.stemarea.sym.nohid[[6]],
     ard.bm1 = warts.stemarea.ard.nohid[[1]],
     ard.bms = warts.stemarea.ard.nohid[[2]],
     ard.ou1 = warts.stemarea.ard.nohid[[3]],
     ard.oum = warts.stemarea.ard.nohid[[4]],
     ard.ouma = warts.stemarea.ard.nohid[[5]],
     ard.oumv = warts.stemarea.ard.nohid[[6]])

bic.index.wasa <- sapply(model.list.warts.stemarea, function(x){tryCatch(x$BIC, error = function(y){NA})})
warts.stemarea.table.hid <- getModelTable(model.list.warts.stemarea[-which(abs(bic.index.wasa) > 1e04)])

if(any(abs(bic.index.wasa) > 1e04, na.rm = TRUE)){
    pars.warts.stemarea <- getModelAvgParams(model.list.warts.stemarea[-which(abs(bic.index.wasa) > 1e04)], force = FALSE)
} else {                                                                              pars.warts.stemarea <- getModelAvgParams(model.list.warts.stemarea, force = FALSE)}

plot.data.hid.warts.stemarea <- reshape2::melt(pars.warts.stemarea)
plot.data.hid.warts.stemarea$discrete <- "Warts"
plot.data.hid.warts.stemarea$continuous <- "Stem Area"
plot.data.hid.warts.stemarea$tip_state <- factor(c("Variable", "Differentiated", "Lost")[as.integer(plot.data.hid.warts.stemarea$tip_state)], levels = c("Variable", "Differentiated", "Lost"))

## ggplot(plot.data.hid.warts.stemarea, aes(x = tip_state, y = value, colour = tip_state)) +
##     geom_point(size = 5, shape = 21) +
##     stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
##     stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
##     theme_classic() +
##     facet_wrap(~ variable, scales = "free")




theta.plot <- rbind(plot.data.hid.domgrow.corleng[grep("theta", plot.data.hid.domgrow.corleng$variable),],
                    plot.data.hid.domgrow.leafarea[grep("theta", plot.data.hid.domgrow.leafarea$variable),],
                    plot.data.hid.domgrow.petleng[grep("theta", plot.data.hid.domgrow.petleng$variable),],
                    plot.data.hid.domgrow.stemarea[grep("theta", plot.data.hid.domgrow.stemarea$variable),],
                    plot.data.hid.holediam.corleng[grep("theta", plot.data.hid.holediam.corleng$variable),],
                    plot.data.hid.holediam.leafarea[grep("theta", plot.data.hid.holediam.leafarea$variable),],
                    plot.data.hid.holediam.petleng[grep("theta", plot.data.hid.holediam.petleng$variable),],
                    plot.data.hid.holediam.stemarea[grep("theta", plot.data.hid.holediam.stemarea$variable),],
                    plot.data.hid.reward.corleng[grep("theta", plot.data.hid.reward.corleng$variable),],
                    plot.data.hid.reward.leafarea[grep("theta", plot.data.hid.reward.leafarea$variable),],
                    plot.data.hid.reward.petleng[grep("theta", plot.data.hid.reward.petleng$variable),],
                    plot.data.hid.reward.stemarea[grep("theta", plot.data.hid.reward.stemarea$variable),],
                    plot.data.hid.warts.corleng[grep("theta", plot.data.hid.warts.corleng$variable),],
                    plot.data.hid.warts.leafarea[grep("theta", plot.data.hid.warts.leafarea$variable),],
                    plot.data.hid.warts.petleng[grep("theta", plot.data.hid.warts.petleng$variable),],
                    plot.data.hid.warts.stemarea[grep("theta", plot.data.hid.warts.stemarea$variable),]
                    )
theta.plot$discrete <- factor(theta.plot$discrete, levels = c("Domatium Growth", "Reward", "Warts", "Hole Diameter"))

alpha.plot <- rbind(plot.data.hid.domgrow.corleng[grep("alpha", plot.data.hid.domgrow.corleng$variable),],
                    plot.data.hid.domgrow.leafarea[grep("alpha", plot.data.hid.domgrow.leafarea$variable),],
                    plot.data.hid.domgrow.petleng[grep("alpha", plot.data.hid.domgrow.petleng$variable),],
                    plot.data.hid.domgrow.stemarea[grep("alpha", plot.data.hid.domgrow.stemarea$variable),],
                    plot.data.hid.holediam.corleng[grep("alpha", plot.data.hid.holediam.corleng$variable),],
                    plot.data.hid.holediam.leafarea[grep("alpha", plot.data.hid.holediam.leafarea$variable),],
                    plot.data.hid.holediam.petleng[grep("alpha", plot.data.hid.holediam.petleng$variable),],
                    plot.data.hid.holediam.stemarea[grep("alpha", plot.data.hid.holediam.stemarea$variable),],
                    plot.data.hid.reward.corleng[grep("alpha", plot.data.hid.reward.corleng$variable),],
                    plot.data.hid.reward.leafarea[grep("alpha", plot.data.hid.reward.leafarea$variable),],
                    plot.data.hid.reward.petleng[grep("alpha", plot.data.hid.reward.petleng$variable),],
                    plot.data.hid.reward.stemarea[grep("alpha", plot.data.hid.reward.stemarea$variable),],
                    plot.data.hid.warts.corleng[grep("alpha", plot.data.hid.warts.corleng$variable),],
                    plot.data.hid.warts.leafarea[grep("alpha", plot.data.hid.warts.leafarea$variable),],
                    plot.data.hid.warts.petleng[grep("alpha", plot.data.hid.warts.petleng$variable),],
                    plot.data.hid.warts.stemarea[grep("alpha", plot.data.hid.warts.stemarea$variable),]
                    )
alpha.plot$discrete <- factor(alpha.plot$discrete, levels = c("Domatium Growth", "Reward", "Warts", "Hole Diameter"))


sigma.plot <- rbind(plot.data.hid.domgrow.corleng[grep("sigma", plot.data.hid.domgrow.corleng$variable),],
                    plot.data.hid.domgrow.leafarea[grep("sigma", plot.data.hid.domgrow.leafarea$variable),],
                    plot.data.hid.domgrow.petleng[grep("sigma", plot.data.hid.domgrow.petleng$variable),],
                    plot.data.hid.domgrow.stemarea[grep("sigma", plot.data.hid.domgrow.stemarea$variable),],
                    plot.data.hid.holediam.corleng[grep("sigma", plot.data.hid.holediam.corleng$variable),],
                    plot.data.hid.holediam.leafarea[grep("sigma", plot.data.hid.holediam.leafarea$variable),],
                    plot.data.hid.holediam.petleng[grep("sigma", plot.data.hid.holediam.petleng$variable),],
                    plot.data.hid.holediam.stemarea[grep("sigma", plot.data.hid.holediam.stemarea$variable),],
                    plot.data.hid.reward.corleng[grep("sigma", plot.data.hid.reward.corleng$variable),],
                    plot.data.hid.reward.leafarea[grep("sigma", plot.data.hid.reward.leafarea$variable),],
                    plot.data.hid.reward.petleng[grep("sigma", plot.data.hid.reward.petleng$variable),],
                    plot.data.hid.reward.stemarea[grep("sigma", plot.data.hid.reward.stemarea$variable),],
                    plot.data.hid.warts.corleng[grep("sigma", plot.data.hid.warts.corleng$variable),],
                    plot.data.hid.warts.leafarea[grep("sigma", plot.data.hid.warts.leafarea$variable),],
                    plot.data.hid.warts.petleng[grep("sigma", plot.data.hid.warts.petleng$variable),],
                    plot.data.hid.warts.stemarea[grep("sigma", plot.data.hid.warts.stemarea$variable),]
                    )
sigma.plot$discrete <- factor(sigma.plot$discrete, levels = c("Domatium Growth", "Reward", "Warts", "Hole Diameter"))


ggplot(theta.plot, aes(x = tip_state, y = value, colour = tip_state, alpha = tip_state)) +
    geom_point(size = 5, shape = 21) +
    stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    scale_colour_manual(values = c(met.brewer("Greek", n = 2, type = "continuous"),
                                   met.brewer("Veronese", n = 3, type = "continuous"),
                                   met.brewer("Hiroshige", n = 2, type = "continuous"),
                                   met.brewer("Isfahan2", n = 3, type = "continuous"))) +
    scale_alpha_manual(values = rep(0.5, 10)) +
    labs(x = "Tip State", y = "Parameter Value (\u03b8)") +
    ggh4x::facet_grid2(c("discrete", "continuous"), scales = "free", independent = "all") +
    cowplot::theme_cowplot() +
    theme(legend.position = "none")
ggsave(here("output/houwie/hidden_state/figs/theta_hidden.png"), bg = "white", dpi = 300, height = 10, width = 18, units = "in")
ggsave(here("output/houwie/hidden_state/figs/theta_hidden.pdf"), bg = "white", height = 10, width = 18, units = "in")

ggplot(alpha.plot, aes(x = tip_state, y = value, colour = tip_state, alpha = tip_state)) +
    geom_point(size = 5, shape = 21) +
    stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    scale_colour_manual(values = c(met.brewer("Greek", n = 2, type = "continuous"),
                                   met.brewer("Veronese", n = 3, type = "continuous"),
                                   met.brewer("Hiroshige", n = 2, type = "continuous"),
                                   met.brewer("Isfahan2", n = 3, type = "continuous"))) +
    scale_alpha_manual(values = rep(0.5, 10)) +
    labs(x = "Tip State", y = "Parameter Value (\u03b1)") +
    ggh4x::facet_grid2(c("discrete", "continuous"), scales = "free", independent = "all") +
    cowplot::theme_cowplot() +
    theme(legend.position = "none")
ggsave(here("manuscript/figs/houwie/alpha_houwie.png"), bg = "white", dpi = 300, height = 10, width = 18, units = "in")
ggsave(here("manuscript/figs/houwie/alpha_houwie.pdf"), bg = "white", height = 10, width = 18, units = "in")


ggplot(sigma.plot, aes(x = tip_state, y = value, colour = tip_state, alpha = tip_state)) +
    geom_point(size = 5, shape = 21) +
    stat_summary(fun = mean, geom = "point", aes(group = 1, size = 2)) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    scale_colour_manual(values = c(met.brewer("Greek", n = 2, type = "continuous"),
                                   met.brewer("Veronese", n = 3, type = "continuous"),
                                   met.brewer("Hiroshige", n = 2, type = "continuous"),
                                   met.brewer("Isfahan2", n = 3, type = "continuous"))) +
    scale_alpha_manual(values = rep(0.5, 10)) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    ggh4x::facet_grid2(c("discrete", "continuous"), scales = "free", independent = "all") +
    cowplot::theme_cowplot() +
    theme(legend.position = "none")

ggsave(here("manuscript/figs/houwie/sigma_houwie.png"), bg = "white", dpi = 300, height = 10, width = 18, units = "in")
ggsave(here("manuscript/figs/houwie/sigma_houwie.pdf"), bg = "white", height = 10, width = 18, units = "in")

## Using dentist package to estimate confidence intervals

domgrow.petleng.er.OU1.dent <- hOUwie.walk(domgrow.petleng.er[[3]], delta = 2, nsteps = 200)
domgrow.petleng.er.OUM.dent <- hOUwie.walk(domgrow.petleng.er[[4]], delta = 2, nsteps = 200)
domgrow.petleng.er.OUMA.dent <- hOUwie.walk(domgrow.petleng.er[[5]], delta = 2, nsteps = 200)
domgrow.petleng.er.OUMV.dent <- hOUwie.walk(domgrow.petleng.er[[6]], delta = 1, nsteps = 200)

domgrow.petleng.sym.OU1.dent <- hOUwie.walk(domgrow.petleng.sym[[3]], delta = 2, nsteps = 200)
domgrow.petleng.sym.OUM.dent <- hOUwie.walk(domgrow.petleng.sym[[4]], delta = 2, nsteps = 200)
domgrow.petleng.sym.OUMA.dent <- hOUwie.walk(domgrow.petleng.sym[[5]], delta = 2, nsteps = 200)
domgrow.petleng.sym.OUMV.dent <- hOUwie.walk(domgrow.petleng.sym[[6]], delta = 2, nsteps = 200, debug = TRUE)

domgrow.petleng.ard.OU1.dent <- hOUwie.walk(domgrow.petleng.ard[[3]], delta = 2, nsteps = 200)
domgrow.petleng.ard.OUM.dent <- hOUwie.walk(domgrow.petleng.ard[[4]], delta = 2, nsteps = 200)
domgrow.petleng.ard.OUMA.dent <- hOUwie.walk(domgrow.petleng.ard[[5]], delta = 2, nsteps = 200)
domgrow.petleng.ard.OUMV.dent <- hOUwie.walk(domgrow.petleng.ard[[6]], delta = 1, nsteps = 200)

domgrow.petleng.er.table.hid <- getModelTable(model.list.domgrow.petleng)

## Theta 1
ci.lower.theta1 <-
    (domgrow.petleng.er.OU1.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[3]) +
    (domgrow.petleng.er.OUM.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[4]) +
    (domgrow.petleng.er.OUMA.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[5]) +
    (domgrow.petleng.er.OUMV.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[6]) +
    (domgrow.petleng.sym.OU1.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[9]) +
    (domgrow.petleng.sym.OUM.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[10]) +
    (domgrow.petleng.sym.OUMA.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[11]) +
    (domgrow.petleng.sym.OUMV.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[12]) +
    (domgrow.petleng.ard.OU1.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[15]) +
    (domgrow.petleng.ard.OUM.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[16]) +
    (domgrow.petleng.ard.OUMA.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[17]) +
    (domgrow.petleng.ard.OUMV.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[18])


ci.upper.theta1 <-
    (domgrow.petleng.er.OU1.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[3]) +
    (domgrow.petleng.er.OUM.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[4]) +
    (domgrow.petleng.er.OUMA.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[5]) +
    (domgrow.petleng.er.OUMV.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[6]) +
    (domgrow.petleng.sym.OU1.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[9]) +
    (domgrow.petleng.sym.OUM.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[10]) +
    (domgrow.petleng.sym.OUMA.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[11]) +
    (domgrow.petleng.sym.OUMV.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[12]) +
    (domgrow.petleng.ard.OU1.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[15]) +
    (domgrow.petleng.ard.OUM.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[16]) +
    (domgrow.petleng.ard.OUMA.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[17]) +
    (domgrow.petleng.ard.OUMV.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[18])

## Theta 2
ci.lower.theta2 <-
    (domgrow.petleng.er.OU1.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[3]) +
    (domgrow.petleng.er.OUM.dent$all_ranges$theta_2[2] * domgrow.petleng.er.table.hid$BICwt[4]) +
    (domgrow.petleng.er.OUMA.dent$all_ranges$theta_2[2] * domgrow.petleng.er.table.hid$BICwt[5]) +
    (domgrow.petleng.er.OUMV.dent$all_ranges$theta_2[2] * domgrow.petleng.er.table.hid$BICwt[6]) +
    (domgrow.petleng.sym.OU1.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[9]) +
    (domgrow.petleng.sym.OUM.dent$all_ranges$theta_2[2] * domgrow.petleng.er.table.hid$BICwt[10]) +
    (domgrow.petleng.sym.OUMA.dent$all_ranges$theta_2[2] * domgrow.petleng.er.table.hid$BICwt[11]) +
    (domgrow.petleng.sym.OUMV.dent$all_ranges$theta_2[2] * domgrow.petleng.er.table.hid$BICwt[12]) +
    (domgrow.petleng.ard.OU1.dent$all_ranges$theta_1[2] * domgrow.petleng.er.table.hid$BICwt[15]) +
    (domgrow.petleng.ard.OUM.dent$all_ranges$theta_2[2] * domgrow.petleng.er.table.hid$BICwt[16]) +
    (domgrow.petleng.ard.OUMA.dent$all_ranges$theta_2[2] * domgrow.petleng.er.table.hid$BICwt[17]) +
    (domgrow.petleng.ard.OUMV.dent$all_ranges$theta_2[2] * domgrow.petleng.er.table.hid$BICwt[18])

ci.upper.theta2 <-
    (domgrow.petleng.er.OU1.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[3]) +
    (domgrow.petleng.er.OUM.dent$all_ranges$theta_2[3] * domgrow.petleng.er.table.hid$BICwt[4]) +
    (domgrow.petleng.er.OUMA.dent$all_ranges$theta_2[3] * domgrow.petleng.er.table.hid$BICwt[5]) +
    (domgrow.petleng.er.OUMV.dent$all_ranges$theta_2[3] * domgrow.petleng.er.table.hid$BICwt[6]) +
    (domgrow.petleng.sym.OU1.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[9]) +
    (domgrow.petleng.sym.OUM.dent$all_ranges$theta_2[3] * domgrow.petleng.er.table.hid$BICwt[10]) +
    (domgrow.petleng.sym.OUMA.dent$all_ranges$theta_2[3] * domgrow.petleng.er.table.hid$BICwt[11]) +
    (domgrow.petleng.sym.OUMV.dent$all_ranges$theta_2[3] * domgrow.petleng.er.table.hid$BICwt[12]) +
    (domgrow.petleng.ard.OU1.dent$all_ranges$theta_1[3] * domgrow.petleng.er.table.hid$BICwt[15]) +
    (domgrow.petleng.ard.OUM.dent$all_ranges$theta_2[3] * domgrow.petleng.er.table.hid$BICwt[16]) +
    (domgrow.petleng.ard.OUMA.dent$all_ranges$theta_2[3] * domgrow.petleng.er.table.hid$BICwt[17]) +
    (domgrow.petleng.ard.OUMV.dent$all_ranges$theta_2[3] * domgrow.petleng.er.table.hid$BICwt[18])

## Alpha 1
ci.lower.alpha1 <-
    (domgrow.petleng.er.OU1.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[3]) +
    (domgrow.petleng.er.OUM.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[4]) +
    (domgrow.petleng.er.OUMA.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[5]) +
    (domgrow.petleng.er.OUMV.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[6]) +
    (domgrow.petleng.sym.OU1.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[9]) +
    (domgrow.petleng.sym.OUM.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[10]) +
    (domgrow.petleng.sym.OUMA.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[11]) +
    (domgrow.petleng.sym.OUMV.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[12]) +
    (domgrow.petleng.ard.OU1.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[15]) +
    (domgrow.petleng.ard.OUM.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[16]) +
    (domgrow.petleng.ard.OUMA.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[17]) +
    (domgrow.petleng.ard.OUMV.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[18])

ci.upper.alpha1 <-
    (domgrow.petleng.er.OU1.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[3]) +
    (domgrow.petleng.er.OUM.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[4]) +
    (domgrow.petleng.er.OUMA.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[5]) +
    (domgrow.petleng.er.OUMV.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[6]) +
    (domgrow.petleng.sym.OU1.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[9]) +
    (domgrow.petleng.sym.OUM.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[10]) +
    (domgrow.petleng.sym.OUMA.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[11]) +
    (domgrow.petleng.sym.OUMV.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[12]) +
    (domgrow.petleng.ard.OU1.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[15]) +
    (domgrow.petleng.ard.OUM.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[16]) +
    (domgrow.petleng.ard.OUMA.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[17]) +
    (domgrow.petleng.ard.OUMV.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[18])

## Alpha 2
ci.lower.alpha2 <-
    (domgrow.petleng.er.OU1.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[3]) +
    (domgrow.petleng.er.OUM.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[4]) +
    (domgrow.petleng.er.OUMA.dent$all_ranges$alpha_2[2] * domgrow.petleng.er.table.hid$BICwt[5]) +
    (domgrow.petleng.er.OUMV.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[6]) +
    (domgrow.petleng.sym.OU1.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[9]) +
    (domgrow.petleng.sym.OUM.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[10]) +
    (domgrow.petleng.sym.OUMA.dent$all_ranges$alpha_2[2] * domgrow.petleng.er.table.hid$BICwt[11]) +
    (domgrow.petleng.sym.OUMV.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[12]) +
    (domgrow.petleng.ard.OU1.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[15]) +
    (domgrow.petleng.ard.OUM.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[16]) +
    (domgrow.petleng.ard.OUMA.dent$all_ranges$alpha_2[2] * domgrow.petleng.er.table.hid$BICwt[17]) +
    (domgrow.petleng.ard.OUMV.dent$all_ranges$alpha_1[2] * domgrow.petleng.er.table.hid$BICwt[18])

ci.upper.alpha2 <-
    (domgrow.petleng.er.OU1.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[3]) +
    (domgrow.petleng.er.OUM.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[4]) +
    (domgrow.petleng.er.OUMA.dent$all_ranges$alpha_2[3] * domgrow.petleng.er.table.hid$BICwt[5]) +
    (domgrow.petleng.er.OUMV.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[6]) +
    (domgrow.petleng.sym.OU1.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[9]) +
    (domgrow.petleng.sym.OUM.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[10]) +
    (domgrow.petleng.sym.OUMA.dent$all_ranges$alpha_2[3] * domgrow.petleng.er.table.hid$BICwt[11]) +
    (domgrow.petleng.sym.OUMV.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[12]) +
    (domgrow.petleng.ard.OU1.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[15]) +
    (domgrow.petleng.ard.OUM.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[16]) +
    (domgrow.petleng.ard.OUMA.dent$all_ranges$alpha_2[3] * domgrow.petleng.er.table.hid$BICwt[17]) +
    (domgrow.petleng.ard.OUMV.dent$all_ranges$alpha_1[3] * domgrow.petleng.er.table.hid$BICwt[18])


## Sigma2 1
ci.lower.sigma21 <-
    (domgrow.petleng.er.OU1.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[3]) +
    (domgrow.petleng.er.OUM.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[4]) +
    (domgrow.petleng.er.OUMA.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[5]) +
    (domgrow.petleng.er.OUMV.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[6]) +
    (domgrow.petleng.sym.OU1.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[9]) +
    (domgrow.petleng.sym.OUM.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[10]) +
    (domgrow.petleng.sym.OUMA.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[11]) +
    (domgrow.petleng.sym.OUMV.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[12]) +
    (domgrow.petleng.ard.OU1.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[15]) +
    (domgrow.petleng.ard.OUM.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[16]) +
    (domgrow.petleng.ard.OUMA.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[17]) +
    (domgrow.petleng.ard.OUMV.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[18])

ci.upper.sigma21 <-
    (domgrow.petleng.er.OU1.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[3]) +
    (domgrow.petleng.er.OUM.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[4]) +
    (domgrow.petleng.er.OUMA.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[5]) +
    (domgrow.petleng.er.OUMV.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[6]) +
    (domgrow.petleng.sym.OU1.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[9]) +
    (domgrow.petleng.sym.OUM.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[10]) +
    (domgrow.petleng.sym.OUMA.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[11]) +
    (domgrow.petleng.sym.OUMV.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[12]) +
    (domgrow.petleng.ard.OU1.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[15]) +
    (domgrow.petleng.ard.OUM.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[16]) +
    (domgrow.petleng.ard.OUMA.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[17]) +
    (domgrow.petleng.ard.OUMV.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[18])

## Sigma2 2
ci.lower.sigma22 <-
    (domgrow.petleng.er.OU1.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[3]) +
    (domgrow.petleng.er.OUM.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[4]) +
    (domgrow.petleng.er.OUMA.dent$all_ranges$sigma2_2[2] * domgrow.petleng.er.table.hid$BICwt[5]) +
    (domgrow.petleng.er.OUMV.dent$all_ranges$sigma2_2[2] * domgrow.petleng.er.table.hid$BICwt[6]) +
    (domgrow.petleng.sym.OU1.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[9]) +
    (domgrow.petleng.sym.OUM.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[10]) +
    (domgrow.petleng.sym.OUMA.dent$all_ranges$sigma2_2[2] * domgrow.petleng.er.table.hid$BICwt[11]) +
    (domgrow.petleng.sym.OUMV.dent$all_ranges$sigma2_2[2] * domgrow.petleng.er.table.hid$BICwt[12]) +
    (domgrow.petleng.ard.OU1.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[15]) +
    (domgrow.petleng.ard.OUM.dent$all_ranges$sigma2_1[2] * domgrow.petleng.er.table.hid$BICwt[16]) +
    (domgrow.petleng.ard.OUMA.dent$all_ranges$sigma2_2[2] * domgrow.petleng.er.table.hid$BICwt[17]) +
    (domgrow.petleng.ard.OUMV.dent$all_ranges$sigma2_2[2] * domgrow.petleng.er.table.hid$BICwt[18])

ci.upper.sigma22 <-
    (domgrow.petleng.er.OU1.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[3]) +
    (domgrow.petleng.er.OUM.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[4]) +
    (domgrow.petleng.er.OUMA.dent$all_ranges$sigma2_2[3] * domgrow.petleng.er.table.hid$BICwt[5]) +
    (domgrow.petleng.er.OUMV.dent$all_ranges$sigma2_2[3] * domgrow.petleng.er.table.hid$BICwt[6]) +
    (domgrow.petleng.sym.OU1.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[9]) +
    (domgrow.petleng.sym.OUM.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[10]) +
    (domgrow.petleng.sym.OUMA.dent$all_ranges$sigma2_2[3] * domgrow.petleng.er.table.hid$BICwt[11]) +
    (domgrow.petleng.sym.OUMV.dent$all_ranges$sigma2_2[3] * domgrow.petleng.er.table.hid$BICwt[12]) +
    (domgrow.petleng.ard.OU1.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[15]) +
    (domgrow.petleng.ard.OUM.dent$all_ranges$sigma2_1[3] * domgrow.petleng.er.table.hid$BICwt[16]) +
    (domgrow.petleng.ard.OUMA.dent$all_ranges$sigma2_2[3] * domgrow.petleng.er.table.hid$BICwt[17]) +
    (domgrow.petleng.ard.OUMV.dent$all_ranges$sigma2_2[3] * domgrow.petleng.er.table.hid$BICwt[18])

