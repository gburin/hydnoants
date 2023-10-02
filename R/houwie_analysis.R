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
library("phylolm")

i_am("R/houwie_analysis.R")

fulldata <- read.csv(here("data/Dataset_traits_Hydnophytinae_CODED_2022.csv"))

names(fulldata)[21:ncol(fulldata)] <- c("strategy", "warts", "stem.area", "leaf.area", "hole.diameter", "architecture", "dom.growth", "reward", "corola.length", "mating.system", "petiole.length", "leaf.struct", "appendages", "holediam.disc")

fulldata$species[grep("Wor2g", fulldata$species)] <- gsub("Wor2g", "Worthing", fulldata$species[grep("Wor2g", fulldata$species)])
fulldata$species[grep("areolata", fulldata$species)] <- gsub("areolata", "aerolata", fulldata$species[grep("areolata", fulldata$species)])
fulldata$hole.diameter[which(fulldata$hole.diameter == "NA (outgroup has no domatia)")] <- NA
fulldata$hole.diameter <- as.numeric(fulldata$hole.diameter)


mcc.tree <- ape::read.tree(here("data/MCCtree135taxa.tre"))

clim.pca <- prcomp(fulldata[, 2:20])
summary(clim.pca)

fulldata <- cbind(fulldata, clim.pca$x[, 1:3])
names(fulldata)[(ncol(fulldata) - 2):ncol(fulldata)] <- paste0("climpc", 1:3)

data.pgls <- comparative.data(mcc.tree, fulldata, names.col = "species", na.omit = FALSE)


## ## Building models to check for significant relationships

## ### Leaf Area
## summary(lm.leafarea.pc1 <- lm(data.pgls$data$leaf.area ~ data.pgls$data$climpc1))
## summary(lm.leafarea.pc2 <- lm(data.pgls$data$leaf.area ~ data.pgls$data$climpc2))
## summary(lm.leafarea.pc3 <- lm(data.pgls$data$leaf.area ~ data.pgls$data$climpc3))

## summary(pgls.leafarea.pc1 <- pgls(leaf.area ~ climpc1, data.pgls, lambda = "ML"))
## summary(pgls.leafarea.pc2 <- pgls(leaf.area ~ climpc2, data.pgls, lambda = "ML"))
## summary(pgls.leafarea.pc3 <- pgls(leaf.area ~ climpc3, data.pgls, lambda = "ML"))

## ### Stem Area
## summary(lm.stemarea.pc1 <- lm(data.pgls$data$stem.area ~ data.pgls$data$climpc1))
## summary(lm.stemarea.pc2 <- lm(data.pgls$data$stem.area ~ data.pgls$data$climpc2))
## summary(lm.stemarea.pc3 <- lm(data.pgls$data$stem.area ~ data.pgls$data$climpc3)) ## p = 0.04825

## summary(pgls.stemarea.pc1 <- pgls(stem.area ~ climpc1, data.pgls, lambda = "ML"))
## summary(pgls.stemarea.pc2 <- pgls(stem.area ~ climpc2, data.pgls, lambda = "ML"))
## summary(pgls.stemarea.pc3 <- pgls(stem.area ~ climpc3, data.pgls, lambda = "ML"))

## ### Hole Diameter
## summary(lm.hole.diameter.pc1 <- lm(data.pgls$data$hole.diameter ~ data.pgls$data$climpc1)) ## p = 0.01772
## summary(lm.hole.diameter.pc2 <- lm(data.pgls$data$hole.diameter ~ data.pgls$data$climpc2)) ## p = 0.03268
## summary(lm.hole.diameter.pc3 <- lm(data.pgls$data$hole.diameter ~ data.pgls$data$climpc3))

## summary(pgls.hole.diameter.pc1 <- pgls(hole.diameter ~ climpc1, data.pgls, lambda = "ML"))
## summary(pgls.hole.diameter.pc2 <- pgls(hole.diameter ~ climpc2, data.pgls, lambda = "ML")) ## p = 0.1422
## summary(pgls.hole.diameter.pc3 <- pgls(hole.diameter ~ climpc3, data.pgls, lambda = "ML"))

## ### Corola Length
## summary(lm.corola.length.pc1 <- lm(data.pgls$data$corola.length ~ data.pgls$data$climpc1))
## summary(lm.corola.length.pc2 <- lm(data.pgls$data$corola.length ~ data.pgls$data$climpc2))
## summary(lm.corola.length.pc3 <- lm(data.pgls$data$corola.length ~ data.pgls$data$climpc3)) ## p = 0.4697

## summary(pgls.corola.length.pc1 <- pgls(corola.length ~ climpc1, data.pgls, lambda = "ML"))
## summary(pgls.corola.length.pc2 <- pgls(corola.length ~ climpc2, data.pgls, lambda = "ML"))
## summary(pgls.corola.length.pc3 <- pgls(corola.length ~ climpc3, data.pgls, lambda = "ML"))

## ### Petiole Length
## summary(lm.petiole.length.pc1 <- lm(data.pgls$data$petiole.length ~ data.pgls$data$climpc1))
## summary(lm.petiole.length.pc2 <- lm(data.pgls$data$petiole.length ~ data.pgls$data$climpc2))
## summary(lm.petiole.length.pc3 <- lm(data.pgls$data$petiole.length ~ data.pgls$data$climpc3))

## summary(pgls.petiole.length.pc1 <- pgls(petiole.length ~ climpc1, data.pgls, lambda = "ML"))
## summary(pgls.petiole.length.pc2 <- pgls(petiole.length ~ climpc2, data.pgls, lambda = "ML"))
## summary(pgls.petiole.length.pc3 <- pgls(petiole.length ~ climpc3, data.pgls, lambda = "ML"))



## fulldata$stem.area.pc3.lm <- summary(lm.stemarea.pc3)$residuals

## fulldata$hole.diameter.pc1.lm[which(!is.na(fulldata$hole.diameter))] <- summary(lm.hole.diameter.pc1)$residuals
## fulldata$hole.diameter.pc2.lm[which(!is.na(fulldata$hole.diameter))] <- summary(lm.hole.diameter.pc2)$residuals
## fulldata$hole.diameter.pc2.pgls[which(!is.na(fulldata$hole.diameter))] <- summary(pgls.hole.diameter.pc2)$residuals[,1]

## fulldata$corola.length.pc3.lm[which(!is.na(fulldata$corola.length))] <- summary(lm.corola.length.pc3)$residuals

## saveRDS(fulldata, file = here("output/fulldata_houwie.RDS"))

fulldata <- readRDS(here("output/fulldata_houwie.RDS"))
fulldata.ingroup <- fulldata[fulldata$strategy != 0,]

## ## Keeping data as originally coded for SIMMAP

## traits.simmap <- fulldata

## ## Recoding trait states for plotting

## app <- c("outgroup", "none", "variable", "spines")
## arch <- c("shrub", "multiple", "single")
## domgrow <- c("outgroup", "diffuse", "apical")
## leafstruc <- c("thick", "variable", "thin", "succulent")
## matsys <- c("heterostylous", "non-heterostylous", "funct_unisexual")
## reward <- c("outgroup", "absent", "present")
## strategy <- c("outgroup", "facultative", "obligate", "lost")
## warts <- c("absent", "variable", "differentiated", "lost")
## holediam.disc <- c("absent", "severallargebase_smallside", "onelargebase_smalldise", "alllarge")

## fulldata$strategy <- strategy[fulldata$strategy + 1]
## fulldata$warts <- warts[fulldata$warts + 1]
## fulldata$reward <- reward[fulldata$reward + 1]
## fulldata$architecture <- arch[fulldata$architecture + 1]
## fulldata$dom.growth <- domgrow[fulldata$dom.growth + 1]
## fulldata$mating.system <- matsys[fulldata$mating.system + 1]
## fulldata$leaf.struct <- leafstruc[fulldata$leaf.struct + 1]
## fulldata$appendages <- app[fulldata$appendages + 1]
## fulldata$holediam.disc <- holediam.disc[fulldata$holediam.disc + 1]

## fulldata.names <- fulldata.ingroup[,1]
## fulldata.ingroup <- fulldata.ingroup[,-1]
## rownames(fulldata.ingroup) <- fulldata.names
fulldata.ingroup <- fulldata.ingroup[, c(1, 21:34)]
disc.traits <- sort(names(fulldata.ingroup)[c(2, 3, 7, 8, 9, 11, 13, 14, 15)])
cont.traits <- sort(names(fulldata.ingroup)[c(4, 5, 10, 12)])

cont.models <- c("BM1", "BMV", "OU1", "OUM", "OUMV", "OUMA")
disc.models <- c("ER", "SYM", "ARD")

registerDoMC(56)

### No hidden states

###############################
#### Appendages + Corola Length
##### ER
## app.corleng.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     app.corleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "appendages", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(app.corleng.er.nohid, file = here("output/houwie/app-corleng-er-nohid.RDS"))
## ##### SYM
## app.corleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     app.corleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "appendages", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(app.corleng.sym.nohid, file = here("output/houwie/app-corleng-sym-nohid.RDS"))
## ##### ARD
## app.corleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     app.corleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "appendages", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(app.corleng.ard.nohid, file = here("output/houwie/app-corleng-ard-nohid.RDS"))

## #### Appendages + Leaf Area
## ##### ER
## app.leafarea.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     app.leafarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "appendages", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(app.leafarea.er.nohid, file = here("output/houwie/app-leafarea-er-nohid.RDS"))
## ##### SYM
## app.leafarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     app.leafarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "appendages", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(app.leafarea.sym.nohid, file = here("output/houwie/app-leafarea-sym-nohid.RDS"))
## ##### ARD
## app.leafarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     app.leafarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "appendages", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(app.leafarea.ard.nohid, file = here("output/houwie/app-leafarea-ard-nohid.RDS"))


## #### Appendages + Petiole Length
## ##### ER
## app.petleng.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     app.petleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "appendages", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(app.petleng.er.nohid, file = here("output/houwie/app-petleng-er-nohid.RDS"))
## ##### SYM
## app.petleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     app.petleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "appendages", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(app.petleng.sym.nohid, file = here("output/houwie/app-petleng-sym-nohid.RDS"))
## ##### ARD
## app.petleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     app.petleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "appendages", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(app.petleng.ard.nohid, file = here("output/houwie/app-petleng-ard-nohid.RDS"))


## #### Appendages + Stem Area
## ##### ER
## app.stemarea.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     app.stemarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "appendages", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(app.stemarea.er.nohid, file = here("output/houwie/app-stemarea-er-nohid.RDS"))
## ##### SYM
## app.stemarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     app.stemarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "appendages", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(app.stemarea.sym.nohid, file = here("output/houwie/app-stemarea-sym-nohid.RDS"))
## ##### ARD
## app.stemarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     app.stemarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "appendages", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(app.stemarea.ard.nohid, file = here("output/houwie/app-stemarea-ard-nohid.RDS"))







###############################
#### Architecture + Corola Length
##### ER
## arch.corleng.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     arch.corleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "architecture", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(arch.corleng.er.nohid, file = here("output/houwie/arch-corleng-er-nohid.RDS"))
## ##### SYM
## arch.corleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     arch.corleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "architecture", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(arch.corleng.sym.nohid, file = here("output/houwie/arch-corleng-sym-nohid.RDS"))
## ##### ARD
## arch.corleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     arch.corleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "architecture", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(arch.corleng.ard.nohid, file = here("output/houwie/arch-corleng-ard-nohid.RDS"))

## #### Architecture + Leaf Area
## ##### ER
## arch.leafarea.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     arch.leafarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "architecture", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(arch.leafarea.er.nohid, file = here("output/houwie/arch-leafarea-er-nohid.RDS"))
## ##### SYM
## arch.leafarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     arch.leafarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "architecture", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(arch.leafarea.sym.nohid, file = here("output/houwie/arch-leafarea-sym-nohid.RDS"))
## ##### ARD
## arch.leafarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     arch.leafarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "architecture", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(arch.leafarea.ard.nohid, file = here("output/houwie/arch-leafarea-ard-nohid.RDS"))


## #### Architecture + Petiole Length
## ##### ER
## arch.petleng.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     arch.petleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "architecture", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(arch.petleng.er.nohid, file = here("output/houwie/arch-petleng-er-nohid.RDS"))
## ##### SYM
## arch.petleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     arch.petleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "architecture", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(arch.petleng.sym.nohid, file = here("output/houwie/arch-petleng-sym-nohid.RDS"))
## ##### ARD
## arch.petleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     arch.petleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "architecture", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(arch.petleng.ard.nohid, file = here("output/houwie/arch-petleng-ard-nohid.RDS"))


## #### Architecture + Stem Area
## ##### ER
## arch.stemarea.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     arch.stemarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "architecture", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(arch.stemarea.er.nohid, file = here("output/houwie/arch-stemarea-er-nohid.RDS"))
## ##### SYM
## arch.stemarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     arch.stemarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "architecture", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(arch.stemarea.sym.nohid, file = here("output/houwie/arch-stemarea-sym-nohid.RDS"))
## ##### ARD
## arch.stemarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     arch.stemarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "architecture", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(arch.stemarea.ard.nohid, file = here("output/houwie/arch-stemarea-ard-nohid.RDS"))







###############################
#### Domatium Growth + Corola Length
##### ER
domgrow.corleng.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.corleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.corleng.er.nohid, file = here("output/houwie/domgrow-corleng-er-nohid.RDS"))
##### SYM
domgrow.corleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.corleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.corleng.sym.nohid, file = here("output/houwie/domgrow-corleng-sym-nohid.RDS"))
##### ARD
domgrow.corleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.corleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.corleng.ard.nohid, file = here("output/houwie/domgrow-corleng-ard-nohid.RDS"))

#### Domatium Growth + Leaf Area
##### ER
domgrow.leafarea.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.leafarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.leafarea.er.nohid, file = here("output/houwie/domgrow-leafarea-er-nohid.RDS"))
##### SYM
domgrow.leafarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.leafarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.leafarea.sym.nohid, file = here("output/houwie/domgrow-leafarea-sym-nohid.RDS"))
##### ARD
domgrow.leafarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.leafarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.leafarea.ard.nohid, file = here("output/houwie/domgrow-leafarea-ard-nohid.RDS"))


#### Domatium Growth + Petiole Length
##### ER
domgrow.petleng.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.petleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.petleng.er.nohid, file = here("output/houwie/domgrow-petleng-er-nohid.RDS"))
##### SYM
domgrow.petleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.petleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.petleng.sym.nohid, file = here("output/houwie/domgrow-petleng-sym-nohid.RDS"))
##### ARD
domgrow.petleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.petleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.petleng.ard.nohid, file = here("output/houwie/domgrow-petleng-ard-nohid.RDS"))


#### Domatium Growth + Stem Area
##### ER
domgrow.stemarea.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.stemarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.stemarea.er.nohid, file = here("output/houwie/domgrow-stemarea-er-nohid.RDS"))
##### SYM
domgrow.stemarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.stemarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.stemarea.sym.nohid, file = here("output/houwie/domgrow-stemarea-sym-nohid.RDS"))
##### ARD
domgrow.stemarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.stemarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.stemarea.ard.nohid, file = here("output/houwie/domgrow-stemarea-ard-nohid.RDS"))







###############################
#### Hole Diameter Discrete + Corola Length
##### ER
holediam.corleng.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.corleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.corleng.er.nohid, file = here("output/houwie/holediam-corleng-er-nohid.RDS"))
##### SYM
holediam.corleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.corleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.corleng.sym.nohid, file = here("output/houwie/holediam-corleng-sym-nohid.RDS"))
##### ARD
holediam.corleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.corleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.corleng.ard.nohid, file = here("output/houwie/holediam-corleng-ard-nohid.RDS"))

#### Hole Diameter Discrete + Leaf Area
##### ER
holediam.leafarea.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.leafarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.leafarea.er.nohid, file = here("output/houwie/holediam-leafarea-er-nohid.RDS"))
##### SYM
holediam.leafarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.leafarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.leafarea.sym.nohid, file = here("output/houwie/holediam-leafarea-sym-nohid.RDS"))
##### ARD
holediam.leafarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.leafarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.leafarea.ard.nohid, file = here("output/houwie/holediam-leafarea-ard-nohid.RDS"))


#### Hole Diameter Discrete + Petiole Length
##### ER
holediam.petleng.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.petleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.petleng.er.nohid, file = here("output/houwie/holediam-petleng-er-nohid.RDS"))
##### SYM
holediam.petleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.petleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.petleng.sym.nohid, file = here("output/houwie/holediam-petleng-sym-nohid.RDS"))
##### ARD
holediam.petleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.petleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.petleng.ard.nohid, file = here("output/houwie/holediam-petleng-ard-nohid.RDS"))


#### Hole Diameter Discrete + Stem Area
##### ER
holediam.stemarea.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.stemarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.stemarea.er.nohid, file = here("output/houwie/holediam-stemarea-er-nohid.RDS"))
##### SYM
holediam.stemarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.stemarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.stemarea.sym.nohid, file = here("output/houwie/holediam-stemarea-sym-nohid.RDS"))
##### ARD
holediam.stemarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.stemarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.stemarea.ard.nohid, file = here("output/houwie/holediam-stemarea-ard-nohid.RDS"))







###############################
#### Leaf Structure + Corola Length
##### ER
## leafstruc.corleng.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     leafstruc.corleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "leaf.struct", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(leafstruc.corleng.er.nohid, file = here("output/houwie/leafstruc-corleng-er-nohid.RDS"))
## ##### SYM
## leafstruc.corleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     leafstruc.corleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "leaf.struct", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(leafstruc.corleng.sym.nohid, file = here("output/houwie/leafstruc-corleng-sym-nohid.RDS"))
## ##### ARD
## leafstruc.corleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     leafstruc.corleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "leaf.struct", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(leafstruc.corleng.ard.nohid, file = here("output/houwie/leafstruc-corleng-ard-nohid.RDS"))

## #### Leaf Structure + Leaf Area
## ##### ER
## leafstruc.leafarea.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     leafstruc.leafarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "leaf.struct", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(leafstruc.leafarea.er.nohid, file = here("output/houwie/leafstruc-leafarea-er-nohid.RDS"))
## ##### SYM
## leafstruc.leafarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     leafstruc.leafarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "leaf.struct", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(leafstruc.leafarea.sym.nohid, file = here("output/houwie/leafstruc-leafarea-sym-nohid.RDS"))
## ##### ARD
## leafstruc.leafarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     leafstruc.leafarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "leaf.struct", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(leafstruc.leafarea.ard.nohid, file = here("output/houwie/leafstruc-leafarea-ard-nohid.RDS"))


## #### Leaf Structure + Petiole Length
## ##### ER
## leafstruc.petleng.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     leafstruc.petleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "leaf.struct", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(leafstruc.petleng.er.nohid, file = here("output/houwie/leafstruc-petleng-er-nohid.RDS"))
## ##### SYM
## leafstruc.petleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     leafstruc.petleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "leaf.struct", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(leafstruc.petleng.sym.nohid, file = here("output/houwie/leafstruc-petleng-sym-nohid.RDS"))
## ##### ARD
## leafstruc.petleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     leafstruc.petleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "leaf.struct", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(leafstruc.petleng.ard.nohid, file = here("output/houwie/leafstruc-petleng-ard-nohid.RDS"))


## #### Leaf Structure + Stem Area
## ##### ER
## leafstruc.stemarea.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     leafstruc.stemarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "leaf.struct", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(leafstruc.stemarea.er.nohid, file = here("output/houwie/leafstruc-stemarea-er-nohid.RDS"))
## ##### SYM
## leafstruc.stemarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     leafstruc.stemarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "leaf.struct", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(leafstruc.stemarea.sym.nohid, file = here("output/houwie/leafstruc-stemarea-sym-nohid.RDS"))
## ##### ARD
## leafstruc.stemarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     leafstruc.stemarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "leaf.struct", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(leafstruc.stemarea.ard.nohid, file = here("output/houwie/leafstruc-stemarea-ard-nohid.RDS"))







###############################
#### Mating System + Corola Length
##### ER
## matsys.corleng.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     matsys.corleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "mating.system", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(matsys.corleng.er.nohid, file = here("output/houwie/matsys-corleng-er-nohid.RDS"))
## ##### SYM
## matsys.corleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     matsys.corleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "mating.system", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(matsys.corleng.sym.nohid, file = here("output/houwie/matsys-corleng-sym-nohid.RDS"))
## ##### ARD
## matsys.corleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     matsys.corleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "mating.system", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(matsys.corleng.ard.nohid, file = here("output/houwie/matsys-corleng-ard-nohid.RDS"))

## #### Mating System + Leaf Area
## ##### ER
## matsys.leafarea.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     matsys.leafarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "mating.system", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(matsys.leafarea.er.nohid, file = here("output/houwie/matsys-leafarea-er-nohid.RDS"))
## ##### SYM
## matsys.leafarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     matsys.leafarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "mating.system", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(matsys.leafarea.sym.nohid, file = here("output/houwie/matsys-leafarea-sym-nohid.RDS"))
## ##### ARD
## matsys.leafarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     matsys.leafarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "mating.system", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(matsys.leafarea.ard.nohid, file = here("output/houwie/matsys-leafarea-ard-nohid.RDS"))


## #### Mating System + Petiole Length
## ##### ER
## matsys.petleng.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     matsys.petleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "mating.system", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(matsys.petleng.er.nohid, file = here("output/houwie/matsys-petleng-er-nohid.RDS"))
## ##### SYM
## matsys.petleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     matsys.petleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "mating.system", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(matsys.petleng.sym.nohid, file = here("output/houwie/matsys-petleng-sym-nohid.RDS"))
## ##### ARD
## matsys.petleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     matsys.petleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "mating.system", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(matsys.petleng.ard.nohid, file = here("output/houwie/matsys-petleng-ard-nohid.RDS"))


## #### Mating System + Stem Area
## ##### ER
## matsys.stemarea.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     matsys.stemarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "mating.system", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(matsys.stemarea.er.nohid, file = here("output/houwie/matsys-stemarea-er-nohid.RDS"))
## ##### SYM
## matsys.stemarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     matsys.stemarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "mating.system", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(matsys.stemarea.sym.nohid, file = here("output/houwie/matsys-stemarea-sym-nohid.RDS"))
## ##### ARD
## matsys.stemarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     matsys.stemarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "mating.system", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(matsys.stemarea.ard.nohid, file = here("output/houwie/matsys-stemarea-ard-nohid.RDS"))







###############################
#### Reward + Corola Length
##### ER
reward.corleng.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.corleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.corleng.er.nohid, file = here("output/houwie/reward-corleng-er-nohid.RDS"))
##### SYM
reward.corleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.corleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.corleng.sym.nohid, file = here("output/houwie/reward-corleng-sym-nohid.RDS"))
##### ARD
reward.corleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.corleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.corleng.ard.nohid, file = here("output/houwie/reward-corleng-ard-nohid.RDS"))

#### Reward + Leaf Area
##### ER
reward.leafarea.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.leafarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.leafarea.er.nohid, file = here("output/houwie/reward-leafarea-er-nohid.RDS"))
##### SYM
reward.leafarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.leafarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.leafarea.sym.nohid, file = here("output/houwie/reward-leafarea-sym-nohid.RDS"))
##### ARD
reward.leafarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.leafarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.leafarea.ard.nohid, file = here("output/houwie/reward-leafarea-ard-nohid.RDS"))


#### Reward + Petiole Length
##### ER
reward.petleng.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.petleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.petleng.er.nohid, file = here("output/houwie/reward-petleng-er-nohid.RDS"))
##### SYM
reward.petleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.petleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.petleng.sym.nohid, file = here("output/houwie/reward-petleng-sym-nohid.RDS"))
##### ARD
reward.petleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.petleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.petleng.ard.nohid, file = here("output/houwie/reward-petleng-ard-nohid.RDS"))


#### Reward + Stem Area
##### ER
reward.stemarea.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.stemarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.stemarea.er.nohid, file = here("output/houwie/reward-stemarea-er-nohid.RDS"))
##### SYM
reward.stemarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.stemarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.stemarea.sym.nohid, file = here("output/houwie/reward-stemarea-sym-nohid.RDS"))
##### ARD
reward.stemarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.stemarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.stemarea.ard.nohid, file = here("output/houwie/reward-stemarea-ard-nohid.RDS"))







###############################
#### Strategy + Corola Length
##### ER
## strategy.corleng.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     strategy.corleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "strategy", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(strategy.corleng.er.nohid, file = here("output/houwie/strategy-corleng-er-nohid.RDS"))
## ##### SYM
## strategy.corleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     strategy.corleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "strategy", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(strategy.corleng.sym.nohid, file = here("output/houwie/strategy-corleng-sym-nohid.RDS"))
## ##### ARD
## strategy.corleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "corola.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     strategy.corleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "strategy", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(strategy.corleng.ard.nohid, file = here("output/houwie/strategy-corleng-ard-nohid.RDS"))

## #### Strategy + Leaf Area
## ##### ER
## strategy.leafarea.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     strategy.leafarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "strategy", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(strategy.leafarea.er.nohid, file = here("output/houwie/strategy-leafarea-er-nohid.RDS"))
## ##### SYM
## strategy.leafarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     strategy.leafarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "strategy", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(strategy.leafarea.sym.nohid, file = here("output/houwie/strategy-leafarea-sym-nohid.RDS"))
## ##### ARD
## strategy.leafarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "leaf.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     strategy.leafarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "strategy", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(strategy.leafarea.ard.nohid, file = here("output/houwie/strategy-leafarea-ard-nohid.RDS"))


## #### Strategy + Petiole Length
## ##### ER
## strategy.petleng.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     strategy.petleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "strategy", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(strategy.petleng.er.nohid, file = here("output/houwie/strategy-petleng-er-nohid.RDS"))
## ##### SYM
## strategy.petleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     strategy.petleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "strategy", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(strategy.petleng.sym.nohid, file = here("output/houwie/strategy-petleng-sym-nohid.RDS"))
## ##### ARD
## strategy.petleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "petiole.length"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     strategy.petleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "strategy", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(strategy.petleng.ard.nohid, file = here("output/houwie/strategy-petleng-ard-nohid.RDS"))


## #### Strategy + Stem Area
## ##### ER
## strategy.stemarea.er.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     strategy.stemarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "strategy", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(strategy.stemarea.er.nohid, file = here("output/houwie/strategy-stemarea-er-nohid.RDS"))
## ##### SYM
## strategy.stemarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     strategy.stemarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "strategy", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(strategy.stemarea.sym.nohid, file = here("output/houwie/strategy-stemarea-sym-nohid.RDS"))
## ##### ARD
## strategy.stemarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
## trait.cont <- "stem.area"
## tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
## for(i in 1:length(cont.models)){
##     strategy.stemarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "strategy", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
## }
## saveRDS(strategy.stemarea.ard.nohid, file = here("output/houwie/strategy-stemarea-ard-nohid.RDS"))







###############################
#### Warts + Corola Length
##### ER
warts.corleng.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.corleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.corleng.er.nohid, file = here("output/houwie/warts-corleng-er-nohid.RDS"))
##### SYM
warts.corleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.corleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.corleng.sym.nohid, file = here("output/houwie/warts-corleng-sym-nohid.RDS"))
##### ARD
warts.corleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.corleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.corleng.ard.nohid, file = here("output/houwie/warts-corleng-ard-nohid.RDS"))

#### Warts + Leaf Area
##### ER
warts.leafarea.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.leafarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.leafarea.er.nohid, file = here("output/houwie/warts-leafarea-er-nohid.RDS"))
##### SYM
warts.leafarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.leafarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.leafarea.sym.nohid, file = here("output/houwie/warts-leafarea-sym-nohid.RDS"))
##### ARD
warts.leafarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.leafarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.leafarea.ard.nohid, file = here("output/houwie/warts-leafarea-ard-nohid.RDS"))


#### Warts + Petiole Length
##### ER
warts.petleng.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.petleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.petleng.er.nohid, file = here("output/houwie/warts-petleng-er-nohid.RDS"))
##### SYM
warts.petleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.petleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.petleng.sym.nohid, file = here("output/houwie/warts-petleng-sym-nohid.RDS"))
##### ARD
warts.petleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.petleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.petleng.ard.nohid, file = here("output/houwie/warts-petleng-ard-nohid.RDS"))


#### Warts + Stem Area
##### ER
warts.stemarea.er.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.stemarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.stemarea.er.nohid, file = here("output/houwie/warts-stemarea-er-nohid.RDS"))
##### SYM
warts.stemarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.stemarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.stemarea.sym.nohid, file = here("output/houwie/warts-stemarea-sym-nohid.RDS"))
##### ARD
warts.stemarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.stemarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 100, ncores = 56, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.stemarea.ard.nohid, file = here("output/houwie/warts-stemarea-ard-nohid.RDS"))







