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

i_am("R/houwie_analysis_hidden_state.R")

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

#fulldata <- readRDS(here("output/fulldata_houwie.RDS"))
fulldata.ingroup <- fulldata[fulldata$strategy != 0,]

fulldata.ingroup <- fulldata.ingroup[, c(1, 21:34)]
disc.traits <- sort(names(fulldata.ingroup)[c(2, 3, 7, 8, 9, 11, 13, 14, 15)])
cont.traits <- sort(names(fulldata.ingroup)[c(4, 5, 10, 12)])

cont.models <- c("BM1", "BMV", "OU1", "OUM", "OUMV", "OUMA")
disc.models <- c("ER", "SYM", "ARD")

registerDoMC(10)

### 1 hidden state

###############################
#### Domatium Growth + Corola Length
##### ER
domgrow.corleng.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.corleng.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.corleng.er.hid, file = here("output/houwie/domgrow-corleng-er-hid.RDS"))

##### SYM
domgrow.corleng.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.corleng.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.corleng.sym.hid, file = here("output/houwie/domgrow-corleng-sym-hid.RDS"))

##### ARD
domgrow.corleng.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.corleng.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.corleng.ard.hid, file = here("output/houwie/domgrow-corleng-ard-hid.RDS"))

#### Domatium Growth + Leaf Area
##### ER
domgrow.leafarea.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.leafarea.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.leafarea.er.hid, file = here("output/houwie/domgrow-leafarea-er-hid.RDS"))

##### SYM
domgrow.leafarea.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.leafarea.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.leafarea.sym.hid, file = here("output/houwie/domgrow-leafarea-sym-hid.RDS"))

##### ARD
domgrow.leafarea.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.leafarea.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.leafarea.ard.hid, file = here("output/houwie/domgrow-leafarea-ard-hid.RDS"))


#### Domatium Growth + Petiole Length
##### ER
domgrow.petleng.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.petleng.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.petleng.er.hid, file = here("output/houwie/domgrow-petleng-er-hid.RDS"))

##### SYM
domgrow.petleng.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.petleng.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.petleng.sym.hid, file = here("output/houwie/domgrow-petleng-sym-hid.RDS"))

##### ARD
domgrow.petleng.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.petleng.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.petleng.ard.hid, file = here("output/houwie/domgrow-petleng-ard-hid.RDS"))


#### Domatium Growth + Stem Area
##### ER
domgrow.stemarea.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.stemarea.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.stemarea.er.hid, file = here("output/houwie/domgrow-stemarea-er-hid.RDS"))

##### SYM
domgrow.stemarea.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.stemarea.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.stemarea.sym.hid, file = here("output/houwie/domgrow-stemarea-sym-hid.RDS"))

##### ARD
domgrow.stemarea.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    domgrow.stemarea.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(domgrow.stemarea.ard.hid, file = here("output/houwie/domgrow-stemarea-ard-hid.RDS"))







###############################
#### Hole Diameter Discrete + Corola Length
##### ER
holediam.corleng.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.corleng.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.corleng.er.hid, file = here("output/houwie/holediam-corleng-er-hid.RDS"))

##### SYM
holediam.corleng.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.corleng.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.corleng.sym.hid, file = here("output/houwie/holediam-corleng-sym-hid.RDS"))

##### ARD
holediam.corleng.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.corleng.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.corleng.ard.hid, file = here("output/houwie/holediam-corleng-ard-hid.RDS"))

#### Hole Diameter Discrete + Leaf Area
##### ER
holediam.leafarea.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.leafarea.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.leafarea.er.hid, file = here("output/houwie/holediam-leafarea-er-hid.RDS"))

##### SYM
holediam.leafarea.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.leafarea.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.leafarea.sym.hid, file = here("output/houwie/holediam-leafarea-sym-hid.RDS"))

##### ARD
holediam.leafarea.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.leafarea.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.leafarea.ard.hid, file = here("output/houwie/holediam-leafarea-ard-hid.RDS"))


#### Hole Diameter Discrete + Petiole Length
##### ER
holediam.petleng.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.petleng.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.petleng.er.hid, file = here("output/houwie/holediam-petleng-er-hid.RDS"))

##### SYM
holediam.petleng.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.petleng.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.petleng.sym.hid, file = here("output/houwie/holediam-petleng-sym-hid.RDS"))

##### ARD
holediam.petleng.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.petleng.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.petleng.ard.hid, file = here("output/houwie/holediam-petleng-ard-hid.RDS"))


#### Hole Diameter Discrete + Stem Area
##### ER
holediam.stemarea.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.stemarea.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.stemarea.er.hid, file = here("output/houwie/holediam-stemarea-er-hid.RDS"))

##### SYM
holediam.stemarea.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.stemarea.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.stemarea.sym.hid, file = here("output/houwie/holediam-stemarea-sym-hid.RDS"))

##### ARD
holediam.stemarea.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    holediam.stemarea.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(holediam.stemarea.ard.hid, file = here("output/houwie/holediam-stemarea-ard-hid.RDS"))



###############################
#### Reward + Corola Length
##### ER
reward.corleng.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.corleng.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.corleng.er.hid, file = here("output/houwie/reward-corleng-er-hid.RDS"))

##### SYM
reward.corleng.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.corleng.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.corleng.sym.hid, file = here("output/houwie/reward-corleng-sym-hid.RDS"))

##### ARD
reward.corleng.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.corleng.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.corleng.ard.hid, file = here("output/houwie/reward-corleng-ard-hid.RDS"))

#### Reward + Leaf Area
##### ER
reward.leafarea.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.leafarea.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.leafarea.er.hid, file = here("output/houwie/reward-leafarea-er-hid.RDS"))

##### SYM
reward.leafarea.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.leafarea.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.leafarea.sym.hid, file = here("output/houwie/reward-leafarea-sym-hid.RDS"))

##### ARD
reward.leafarea.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.leafarea.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.leafarea.ard.hid, file = here("output/houwie/reward-leafarea-ard-hid.RDS"))


#### Reward + Petiole Length
##### ER
reward.petleng.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.petleng.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.petleng.er.hid, file = here("output/houwie/reward-petleng-er-hid.RDS"))

##### SYM
reward.petleng.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.petleng.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.petleng.sym.hid, file = here("output/houwie/reward-petleng-sym-hid.RDS"))

##### ARD
reward.petleng.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.petleng.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.petleng.ard.hid, file = here("output/houwie/reward-petleng-ard-hid.RDS"))


#### Reward + Stem Area
##### ER
reward.stemarea.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.stemarea.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.stemarea.er.hid, file = here("output/houwie/reward-stemarea-er-hid.RDS"))

##### SYM
reward.stemarea.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.stemarea.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.stemarea.sym.hid, file = here("output/houwie/reward-stemarea-sym-hid.RDS"))

##### ARD
reward.stemarea.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    reward.stemarea.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(reward.stemarea.ard.hid, file = here("output/houwie/reward-stemarea-ard-hid.RDS"))



###############################
#### Warts + Corola Length
##### ER
warts.corleng.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.corleng.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.corleng.er.hid, file = here("output/houwie/warts-corleng-er-hid.RDS"))

##### SYM
warts.corleng.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.corleng.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.corleng.sym.hid, file = here("output/houwie/warts-corleng-sym-hid.RDS"))

##### ARD
warts.corleng.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "corola.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.corleng.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.corleng.ard.hid, file = here("output/houwie/warts-corleng-ard-hid.RDS"))

#### Warts + Leaf Area
##### ER
warts.leafarea.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.leafarea.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.leafarea.er.hid, file = here("output/houwie/warts-leafarea-er-hid.RDS"))

##### SYM
warts.leafarea.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.leafarea.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.leafarea.sym.hid, file = here("output/houwie/warts-leafarea-sym-hid.RDS"))

##### ARD
warts.leafarea.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "leaf.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.leafarea.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.leafarea.ard.hid, file = here("output/houwie/warts-leafarea-ard-hid.RDS"))


#### Warts + Petiole Length
##### ER
warts.petleng.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.petleng.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.petleng.er.hid, file = here("output/houwie/warts-petleng-er-hid.RDS"))

##### SYM
warts.petleng.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.petleng.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.petleng.sym.hid, file = here("output/houwie/warts-petleng-sym-hid.RDS"))

##### ARD
warts.petleng.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "petiole.length"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.petleng.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.petleng.ard.hid, file = here("output/houwie/warts-petleng-ard-hid.RDS"))


#### Warts + Stem Area
##### ER
warts.stemarea.er.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.stemarea.er.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.stemarea.er.hid, file = here("output/houwie/warts-stemarea-er-hid.RDS"))

##### SYM
warts.stemarea.sym.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.stemarea.sym.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.stemarea.sym.hid, file = here("output/houwie/warts-stemarea-sym-hid.RDS"))

##### ARD
warts.stemarea.ard.hid <- vector(mode = "list", length = length(cont.models))
trait.cont <- "stem.area"
tree.pruned <- drop.tip(mcc.tree, tip = mcc.tree$tip.label[!(mcc.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
for(i in 1:length(cont.models)){
    warts.stemarea.ard.hid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 10, n_starts = 56), error = function(x){NA})
}
saveRDS(warts.stemarea.ard.hid, file = here("output/houwie/warts-stemarea-ard-hid.RDS"))
