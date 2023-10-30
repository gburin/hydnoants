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

i_am("R/houwie_analysis_posterior.R")

fulldata <- read.csv(here("data/Dataset_traits_Hydnophytinae_CODED_2022.csv"))

names(fulldata)[21:ncol(fulldata)] <- c("strategy", "warts", "stem.area", "leaf.area", "hole.diameter", "architecture", "dom.growth", "reward", "corola.length", "mating.system", "petiole.length", "leaf.struct", "appendages", "holediam.disc")

fulldata$species[grep("Wor2g", fulldata$species)] <- gsub("Wor2g", "Worthing", fulldata$species[grep("Wor2g", fulldata$species)])
fulldata$species[grep("areolata", fulldata$species)] <- gsub("areolata", "aerolata", fulldata$species[grep("areolata", fulldata$species)])
fulldata$hole.diameter[which(fulldata$hole.diameter == "NA (outgroup has no domatia)")] <- NA
fulldata$hole.diameter <- as.numeric(fulldata$hole.diameter)


## mcc.tree <- ape::read.tree(here("data/MCCtree135taxa.tre"))
trees.sample <- readRDS(here("output/sampled_trees_posterior.RDS"))
posterior.trees <- read.nexus(here("data/posterior_trimmed_135taxa.tre"))

#fulldata <- readRDS(here("output/fulldata_houwie.RDS"))
fulldata.ingroup <- fulldata[fulldata$strategy != 0,]

fulldata.ingroup <- fulldata.ingroup[, c(1, 21:34)]
disc.traits <- sort(names(fulldata.ingroup)[c(2, 3, 7, 8, 9, 11, 13, 14, 15)])
cont.traits <- sort(names(fulldata.ingroup)[c(4, 5, 10, 12)])

cont.models <- c("BM1", "BMV", "OU1", "OUM", "OUMV", "OUMA")
disc.models <- c("ER", "SYM", "ARD")

registerDoMC(50)

### No hidden states

###############################
#### Domatium Growth + Corola Length
##### ER
for(i = 1:20) %dopar% {
    domgrow.corleng.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "corola.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        domgrow.corleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(domgrow.corleng.er.nohid, file = here(paste0("output/houwie/domgrow-corleng-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    domgrow.corleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "corola.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        domgrow.corleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(domgrow.corleng.sym.nohid, file = here(paste0("output/houwie/domgrow-corleng-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    domgrow.corleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "corola.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        domgrow.corleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(domgrow.corleng.ard.nohid, file = here(paste0("output/houwie/domgrow-corleng-ard-nohid-tree", i, ".RDS")))
}

#### Domatium Growth + Leaf Area
##### ER
for(i = 1:20) %dopar% {
    domgrow.leafarea.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "leaf.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        domgrow.leafarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(domgrow.leafarea.er.nohid, file = here(paste0("output/houwie/domgrow-leafarea-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    domgrow.leafarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "leaf.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        domgrow.leafarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(domgrow.leafarea.sym.nohid, file = here(paste0("output/houwie/domgrow-leafarea-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    domgrow.corleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "leaf.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        domgrow.corleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(domgrow.corleng.ard.nohid, file = here(paste0("output/houwie/domgrow-corleng-ard-nohid-tree", i, ".RDS")))
}


#### Domatium Growth + Petiole Length
##### ER
for(i = 1:20) %dopar% {
    domgrow.petleng.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "petiole.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        domgrow.petleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(domgrow.petleng.er.nohid, file = here(paste0("output/houwie/domgrow-petleng-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    domgrow.petleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "petiole.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        domgrow.petleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(domgrow.petleng.sym.nohid, file = here(paste0("output/houwie/domgrow-petleng-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    domgrow.petleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "petiole.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        domgrow.petleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(domgrow.petleng.ard.nohid, file = here(paste0("output/houwie/domgrow-petleng-ard-nohid-tree", i, ".RDS")))
}


#### Domatium Growth + Stem Area
##### ER
for(i = 1:20) %dopar% {
    domgrow.stemarea.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "stem.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        domgrow.stemarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(domgrow.stemarea.er.nohid, file = here(paste0("output/houwie/domgrow-stemarea-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    domgrow.stemarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "stem.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        domgrow.stemarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(domgrow.stemarea.sym.nohid, file = here(paste0("output/houwie/domgrow-stemarea-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    domgrow.stemarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "stem.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        domgrow.stemarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "dom.growth", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(domgrow.stemarea.ard.nohid, file = here(paste0("output/houwie/domgrow-stemarea-ard-nohid-tree", i, ".RDS")))
}







###############################
#### Hole Diameter Discrete + Corola Length
##### ER
for(i = 1:20) %dopar% {
    holediam.corleng.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "corola.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        holediam.corleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(holediam.corleng.er.nohid, file = here(paste0("output/houwie/holediam-corleng-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    holediam.corleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "corola.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        holediam.corleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(holediam.corleng.sym.nohid, file = here(paste0("output/houwie/holediam-corleng-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    holediam.corleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "corola.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        holediam.corleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(holediam.corleng.ard.nohid, file = here(paste0("output/houwie/holediam-corleng-ard-nohid-tree", i, ".RDS")))
}

#### Hole Diameter Discrete + Leaf Area
##### ER
for(i = 1:20) %dopar% {
    holediam.leafarea.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "leaf.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        holediam.leafarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(holediam.leafarea.er.nohid, file = here(paste0("output/houwie/holediam-leafarea-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    holediam.leafarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "leaf.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        holediam.leafarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(holediam.leafarea.sym.nohid, file = here(paste0("output/houwie/holediam-leafarea-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    holediam.leafarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "leaf.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        holediam.leafarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(holediam.leafarea.ard.nohid, file = here(paste0("output/houwie/holediam-leafarea-ard-nohid-tree", i, ".RDS")))
}


#### Hole Diameter Discrete + Petiole Length
##### ER
for(i = 1:20) %dopar% {
    holediam.petleng.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "petiole.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        holediam.petleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(holediam.petleng.er.nohid, file = here(paste0("output/houwie/holediam-petleng-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    holediam.petleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "petiole.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        holediam.petleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(holediam.petleng.sym.nohid, file = here(paste0("output/houwie/holediam-petleng-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    holediam.petleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "petiole.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        holediam.petleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(holediam.petleng.ard.nohid, file = here(paste0("output/houwie/holediam-petleng-ard-nohid-tree", i, ".RDS")))
}


#### Hole Diameter Discrete + Stem Area
##### ER
for(i = 1:20) %dopar% {
    holediam.stemarea.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "stem.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        holediam.stemarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(holediam.stemarea.er.nohid, file = here(paste0("output/houwie/holediam-stemarea-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    holediam.stemarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "stem.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        holediam.stemarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(holediam.stemarea.sym.nohid, file = here(paste0("output/houwie/holediam-stemarea-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    holediam.stemarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "stem.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        holediam.stemarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "holediam.disc", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(holediam.stemarea.ard.nohid, file = here(paste0("output/houwie/holediam-stemarea-ard-nohid-tree", i, ".RDS")))
}




###############################
#### Reward + Corola Length
##### ER
for(i = 1:20) %dopar% {
    reward.corleng.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "corola.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        reward.corleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(reward.corleng.er.nohid, file = here(paste0("output/houwie/reward-corleng-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    reward.corleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "corola.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        reward.corleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(reward.corleng.sym.nohid, file = here(paste0("output/houwie/reward-corleng-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    reward.corleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "corola.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        reward.corleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(reward.corleng.ard.nohid, file = here(paste0("output/houwie/reward-corleng-ard-nohid-tree", i, ".RDS")))
}

#### Reward + Leaf Area
##### ER
for(i = 1:20) %dopar% {
    reward.leafarea.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "leaf.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        reward.leafarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(reward.leafarea.er.nohid, file = here(paste0("output/houwie/reward-leafarea-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    reward.leafarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "leaf.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        reward.leafarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(reward.leafarea.sym.nohid, file = here(paste0("output/houwie/reward-leafarea-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    reward.leafarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "leaf.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        reward.leafarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(reward.leafarea.ard.nohid, file = here(paste0("output/houwie/reward-leafarea-ard-nohid-tree", i, ".RDS")))
}


#### Reward + Petiole Length
##### ER
for(i = 1:20) %dopar% {
    reward.petleng.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "petiole.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        reward.petleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(reward.petleng.er.nohid, file = here(paste0("output/houwie/reward-petleng-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    reward.petleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "petiole.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        reward.petleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(reward.petleng.sym.nohid, file = here(paste0("output/houwie/reward-petleng-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    reward.petleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "corola.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        reward.petleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(reward.petleng.ard.nohid, file = here(paste0("output/houwie/reward-petleng-ard-nohid-tree", i, ".RDS")))
}


#### Reward + Stem Area
##### ER
for(i = 1:20) %dopar% {
    reward.stemarea.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "stem.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        reward.stemarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(reward.stemarea.er.nohid, file = here(paste0("output/houwie/reward-stemarea-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    reward.stemarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "stem.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        reward.stemarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(reward.stemarea.sym.nohid, file = here(paste0("output/houwie/reward-stemarea-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    reward.stemarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "stem.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        reward.stemarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "reward", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(reward.stemarea.ard.nohid, file = here(paste0("output/houwie/reward-stemarea-ard-nohid-tree", i, ".RDS")))
}



###############################
#### Warts + Corola Length
##### ER
for(i = 1:20) %dopar% {
    warts.corleng.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "corola.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        warts.corleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(warts.corleng.er.nohid, file = here(paste0("output/houwie/warts-corleng-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    warts.corleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "corola.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        warts.corleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(warts.corleng.sym.nohid, file = here(paste0("output/houwie/warts-corleng-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    warts.corleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "corola.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        warts.corleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(warts.corleng.ard.nohid, file = here(paste0("output/houwie/warts-corleng-ard-nohid-tree", i, ".RDS")))
}

#### Warts + Leaf Area
##### ER
for(i = 1:20) %dopar% {
    warts.leafarea.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "leaf.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        warts.leafarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(warts.leafarea.er.nohid, file = here(paste0("output/houwie/warts-leafarea-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    warts.leafarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "leaf.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        warts.leafarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(warts.leafarea.sym.nohid, file = here(paste0("output/houwie/warts-leafarea-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    warts.leafarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "leaf.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        warts.leafarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(warts.leafarea.ard.nohid, file = here(paste0("output/houwie/warts-leafarea-ard-nohid-tree", i, ".RDS")))
}


#### Warts + Petiole Length
##### ER
for(i = 1:20) %dopar% {
    warts.petleng.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "petiole.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        warts.petleng.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(warts.petleng.er.nohid, file = here(paste0("output/houwie/warts-petleng-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    warts.petleng.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "petiole.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        warts.petleng.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(warts.petleng.sym.nohid, file = here(paste0("output/houwie/warts-petleng-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    warts.petleng.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "petiole.length"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        warts.petleng.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(warts.petleng.ard.nohid, file = here(paste0("output/houwie/warts-petleng-ard-nohid-tree", i, ".RDS")))
}


#### Warts + Stem Area
##### ER
for(i = 1:20) %dopar% {
    warts.stemarea.er.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "stem.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        warts.stemarea.er.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ER", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(warts.stemarea.er.nohid, file = here(paste0("output/houwie/warts-stemarea-er-nohid-tree", i, ".RDS")))
}

##### SYM
for(i = 1:20) %dopar% {
    warts.stemarea.sym.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "stem.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        warts.stemarea.sym.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "SYM", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(warts.stemarea.sym.nohid, file = here(paste0("output/houwie/warts-stemarea-sym-nohid-tree", i, ".RDS")))
}

##### ARD
for(i = 1:20) %dopar% {
    warts.stemarea.ard.nohid <- vector(mode = "list", length = length(cont.models))
    trait.cont <- "stem.area"
    sampled.tree <- posterior.trees[[trees.sample[i]]]
    tree.pruned <- drop.tip(sampled.tree, tip = sampled.tree$tip.label[!(sampled.tree$tip.label %in% fulldata.ingroup$species[!is.na(fulldata.ingroup[, trait.cont])])])
    for(i in 1:length(cont.models)){
        warts.stemarea.ard.nohid[[i]] <- tryCatch(hOUwie(tree.pruned, na.omit(fulldata.ingroup[, c("species", "warts", trait.cont)]), rate.cat = 1, discrete_model = "ARD", continuous_model = cont.models[i], nSim = 10, ncores = 2, n_starts = 10), error = function(x){NA})
    }
    saveRDS(warts.stemarea.ard.nohid, file = here(paste0("output/houwie/warts-stemarea-ard-nohid-tree", i, ".RDS")))
}
