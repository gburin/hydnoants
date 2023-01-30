#library("tidyverse")
library("ape")
library("phytools")
library("foreach")
library("doMC")
library("OUwie")
here::i_am("R/data_preparation_chelsea.R")

fulldata <- read.csv(here::here("data/Dataset_traits_Hydnophytinae_CODED_apr2021.csv"))

names(fulldata)[21:ncol(fulldata)] <- c("strategy", "warts", "stem.area", "leaf.area", "hole.diameter", "architecture", "dom.growth", "reward", "corola.length", "mating.system", "petiole.length", "leaf.struct", "appendages")

fulldata$species[grep("Wor2g", fulldata$species)] <- gsub("Wor2g", "Worthing", fulldata$species[grep("Wor2g", fulldata$species)])
fulldata$species[grep("areolata", fulldata$species)] <- gsub("areolata", "aerolata", fulldata$species[grep("areolata", fulldata$species)])
fulldata$hole.diameter[which(fulldata$hole.diameter == "NA (outgroup has no domatia)")] <- NA
fulldata$hole.diameter <- as.numeric(fulldata$hole.diameter)


mcc.tree <- ape::read.tree(here::here("data/MCCtree135taxa.tre"))

clim.pca <- prcomp(fulldata[, 2:20])
summary(clim.pca)

fulldata <- cbind(fulldata, clim.pca$x)
rownames(fulldata) <- fulldata$species

fulldata$leaf.area.pc1 <- nlme::gls(leaf.area ~ PC1, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$leaf.area.pc1, "std") <- NULL
fulldata$leaf.area.pc2 <- nlme::gls(leaf.area ~ PC2, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$leaf.area.pc2, "std") <- NULL
fulldata$leaf.area.pc3 <- nlme::gls(leaf.area ~ PC3, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$leaf.area.pc3, "std") <- NULL

fulldata$stem.area.pc1 <- nlme::gls(stem.area ~ PC1, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$stem.area.pc1, "std") <- NULL
fulldata$stem.area.pc2 <- nlme::gls(stem.area ~ PC2, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$stem.area.pc2, "std") <- NULL
fulldata$stem.area.pc3 <- nlme::gls(stem.area ~ PC3, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$stem.area.pc3, "std") <- NULL

fulldata$hole.diameter.pc1 <- nlme::gls(hole.diameter ~ PC1, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$hole.diameter.pc1, "std") <- NULL
fulldata$hole.diameter.pc2 <- nlme::gls(hole.diameter ~ PC2, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$hole.diameter.pc2, "std") <- NULL
fulldata$hole.diameter.pc3 <- nlme::gls(hole.diameter ~ PC3, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$hole.diameter.pc3, "std") <- NULL

fulldata$corola.length.pc1 <- nlme::gls(corola.length ~ PC1, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$corola.length.pc1, "std") <- NULL
fulldata$corola.length.pc2 <- nlme::gls(corola.length ~ PC2, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$corola.length.pc2, "std") <- NULL
fulldata$corola.length.pc3 <- nlme::gls(corola.length ~ PC3, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$corola.length.pc3, "std") <- NULL

fulldata$petiole.length.pc1 <- nlme::gls(petiole.length ~ PC1, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$petiole.length.pc1, "std") <- NULL
fulldata$petiole.length.pc2 <- nlme::gls(petiole.length ~ PC2, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$petiole.length.pc2, "std") <- NULL
fulldata$petiole.length.pc3 <- nlme::gls(petiole.length ~ PC3, data = fulldata, correlation = corBrownian(phy = mcc.tree), method = "ML")$residuals
attr(fulldata$petiole.length.pc3, "std") <- NULL


## Keeping data as originally coded for SIMMAP

traits.simmap <- fulldata

## Recoding trait states for plotting

strategy <- c("outgroup", "facultative", "obligate", "lost")
warts <- c("absent", "variable", "differentiated", "lost")
reward <- c("outgroup", "absent", "present")
archit <- c("shrub", "multiple", "single")
dom.growth <- c("outgroup", "diffuse", "apical")
mating.system <- c("heterostylous", "non-heterostylous", "funct_unisexual")
leaf.struct <- c("thick", "variable", "thin", "succulent")
appendages <- c("outgroup", "none", "variable", "spines")

fulldata$strategy <- strategy[fulldata$strategy + 1]
fulldata$warts <- warts[fulldata$warts + 1]
fulldata$reward <- reward[fulldata$reward + 1]
fulldata$architecture <- archit[fulldata$architecture + 1]
fulldata$dom.growth <- dom.growth[fulldata$dom.growth + 1]
fulldata$mating.system <- mating.system[fulldata$mating.system + 1]
fulldata$leaf.struct <- leaf.struct[fulldata$leaf.struct + 1]
fulldata$appendages <- appendages[fulldata$appendages + 1]

fulldata.names <- fulldata[,1]
fulldata <- fulldata[,-1]
rownames(fulldata) <- fulldata.names

discrete.traits <- c("strategy", "warts", "reward", "architecture", "dom.growth", "mating.system", "leaf.struct", "appendages")
# Strategy as discrete character

best.models <- vector(mode = "character", length = length(discrete.traits))

for(i in 1:length(best.models)){
    print(i)
    ard <- make.simmap(tree = mcc.tree, x = factor(setNames(traits.simmap[, discrete.traits[i]], traits.simmap$species)), model = "ARD", nsim = 1, Q = "empirical", pi = "estimated", .parallel = TRUE)
    sym <- make.simmap(tree = mcc.tree, x = factor(setNames(traits.simmap[, discrete.traits[i]], traits.simmap$species)), model = "SYM", nsim = 1, Q = "empirical", pi = "estimated", .parallel = TRUE)
    er <- make.simmap(tree = mcc.tree, x = factor(setNames(traits.simmap[, discrete.traits[i]], traits.simmap$species)), model = "ER", nsim = 1, Q = "empirical", pi = "estimated", .parallel = TRUE)
    aic.tab <- AIC(er, sym, ard)
    best.models[i] <- stringr::str_to_upper(rownames(aic.tab)[which.min(aic.tab$AIC)])
}

registerDoMC(50)

for(i in 1:length(discrete.traits)){
    print(i)
    raw.simmap <- foreach(j = 1:1000) %dopar% {
        make.simmap(tree = mcc.tree, x = factor(setNames(traits.simmap[, discrete.traits[i]], traits.simmap$species)), model = best.models[i], nsim = 1, Q = "empirical", pi = "estimated", .parallel = TRUE)
    simmap.temp <- lapply(raw.simmap, drop.tip.simmap, tip = traits.simmap$species[traits.simmap$strategy == 0])
    class(simmap.temp) <- c("multiPhylo", "multiSimmap")
    saveRDS(simmap.temp, file = paste0("../output/simmaps_", discrete.traits[i], ".RDS"))
    }


## Generating Stochastic Maps
## Code here is a little "clumsy" because it is written to be self-dependent and also to run in parallel. Let me know if you have problems understanding it!


for(disc in 1:length(discrete.traits)){

    ## registerDoMC(56)
    assign(paste0("simmap.list.", discrete.traits[disc]), readRDS(paste0("../output/simmaps_", discrete.traits[disc], ".RDS")))
    traits.simmap.in <- traits.simmap[traits.simmap$strategy != 0,]

    ## Fitting trait evolution models

    models <- c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")
    for(k in 1:3){
        
        assign(paste0("fit.mcc.", discrete.traits[disc], ".stemarea.pc", k), vector(mode = "list", length = 7))
        assign(paste0("fit.mcc.", discrete.traits[disc], ".leafarea.pc", k), vector(mode = "list", length = 7))
        assign(paste0("fit.mcc.", discrete.traits[disc], ".corleng.pc", k), vector(mode = "list", length = 7))
        assign(paste0("fit.mcc.", discrete.traits[disc], ".petleng.pc", k), vector(mode = "list", length = 7))
        assign(paste0("fit.mcc.", discrete.traits[disc], ".holediam.pc", k), vector(mode = "list", length = 7))

        simmap.sample <- sort(sample(1:1000, 100))
        saveRDS(simmap.sample, file = paste0("../output/simmap_sample_", discrete.traits[disc], "_pc", k, ".RDS"))
        
### Stem area

        assign(paste0("tree.stem.", discrete.traits[disc]), lapply(get(paste0("simmap.list.", discrete.traits[disc])), drop.tip, tip = traits.simmap.in$species[is.na(traits.simmap.in$stem.area)]))
        for(j in 1:length(models)){
            res.stemarea <- foreach(i = 1:length(simmap.sample)) %dopar% {
                print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Stem Area"))
                tryCatch(OUwie(phy = eval(parse(text = paste0("tree.stem.", discrete.traits[disc], "[[simmap.sample[", i, "]]]"))), data = traits.simmap.in[!is.na(traits.simmap.in$stem.area), c("species", discrete.traits[disc], paste0("stem.area.pc", k))], model = models[j], simmap.tree = TRUE), error = function(x){NA})
            }
            eval(parse(text = paste0("fit.mcc.", discrete.traits[disc], ".stemarea.pc", k, "[[", j, "]] <- res.stemarea")))
        }

### Leaf area

        assign(paste0("tree.leaf.", discrete.traits[disc]), lapply(get(paste0("simmap.list.", discrete.traits[disc])), drop.tip, tip = traits.simmap.in$species[is.na(traits.simmap.in$leaf.area)]))
        for(j in 1:length(models)){
            res.leafarea <- foreach(i = 1:length(simmap.sample)) %dopar% {
                print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Leaf Area"))
                tryCatch(OUwie(phy = eval(parse(text = paste0("tree.leaf.", discrete.traits[disc], "[[simmap.sample[", i, "]]]"))), data = traits.simmap.in[!is.na(traits.simmap.in$leaf.area), c("species", discrete.traits[disc], paste0("leaf.area.pc", k))], model = models[j], simmap.tree = TRUE), error = function(x){NA})
            }
            eval(parse(text = paste0("fit.mcc.", discrete.traits[disc], ".leafarea.pc", k, "[[", j, "]] <- res.leafarea")))
        }


### Corola length

        assign(paste0("tree.corola.", discrete.traits[disc]), lapply(get(paste0("simmap.list.", discrete.traits[disc])), drop.tip, tip = traits.simmap.in$species[is.na(traits.simmap.in$corola.length)]))
        for(j in 1:length(models)){
            res.corleng <- foreach(i = 1:length(simmap.sample)) %dopar% {
                print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Corola Length"))
                tryCatch(OUwie(phy = eval(parse(text = paste0("tree.corola.", discrete.traits[disc], "[[simmap.sample[", i, "]]]"))), data = traits.simmap.in[!is.na(traits.simmap.in$corola.length), c("species", discrete.traits[disc], paste0("corola.length.pc", k))], model = models[j], simmap.tree = TRUE), error = function(x){NA})
            }
            eval(parse(text = paste0("fit.mcc.", discrete.traits[disc], ".corleng.pc", k, "[[", j, "]] <- res.corleng")))
        }

### Petiole length

        assign(paste0("tree.petiole.", discrete.traits[disc]), lapply(get(paste0("simmap.list.", discrete.traits[disc])), drop.tip, tip = traits.simmap.in$species[is.na(traits.simmap.in$petiole.length)]))
        for(j in 1:length(models)){
            res.petleng <- foreach(i = 1:length(simmap.sample)) %dopar% {
                print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Petiole Length"))
                tryCatch(OUwie(phy = eval(parse(text = paste0("tree.petiole.", discrete.traits[disc], "[[simmap.sample[", i, "]]]"))), data = traits.simmap.in[!is.na(traits.simmap.in$petiole.length), c("species", discrete.traits[disc], paste0("petiole.length.pc", k))], model = models[j], simmap.tree = TRUE), error = function(x){NA})
            }
            eval(parse(text = paste0("fit.mcc.", discrete.traits[disc], ".petleng.pc", k, "[[", j, "]] <- res.petleng")))
        }

### Hole diameter

        assign(paste0("tree.hole.", discrete.traits[disc]), lapply(get(paste0("simmap.list.", discrete.traits[disc])), drop.tip, tip = traits.simmap.in$species[is.na(traits.simmap.in$hole.diameter)]))
        for(j in 1:length(models)){
            res.holediam <- foreach(i = 1:length(simmap.sample)) %dopar% {
                print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Hole Diameter"))
                tryCatch(OUwie(phy = eval(parse(text = paste0("tree.hole.", discrete.traits[disc], "[[simmap.sample[", i, "]]]"))), data = traits.simmap.in[!is.na(traits.simmap.in$hole.diameter), c("species", discrete.traits[disc], paste0("hole.diameter.pc", k))], model = models[j], simmap.tree = TRUE), error = function(x){NA})
            }
            eval(parse(text = paste0("fit.mcc.", discrete.traits[disc], ".holediam.pc", k, "[[", j, "]] <- res.holediam")))
        }

        saveRDS(list(simmap.sample,
                     get(paste0("simmap.list.", discrete.traits[disc])),
                     get(paste0("fit.mcc.", discrete.traits[disc], ".stemarea.pc", k)),
                     get(paste0("fit.mcc.", discrete.traits[disc], ".leafarea.pc", k)),
                     get(paste0("fit.mcc.", discrete.traits[disc], ".corleng.pc", k)),
                     get(paste0("fit.mcc.", discrete.traits[disc], ".petleng.pc", k)),
                     get(paste0("fit.mcc.", discrete.traits[disc], ".holediam.pc", k))),
                file = paste0("../output/ouwie_mcc_", discrete.traits[disc], "_fit_residuals_pc", k, ".RDS"))
        rm(list = ls(pattern = "fit.mcc."))
    }
}

## save.image("../output/ouwie_mcc_full_fit_residuals.RData")
