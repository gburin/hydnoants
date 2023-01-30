#library("tidyverse")
library("ape")
library("phytools")
library("foreach")
library("doMC")
library("OUwie")
library("caper")

here::i_am("R/run_ouwie_climpc_discholediam.R")

fulldata <- read.csv(here::here("data/Dataset_traits_Hydnophytinae_CODED_2022.csv"))

names(fulldata)[21:ncol(fulldata)] <- c("strategy", "warts", "stem.area", "leaf.area", "hole.diameter", "architecture", "dom.growth", "reward", "corola.length", "mating.system", "petiole.length", "leaf.struct", "appendages", "holediam.disc")

fulldata$species[grep("Wor2g", fulldata$species)] <- gsub("Wor2g", "Worthing", fulldata$species[grep("Wor2g", fulldata$species)])
fulldata$species[grep("areolata", fulldata$species)] <- gsub("areolata", "aerolata", fulldata$species[grep("areolata", fulldata$species)])
fulldata$hole.diameter[which(fulldata$hole.diameter == "NA (outgroup has no domatia)")] <- NA
fulldata$hole.diameter <- as.numeric(fulldata$hole.diameter)


mcc.tree <- ape::read.tree(here::here("data/MCCtree135taxa.tre"))

clim.pca <- prcomp(fulldata[, 2:20])
summary(clim.pca)

fulldata <- cbind(fulldata, clim.pca$x[, 1:3])
names(fulldata)[(ncol(fulldata) - 2):ncol(fulldata)] <- paste0("climpc", 1:3)

data.pgls <- comparative.data(mcc.tree, fulldata, names.col = "species", na.omit = FALSE)


## Building models to check for significant relationships

### Leaf Area
summary(lm.leafarea.pc1 <- lm(data.pgls$data$leaf.area ~ data.pgls$data$climpc1))
summary(lm.leafarea.pc2 <- lm(data.pgls$data$leaf.area ~ data.pgls$data$climpc2))
summary(lm.leafarea.pc3 <- lm(data.pgls$data$leaf.area ~ data.pgls$data$climpc3))

summary(pgls.leafarea.pc1 <- pgls(leaf.area ~ climpc1, data.pgls, lambda = "ML"))
summary(pgls.leafarea.pc2 <- pgls(leaf.area ~ climpc2, data.pgls, lambda = "ML"))
summary(pgls.leafarea.pc3 <- pgls(leaf.area ~ climpc3, data.pgls, lambda = "ML"))

### Stem Area
summary(lm.stemarea.pc1 <- lm(data.pgls$data$stem.area ~ data.pgls$data$climpc1))
summary(lm.stemarea.pc2 <- lm(data.pgls$data$stem.area ~ data.pgls$data$climpc2))
summary(lm.stemarea.pc3 <- lm(data.pgls$data$stem.area ~ data.pgls$data$climpc3)) ## p = 0.04825

summary(pgls.stemarea.pc1 <- pgls(stem.area ~ climpc1, data.pgls, lambda = "ML"))
summary(pgls.stemarea.pc2 <- pgls(stem.area ~ climpc2, data.pgls, lambda = "ML"))
summary(pgls.stemarea.pc3 <- pgls(stem.area ~ climpc3, data.pgls, lambda = "ML"))

### Hole Diameter
summary(lm.hole.diameter.pc1 <- lm(data.pgls$data$hole.diameter ~ data.pgls$data$climpc1)) ## p = 0.01772
summary(lm.hole.diameter.pc2 <- lm(data.pgls$data$hole.diameter ~ data.pgls$data$climpc2)) ## p = 0.03268
summary(lm.hole.diameter.pc3 <- lm(data.pgls$data$hole.diameter ~ data.pgls$data$climpc3))

summary(pgls.hole.diameter.pc1 <- pgls(hole.diameter ~ climpc1, data.pgls, lambda = "ML"))
summary(pgls.hole.diameter.pc2 <- pgls(hole.diameter ~ climpc2, data.pgls, lambda = "ML")) ## p = 0.1422
summary(pgls.hole.diameter.pc3 <- pgls(hole.diameter ~ climpc3, data.pgls, lambda = "ML"))

### Corola Length
summary(lm.corola.length.pc1 <- lm(data.pgls$data$corola.length ~ data.pgls$data$climpc1))
summary(lm.corola.length.pc2 <- lm(data.pgls$data$corola.length ~ data.pgls$data$climpc2))
summary(lm.corola.length.pc3 <- lm(data.pgls$data$corola.length ~ data.pgls$data$climpc3)) ## p = 0.4697

summary(pgls.corola.length.pc1 <- pgls(corola.length ~ climpc1, data.pgls, lambda = "ML"))
summary(pgls.corola.length.pc2 <- pgls(corola.length ~ climpc2, data.pgls, lambda = "ML"))
summary(pgls.corola.length.pc3 <- pgls(corola.length ~ climpc3, data.pgls, lambda = "ML"))

### Petiole Length
summary(lm.petiole.length.pc1 <- lm(data.pgls$data$petiole.length ~ data.pgls$data$climpc1))
summary(lm.petiole.length.pc2 <- lm(data.pgls$data$petiole.length ~ data.pgls$data$climpc2))
summary(lm.petiole.length.pc3 <- lm(data.pgls$data$petiole.length ~ data.pgls$data$climpc3))

summary(pgls.petiole.length.pc1 <- pgls(petiole.length ~ climpc1, data.pgls, lambda = "ML"))
summary(pgls.petiole.length.pc2 <- pgls(petiole.length ~ climpc2, data.pgls, lambda = "ML"))
summary(pgls.petiole.length.pc3 <- pgls(petiole.length ~ climpc3, data.pgls, lambda = "ML"))



fulldata$stem.area.pc3.lm <- summary(lm.stemarea.pc3)$residuals

fulldata$hole.diameter.pc1.lm[which(!is.na(fulldata$hole.diameter))] <- summary(lm.hole.diameter.pc1)$residuals
fulldata$hole.diameter.pc2.lm[which(!is.na(fulldata$hole.diameter))] <- summary(lm.hole.diameter.pc2)$residuals
fulldata$hole.diameter.pc2.pgls[which(!is.na(fulldata$hole.diameter))] <- summary(pgls.hole.diameter.pc2)$residuals[,1]

fulldata$corola.length.pc3.lm[which(!is.na(fulldata$corola.length))] <- summary(lm.corola.length.pc3)$residuals


## Keeping data as originally coded for SIMMAP

traits.simmap <- fulldata

## Recoding trait states for plotting

app <- c("outgroup", "none", "variable", "spines")
arch <- c("shrub", "multiple", "single")
domgrow <- c("outgroup", "diffuse", "apical")
leafstruc <- c("thick", "variable", "thin", "succulent")
matsys <- c("heterostylous", "non-heterostylous", "funct_unisexual")
reward <- c("outgroup", "absent", "present")
strategy <- c("outgroup", "facultative", "obligate", "lost")
warts <- c("absent", "variable", "differentiated", "lost")
holediam.disc <- c("absent", "severallargebase_smallside", "onelargebase_smalldise", "alllarge")

fulldata$strategy <- strategy[fulldata$strategy + 1]
fulldata$warts <- warts[fulldata$warts + 1]
fulldata$reward <- reward[fulldata$reward + 1]
fulldata$architecture <- arch[fulldata$architecture + 1]
fulldata$dom.growth <- domgrow[fulldata$dom.growth + 1]
fulldata$mating.system <- matsys[fulldata$mating.system + 1]
fulldata$leaf.struct <- leafstruc[fulldata$leaf.struct + 1]
fulldata$appendages <- app[fulldata$appendages + 1]
fulldata$holediam.disc <- holediam.disc[fulldata$holediam.disc + 1]

fulldata.names <- fulldata[,1]
fulldata <- fulldata[,-1]
rownames(fulldata) <- fulldata.names

discrete.traits <- c("holediam.disc")
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
        make.simmap(tree = mcc.tree, x = factor(setNames(traits.simmap[, discrete.traits[i]], traits.simmap$species)), model = best.models[i], nsim = 1, Q = "empirical", pi = "estimated", .parallel = TRUE)}
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
