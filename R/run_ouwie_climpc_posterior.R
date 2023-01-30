#library("tidyverse")
library("ape")
library("phytools")
library("foreach")
library("doMC")
library("OUwie")
library("caper")

here::i_am("R/run_ouwie_climpc_posterior.R")

fulldata <- read.csv(here::here("data/Dataset_traits_Hydnophytinae_CODED_2022.csv"))

names(fulldata)[21:ncol(fulldata)] <- c("strategy", "warts", "stem.area", "leaf.area", "hole.diameter", "architecture", "dom.growth", "reward", "corola.length", "mating.system", "petiole.length", "leaf.struct", "appendages", "holediam.disc")

fulldata$species[grep("Wor2g", fulldata$species)] <- gsub("Wor2g", "Worthing", fulldata$species[grep("Wor2g", fulldata$species)])
fulldata$species <- gsub("aerolata", "areolata", fulldata$species)
fulldata$hole.diameter[which(fulldata$hole.diameter == "NA (outgroup has no domatia)")] <- NA
fulldata$hole.diameter <- as.numeric(fulldata$hole.diameter)


post.trees <- read.nexus(here::here("data/posterior_trimmed_135taxa.tre"))

clim.pca <- prcomp(fulldata[, 2:20])
## summary(clim.pca)

fulldata <- cbind(fulldata, clim.pca$x[, 1:3])
names(fulldata)[(ncol(fulldata) - 2):ncol(fulldata)] <- paste0("climpc", 1:3)


## sampled.trees <- sort(sample(1:length(post.trees), 20))
## saveRDS(sampled.trees, file = here::here("output/sampled_trees_posterior.RDS"))

sampled.trees <- readRDS(here::here("output/sampled_trees_posterior.RDS"))

for(i in 1:length(sampled.trees)){
    print(i)

    data.pgls <- comparative.data(post.trees[[sampled.trees[i]]], fulldata, names.col = "species", na.omit = FALSE)

    ## Building models to check for significant relationships

### Leaf Area
    fulldata[which(!is.na(fulldata$leaf.area)), paste0("leafarea.pc1.tree", i)] <- tryCatch(tryCatch(summary(pgls(leaf.area ~ climpc1, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(leaf.area ~ climpc1, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})
    fulldata[which(!is.na(fulldata$leaf.area)), paste0("leafarea.pc2.tree", i)] <- tryCatch(tryCatch(summary(pgls(leaf.area ~ climpc2, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(leaf.area ~ climpc2, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})
    fulldata[which(!is.na(fulldata$leaf.area)), paste0("leafarea.pc3.tree", i)] <- tryCatch(tryCatch(summary(pgls(leaf.area ~ climpc3, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(leaf.area ~ climpc3, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})

### Stem Area
    fulldata[which(!is.na(fulldata$stem.area)), paste0("stemarea.pc1.tree", i)] <- tryCatch(tryCatch(summary(pgls(stem.area ~ climpc1, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(stem.area ~ climpc1, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})
    fulldata[which(!is.na(fulldata$stem.area)), paste0("stemarea.pc2.tree", i)] <- tryCatch(tryCatch(summary(pgls(stem.area ~ climpc2, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(stem.area ~ climpc2, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})
    fulldata[which(!is.na(fulldata$stem.area)), paste0("stemarea.pc3.tree", i)] <- tryCatch(tryCatch(summary(pgls(stem.area ~ climpc3, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(stem.area ~ climpc3, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})

### Hole Diameter
    fulldata[which(!is.na(fulldata$hole.diameter)), paste0("holediam.pc1.tree", i)] <- tryCatch(tryCatch(summary(pgls(hole.diameter ~ climpc1, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(hole.diameter ~ climpc1, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})
    fulldata[which(!is.na(fulldata$hole.diameter)), paste0("holediam.pc2.tree", i)] <- tryCatch(tryCatch(summary(pgls(hole.diameter ~ climpc2, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(hole.diameter ~ climpc2, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})
    fulldata[which(!is.na(fulldata$hole.diameter)), paste0("holediam.pc3.tree", i)] <- tryCatch(tryCatch(summary(pgls(hole.diameter ~ climpc3, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(hole.diameter ~ climpc3, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})

### Corola Length
    fulldata[which(!is.na(fulldata$corola.length)), paste0("corola.length.pc1.tree", i)] <- tryCatch(tryCatch(summary(pgls(corola.length ~ climpc1, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(corola.length ~ climpc1, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})
    fulldata[which(!is.na(fulldata$corola.length)), paste0("corola.length.pc2.tree", i)] <- tryCatch(tryCatch(summary(pgls(corola.length ~ climpc2, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(corola.length ~ climpc2, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})
    fulldata[which(!is.na(fulldata$corola.length)), paste0("corola.length.pc3.tree", i)] <- tryCatch(tryCatch(summary(pgls(corola.length ~ climpc3, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(corola.length ~ climpc3, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})

### Petiole Length
    fulldata[which(!is.na(fulldata$petiole.length)), paste0("petiole.length.pc1.tree", i)] <- tryCatch(tryCatch(summary(pgls(petiole.length ~ climpc1, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(petiole.length ~ climpc1, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})
    fulldata[which(!is.na(fulldata$petiole.length)), paste0("petiole.length.pc2.tree", i)] <- tryCatch(tryCatch(summary(pgls(petiole.length ~ climpc2, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(petiole.length ~ climpc2, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})
    fulldata[which(!is.na(fulldata$petiole.length)), paste0("petiole.length.pc3.tree", i)] <- tryCatch(tryCatch(summary(pgls(petiole.length ~ climpc3, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 1))))$residuals, error = summary(pgls(petiole.length ~ climpc3, data.pgls, lambda = "ML", bounds = list(lambda = c(0.1, 0.5))))$residuals), error = function(x){NA})

}


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

fulldata$strategy <- strategy[fulldata$strategy + 1]
fulldata$warts <- warts[fulldata$warts + 1]
fulldata$reward <- reward[fulldata$reward + 1]
fulldata$architecture <- arch[fulldata$architecture + 1]
fulldata$dom.growth <- domgrow[fulldata$dom.growth + 1]
fulldata$mating.system <- matsys[fulldata$mating.system + 1]
fulldata$leaf.struct <- leafstruc[fulldata$leaf.struct + 1]
fulldata$appendages <- app[fulldata$appendages + 1]

fulldata.names <- fulldata[,1]
fulldata <- fulldata[,-1]
rownames(fulldata) <- fulldata.names

discrete.traits <- c("strategy", "warts", "reward", "architecture", "dom.growth", "mating.system", "leaf.struct", "appendages", "holediam.disc")
# Strategy as discrete character

registerDoMC(length(sampled.trees))

best.models <- foreach(i = 1:length(sampled.trees)) %dopar% {
    print(i)
    res <- list()
    for(j in 1:length(discrete.traits)){
        ard <- make.simmap(tree = post.trees[[sampled.trees[i]]], x = factor(setNames(traits.simmap[, discrete.traits[j]], traits.simmap$species)), model = "ARD", nsim = 1, Q = "empirical", pi = "estimated", .parallel = TRUE)
        sym <- make.simmap(tree = post.trees[[sampled.trees[i]]], x = factor(setNames(traits.simmap[, discrete.traits[j]], traits.simmap$species)), model = "SYM", nsim = 1, Q = "empirical", pi = "estimated", .parallel = TRUE)
        er <- make.simmap(tree = post.trees[[sampled.trees[i]]], x = factor(setNames(traits.simmap[, discrete.traits[j]], traits.simmap$species)), model = "ER", nsim = 1, Q = "empirical", pi = "estimated", .parallel = TRUE)
        aic.tab <- AIC(er, sym, ard)
        res[[j]] <- stringr::str_to_upper(rownames(aic.tab)[which.min(aic.tab$AIC)])
    }
    res
}

registerDoMC(50)

for(k in 1:length(sampled.trees)){
    for(i in 1:length(discrete.traits)){
        print(paste0("Tree ", k, " Trait ", i))
        raw.simmap <- foreach(j = 1:100) %dopar% {
            make.simmap(tree = post.trees[[sampled.trees[k]]], x = factor(setNames(traits.simmap[, discrete.traits[i]], traits.simmap$species)), model = best.models[[k]][[i]], nsim = 1, Q = "empirical", pi = "estimated", .parallel = TRUE)
            }
        simmap.temp <- lapply(raw.simmap, drop.tip.simmap, tip = traits.simmap$species[traits.simmap$strategy == 0])
        class(simmap.temp) <- c("multiPhylo", "multiSimmap")
        saveRDS(simmap.temp, file = here::here(paste0("output/posterior_trees/simmaps_tree_", k, "_", discrete.traits[i], ".RDS")))
    }
}

## Generating Stochastic Maps
## Code here is a little "clumsy" because it is written to be self-dependent and also to run in parallel. Let me know if you have problems understanding it!


for(disc in 1:length(discrete.traits)){
    for(p in 1:length(sampled.trees)){

        ## registerDoMC(56)
        assign(paste0("simmap.list.tree.", p, ".", discrete.traits[disc]), readRDS(here::here(paste0("output/posterior_trees/simmaps_tree_", p, "_", discrete.traits[disc], ".RDS"))))
        traits.simmap.in <- traits.simmap[traits.simmap$strategy != 0,]

        ## Fitting trait evolution models

        models <- c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")
        for(k in 1:3){
            
            assign(paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".stemarea.pc", k), vector(mode = "list", length = 7))
            assign(paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".leafarea.pc", k), vector(mode = "list", length = 7))
            assign(paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".corleng.pc", k), vector(mode = "list", length = 7))
            assign(paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".petleng.pc", k), vector(mode = "list", length = 7))
            assign(paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".holediam.pc", k), vector(mode = "list", length = 7))

            simmap.sample <- sort(sample(1:100, 10))
            saveRDS(simmap.sample, file = here::here(paste0("output/posterior_trees/simmap_sample_tree_", p, "_", discrete.traits[disc], "_pc", k, ".RDS")))
            
### Stem area

            assign(paste0("tree", p, ".stem.", discrete.traits[disc]), lapply(get(paste0("simmap.list.tree.", p, ".", discrete.traits[disc])), drop.tip, tip = traits.simmap.in$species[is.na(traits.simmap.in$stem.area)]))
            for(j in 1:length(models)){
                res.stemarea <- foreach(i = 1:length(simmap.sample)) %dopar% {
                    print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Stem Area"))
                    tryCatch(OUwie(phy = eval(parse(text = paste0("tree", p, ".stem.", discrete.traits[disc])))[[i]], data = traits.simmap.in[!is.na(traits.simmap.in$stem.area), c("species", discrete.traits[disc], paste0("stemarea.pc", k, ".tree", p))], model = models[j], simmap.tree = TRUE), error = function(x){NA})
                }
                eval(parse(text = paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".stemarea.pc", k, "[[", j, "]] <- res.stemarea")))
            }

### Leaf area

            assign(paste0("tree", p, ".leaf.", discrete.traits[disc]), lapply(get(paste0("simmap.list.tree.", p, ".", discrete.traits[disc])), drop.tip, tip = traits.simmap.in$species[is.na(traits.simmap.in$leaf.area)]))
            for(j in 1:length(models)){
                res.leafarea <- foreach(i = 1:length(simmap.sample)) %dopar% {
                    print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Leaf Area"))
                    tryCatch(OUwie(phy = eval(parse(text = paste0("tree", p, ".leaf.", discrete.traits[disc])))[[i]], data = traits.simmap.in[!is.na(traits.simmap.in$leaf.area), c("species", discrete.traits[disc], paste0("leafarea.pc", k, ".tree", p))], model = models[j], simmap.tree = TRUE), error = function(x){NA})
                }
                eval(parse(text = paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".leafarea.pc", k, "[[", j, "]] <- res.leafarea")))
            }


### Corola length

            assign(paste0("tree", p, ".corola.", discrete.traits[disc]), lapply(get(paste0("simmap.list.tree.", p, ".", discrete.traits[disc])), drop.tip, tip = traits.simmap.in$species[is.na(traits.simmap.in$corola.length)]))
            for(j in 1:length(models)){
                res.corola.length <- foreach(i = 1:length(simmap.sample)) %dopar% {
                    print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Corola Length"))
                    tryCatch(OUwie(phy = eval(parse(text = paste0("tree", p, ".corola.", discrete.traits[disc])))[[i]], data = traits.simmap.in[!is.na(traits.simmap.in$corola.length), c("species", discrete.traits[disc], paste0("corola.length.pc", k, ".tree", p))], model = models[j], simmap.tree = TRUE), error = function(x){NA})
                }
                eval(parse(text = paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".corleng.pc", k, "[[", j, "]] <- res.corola.length")))
            }


### Petiole length

            assign(paste0("tree", p, ".petiole.", discrete.traits[disc]), lapply(get(paste0("simmap.list.tree.", p, ".", discrete.traits[disc])), drop.tip, tip = traits.simmap.in$species[is.na(traits.simmap.in$petiole.length)]))
            for(j in 1:length(models)){
                res.petiole.length <- foreach(i = 1:length(simmap.sample)) %dopar% {
                    print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Petiole Length"))
                    tryCatch(OUwie(phy = eval(parse(text = paste0("tree", p, ".petiole.", discrete.traits[disc])))[[i]], data = traits.simmap.in[!is.na(traits.simmap.in$petiole.length), c("species", discrete.traits[disc], paste0("petiole.length.pc", k, ".tree", p))], model = models[j], simmap.tree = TRUE), error = function(x){NA})
                }
                eval(parse(text = paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".petleng.pc", k, "[[", j, "]] <- res.petiole.length")))
            }

### Hole diameter

            assign(paste0("tree", p, ".hole.", discrete.traits[disc]), lapply(get(paste0("simmap.list.tree.", p, ".", discrete.traits[disc])), drop.tip, tip = traits.simmap.in$species[is.na(traits.simmap.in$hole.diameter)]))
            for(j in 1:length(models)){
                res.hole.diameter <- foreach(i = 1:length(simmap.sample)) %dopar% {
                    print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Hole Length"))
                    tryCatch(OUwie(phy = eval(parse(text = paste0("tree", p, ".hole.", discrete.traits[disc])))[[i]], data = traits.simmap.in[!is.na(traits.simmap.in$hole.diameter), c("species", discrete.traits[disc], paste0("holediam.pc", k, ".tree", p))], model = models[j], simmap.tree = TRUE), error = function(x){NA})
                }
                eval(parse(text = paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".holediam.pc", k, "[[", j, "]] <- res.hole.diameter")))
            }

            saveRDS(list(simmap.sample,
                         get(paste0("simmap.list.tree.", p, ".", discrete.traits[disc]))[[i]],
                         get(paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".stemarea.pc", k)),
                         get(paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".leafarea.pc", k)),
                         get(paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".corleng.pc", k)),
                         get(paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".petleng.pc", k)),
                         get(paste0("fit.mcc.tree.", p, ".", discrete.traits[disc], ".holediam.pc", k))),
                    file = here::here(paste0("output/posterior_trees/ouwie_mcc_tree", p, "_", discrete.traits[disc], "_fit_residuals_pc", k, ".RDS")))
            rm(list = ls(pattern = "fit.mcc.tree."))
        }
    }
}

## save.image(here::here("output/posterior_trees/ouwie_mcc_full_fit_residuals.RData")
