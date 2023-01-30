library("ape")
library("phytools")
library("OUwie")
library("tidyverse")
library("plyr")
library("phangorn")
library("caper")
library("picante")
library("foreach")
library("doMC")

## Importing trees

mcc.tree <- read.tree("../data/MCCtree135taxa.tre")
mcc.tree <- force.ultrametric(mcc.tree, method = "nnls")

traits <- read.csv("../data/Dataset_traits_Hydnophytinae_CODED_nov2020.csv", as.is = TRUE)
names(traits)[1:14] <- c("species", "strategy", "warts", "stem.area", "leaf.area", "hole.diameter", "architecture", "dom.growth", "reward", "corola.length", "mating.system", "petiole.length", "leaf.struct", "appendages")

## Fixing typo in species name

traits$species[grep("Wor2g", traits$species)] <- gsub("Wor2g", "Worthing", traits$species[grep("Wor2g", traits$species)])
traits$species[grep("areolata", traits$species)] <- gsub("areolata", "aerolata", traits$species[grep("areolata", traits$species)])
traits$hole.diameter[which(traits$hole.diameter == "NA (outgroup has no domatia)")] <- NA
traits$hole.diameter <- as.numeric(traits$hole.diameter)

## Keeping data as originally coded for SIMMAP

traits.simmap <- traits

## Recoding trait states for plotting

strategy <- c("outgroup", "facultative", "obligate", "lost")
warts <- c("absent", "variable", "differentiated", "lost")
reward <- c("outgroup", "absent", "present")
archit <- c("shrub", "multiple", "single")
dom.growth <- c("outgroup", "diffuse", "apical")
mating.system <- c("heterostylous", "non-heterostylous", "funct_unisexual")
leaf.struct <- c("thick", "variable", "thin", "succulent")
appendages <- c("outgroup", "none", "variable", "spines")

traits$strategy <- strategy[traits$strategy + 1]
traits$warts <- warts[traits$warts + 1]
traits$reward <- reward[traits$reward + 1]
traits$architecture <- archit[traits$architecture + 1]
traits$dom.growth <- dom.growth[traits$dom.growth + 1]
traits$mating.system <- mating.system[traits$mating.system + 1]
traits$leaf.struct <- leaf.struct[traits$leaf.struct + 1]
traits$appendages <- appendages[traits$appendages + 1]

traits.names <- traits[,1]
traits <- traits[,-1]
rownames(traits) <- traits.names

## Removing outgroup species


discrete.traits <- c("strategy", "warts", "reward", "architecture", "dom.growth", "mating.system", "leaf.struct", "appendages")
# Strategy as discrete character

## Generating Stochastic Maps
## Code here is a little "clumsy" because it is written to be self-dependent and also to run in parallel. Let me know if you have problems understanding it!

registerDoMC(56)

for(disc in 1:length(discrete.traits)){    
    raw.simmap <- foreach(a = 1:1000) %dopar% {
        make.simmap(tree = mcc.tree, x = factor(setNames(traits.simmap[, discrete.traits[disc]], traits.simmap$species)), model = "ARD", nsim = 1, Q = "mcmc", pi = "estimated", .parallel = TRUE)
    }
    simmap.temp <- lapply(raw.simmap, drop.tip.simmap, tip = traits.simmap$species[traits.simmap$strategy == 0])
    class(simmap.temp) <- c("multiPhylo", "multiSimmap")
    assign(paste0("simmap.list.", discrete.traits[disc]), simmap.temp)

    traits.simmap.in <- traits.simmap[traits.simmap$strategy != 0,]

    #mcc.tree.in <- drop.tip(mcc.tree, traits.simmap$species[traits.simmap$strategy == 0])

    ## Plotting Resulting Maps

    ## colors <- setNames(brewer.pal(3, "Set1"), 1:3)
    ## sum.simmap <- summary(simmap.list.strategy)
    ## plot(sum.simmap, colors = colors)
    ## add.simmap.legend(leg = c("facultative", "obligate", "lost"), colors = colors, x = 2, y = 4, prompt = FALSE)

    ## Fitting trait evolution models

    models <- c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")

    assign(paste0("fit.mcc.", discrete.traits[disc], ".stemarea"), vector(mode = "list", length = 7))
    assign(paste0("fit.mcc.", discrete.traits[disc], ".leafarea"), vector(mode = "list", length = 7))
    assign(paste0("fit.mcc.", discrete.traits[disc], ".corleng"), vector(mode = "list", length = 7))
    assign(paste0("fit.mcc.", discrete.traits[disc], ".petleng"), vector(mode = "list", length = 7))
    assign(paste0("fit.mcc.", discrete.traits[disc], ".holediam"), vector(mode = "list", length = 7))

    simmap.sample <- sort(sample(1:1000, 100))

    registerDoMC(56)
    
### Stem area

    for(j in 1:length(models)){
        res.stemarea <- foreach(i = 1:length(simmap.sample)) %dopar% {
            print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Stem Area"))
            tryCatch(OUwie(phy = eval(parse(text = paste0("simmap.list.", discrete.traits[disc], "[[", simmap.sample[i], "]]"))), data = traits.simmap.in[, c("species", discrete.traits[disc], "stem.area")], model = models[j], simmap.tree = TRUE), error = function(x){NA})
        }
        eval(parse(text = paste0("fit.mcc.", discrete.traits[disc], ".stemarea", "[[", j, "]] <- res.stemarea")))
    }

### Leaf area

    for(j in 1:length(models)){
        res.leafarea <- foreach(i = 1:length(simmap.sample)) %dopar% {
            print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Leaf Area"))
            tryCatch(OUwie(phy = eval(parse(text = paste0("simmap.list.", discrete.traits[disc], "[[", simmap.sample[i], "]]"))), data = traits.simmap.in[, c("species", discrete.traits[disc], "leaf.area")], model = models[j], simmap.tree = TRUE), error = function(x){NA})
        }
        eval(parse(text = paste0("fit.mcc.", discrete.traits[disc], ".leafarea", "[[", j, "]] <- res.leafarea")))
    }

### Corola length

    assign(paste0("tree.corola.", discrete.traits[disc]), lapply(get(paste0("simmap.list.", discrete.traits[disc])), drop.tip, tip = traits.simmap.in$species[is.na(traits.simmap.in$corola.length)]))
    for(j in 1:length(models)){
        res.corleng <- foreach(i = 1:length(simmap.sample)) %dopar% {
            print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Corola Length"))
            tryCatch(OUwie(phy = eval(parse(text = paste0("tree.corola.", discrete.traits[disc], "[[simmap.sample[", i, "]]]"))), data = traits.simmap.in[!is.na(traits.simmap.in$corola.length), c("species", discrete.traits[disc], "corola.length")], model = models[j], simmap.tree = TRUE), error = function(x){NA})
        }
        eval(parse(text = paste0("fit.mcc.", discrete.traits[disc], ".corleng", "[[", j, "]] <- res.corleng")))
    }

### Petiole length

    for(j in 1:length(models)){
        res.petleng <- foreach(i = 1:length(simmap.sample)) %dopar% {
            print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Stem Area"))
            tryCatch(OUwie(phy = eval(parse(text = paste0("simmap.list.", discrete.traits[disc], "[[", simmap.sample[i], "]]"))), data = traits.simmap.in[, c("species", discrete.traits[disc], "petiole.length")], model = models[j], simmap.tree = TRUE), error = function(x){NA})
        }
        eval(parse(text = paste0("fit.mcc.", discrete.traits[disc], ".petleng", "[[", j, "]] <- res.petleng")))
    }

### Hole diameter

    for(j in 1:length(models)){
        res.holediam <- foreach(i = 1:length(simmap.sample)) %dopar% {
            print(paste(models[j], "from SIMMAP", i, "for", discrete.traits[disc], "+ Stem Area"))
            tryCatch(OUwie(phy = eval(parse(text = paste0("simmap.list.", discrete.traits[disc], "[[", simmap.sample[i], "]]"))), data = traits.simmap.in[, c("species", discrete.traits[disc], "hole.diameter")], model = models[j], simmap.tree = TRUE), error = function(x){NA})
        }
        eval(parse(text = paste0("fit.mcc.", discrete.traits[disc], ".holediam", "[[", j, "]] <- res.holediam")))
    }

    saveRDS(list(simmap.sample,
                 get(paste0("simmap.list.", discrete.traits[disc])),
                 get(paste0("fit.mcc.", discrete.traits[disc], ".stemarea")),
                 get(paste0("fit.mcc.", discrete.traits[disc], ".leafarea")),
                 get(paste0("fit.mcc.", discrete.traits[disc], ".corleng")),
                 get(paste0("fit.mcc.", discrete.traits[disc], ".petleng")),
                 get(paste0("fit.mcc.", discrete.traits[disc], ".holediam"))),
            file = paste0("../output/ouwie_mcc_", discrete.traits[disc], "_fit.RDS"))
}

save.image("../output/ouwie_mcc_full_fit.RData")
