library("phytools")
library("RColorBrewer")
library("plyr")
library("patchwork")
library("tidyverse")

here::i_am("R/ouwie_summary_mcc_climpc_posterior.R")

## Reorganizing results

par.extr <- function(sol){
    if(is.na(sol$solution[1, 1])){
        if(ncol(sol$solution) == 1){
            sigma <- rep(sol$solution[2,], 4)
            alpha <- rep(NA, 4)
            theta <- rep(NA, 4)
        } else {
            sigma <- c(sol$solution[2,], rep(NA, 4 - ncol(sol$solution)))
            alpha <- rep(NA, 4)
            theta <- rep(NA, 4)
        }
    } else {
        ifelse(length(unique(sol$solution[2,])) == 1, assign("sigma", rep(unique(sol$solution[2,]), 4)), assign("sigma", c(sol$solution[2,], rep(NA, 4 - ncol(sol$solution)))))
        ifelse(length(unique(sol$solution[1,])) == 1, assign("alpha", rep(unique(sol$solution[1,]), 4)), assign("alpha", c(sol$solution[1,], rep(NA, 4 - ncol(sol$solution)))))
        ifelse(length(unique(sol$theta[,1])) == 1, assign("theta", rep(unique(sol$theta[,1]), 4)), assign("theta", c(sol$theta[,1], rep(NA, 4 - nrow(sol$theta)))))
    }
    res <- c(sigma, alpha, theta)
    names(res) <- paste0(rep(c("sigma", "alpha", "theta"), each = 4), rep(1:4, 3))
    return(res)
}

result.reorg <- function(x, pc, nsimmap){
    traits <- c("stemarea", "leafarea", "corleng", "petleng", "holediam")
    data <- get(paste0("mcc.", x, ".pc", pc))
    res <- vector(mode = "list", length = nsimmap)
    for(k in 1:length(traits)){
        for(i in 1:nsimmap){
            res[[i]] <- vector(mode = "list", length = 6)
            for(j in 1:6){
                res[[i]][[j]] <- data[[(k + 2)]][[j]][[i]]
            }
        }
        assign(paste0(traits[k], ".res"), res)
    }
    return(list(stemarea.res, leafarea.res, corleng.res, petleng.res, holediam.res))
}

par.averag <- function(res, nsimmap){
    pars <- llply(res, function(y){suppressMessages(bind_rows(llply(y, par.extr)))})
    aicc <- ldply(res, function(y){t(ldply(y, "[[", "AICc"))})
    daicc <- apply(aicc, 1, function(x){x - min(x)})
    aicc.w <- exp(-0.5 * daicc)
    mod.av <- as.data.frame(matrix(NA, ncol = 12, nrow = nsimmap))
    names(mod.av) <- paste0(rep(c("alpha", "sigma", "theta"), each = 4), 1:4)
    for(i in 1:nsimmap){
        pars[[i]]$aic.w <- aicc.w[, i]/sum(aicc.w[, i])
        temp <- suppressMessages(bind_cols(llply(pars[[i]][, 1:(ncol(pars[[i]]) - 1)], function(x){x * pars[[i]][, ncol(pars[[i]])]})))
        mod.av[i, match(names(pars[[i]])[-ncol(pars[[i]])], names(mod.av))] <- apply(temp, 2, sum, na.rm = TRUE)
    }
    return(mod.av)
}






#### PC1
## Appendages
list.app.pc1.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.app.pc1 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_appendages_fit_residuals_pc1.RDS")))
    names(mcc.app.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    app.pc1 <- result.reorg("app", pc = 1, nsimmap = 10)
    temp.app.pc1.av <- vector(mode = "list", length = length(app.pc1))
    for(j in 1:length(app.pc1)){
        temp.app.pc1.av[[j]] <- tryCatch(par.averag(app.pc1[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.app.pc1.av <- bind_rows(temp.app.pc1.av)
    list.app.pc1.av[[i]] <- cbind(temp.app.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.app.pc1.av <- bind_rows(list.app.pc1.av)
post.app.pc1.av$tree <- rep(1:20, each = 50)

## Architecture
list.arch.pc1.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.arch.pc1 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_architecture_fit_residuals_pc1.RDS")))
    names(mcc.arch.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    arch.pc1 <- result.reorg("arch", pc = 1, nsimmap = 10)
    temp.arch.pc1.av <- vector(mode = "list", length = length(arch.pc1))
    for(j in 1:length(arch.pc1)){
        temp.arch.pc1.av[[j]] <- tryCatch(par.averag(arch.pc1[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.arch.pc1.av <- bind_rows(temp.arch.pc1.av)
    list.arch.pc1.av[[i]] <- cbind(temp.arch.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.arch.pc1.av <- bind_rows(list.arch.pc1.av)
post.arch.pc1.av$tree <- rep(1:20, each = 50)

## Domatium Growth
list.domgrow.pc1.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.domgrow.pc1 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_dom.growth_fit_residuals_pc1.RDS")))
    names(mcc.domgrow.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    domgrow.pc1 <- result.reorg("domgrow", pc = 1, nsimmap = 10)
    temp.domgrow.pc1.av <- vector(mode = "list", length = length(domgrow.pc1))
    for(j in 1:length(domgrow.pc1)){
        temp.domgrow.pc1.av[[j]] <- tryCatch(par.averag(domgrow.pc1[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.domgrow.pc1.av <- bind_rows(temp.domgrow.pc1.av)
    list.domgrow.pc1.av[[i]] <- cbind(temp.domgrow.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.domgrow.pc1.av <- bind_rows(list.domgrow.pc1.av)
post.domgrow.pc1.av$tree <- rep(1:20, each = 50)

## Leaf Structure
list.leafstruc.pc1.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.leafstruc.pc1 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_leaf.struct_fit_residuals_pc1.RDS")))
    names(mcc.leafstruc.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    leafstruc.pc1 <- result.reorg("leafstruc", pc = 1, nsimmap = 10)
    temp.leafstruc.pc1.av <- vector(mode = "list", length = length(leafstruc.pc1))
    for(j in 1:length(leafstruc.pc1)){
        temp.leafstruc.pc1.av[[j]] <- tryCatch(par.averag(leafstruc.pc1[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.leafstruc.pc1.av <- bind_rows(temp.leafstruc.pc1.av)
    list.leafstruc.pc1.av[[i]] <- cbind(temp.leafstruc.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.leafstruc.pc1.av <- bind_rows(list.leafstruc.pc1.av)
post.leafstruc.pc1.av$tree <- rep(1:20, each = 50)


## Mating System
list.matsys.pc1.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.matsys.pc1 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_mating.system_fit_residuals_pc1.RDS")))
    names(mcc.matsys.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    matsys.pc1 <- result.reorg("matsys", pc = 1, nsimmap = 10)
    temp.matsys.pc1.av <- vector(mode = "list", length = length(matsys.pc1))
    for(j in 1:length(matsys.pc1)){
        temp.matsys.pc1.av[[j]] <- tryCatch(par.averag(matsys.pc1[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.matsys.pc1.av <- bind_rows(temp.matsys.pc1.av)
    list.matsys.pc1.av[[i]] <- cbind(temp.matsys.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.matsys.pc1.av <- bind_rows(list.matsys.pc1.av)
post.matsys.pc1.av$tree <- rep(1:20, each = 50)

## Reward
list.reward.pc1.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.reward.pc1 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_reward_fit_residuals_pc1.RDS")))
    names(mcc.reward.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    reward.pc1 <- result.reorg("reward", pc = 1, nsimmap = 10)
    temp.reward.pc1.av <- vector(mode = "list", length = length(reward.pc1))
    for(j in 1:length(reward.pc1)){
        temp.reward.pc1.av[[j]] <- tryCatch(par.averag(reward.pc1[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.reward.pc1.av <- bind_rows(temp.reward.pc1.av)
    list.reward.pc1.av[[i]] <- cbind(temp.reward.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.reward.pc1.av <- bind_rows(list.reward.pc1.av)
post.reward.pc1.av$tree <- rep(1:20, each = 50)

## Strategy
list.strategy.pc1.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.strategy.pc1 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_strategy_fit_residuals_pc1.RDS")))
    names(mcc.strategy.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    strategy.pc1 <- result.reorg("strategy", pc = 1, nsimmap = 10)
    temp.strategy.pc1.av <- vector(mode = "list", length = length(strategy.pc1))
    for(j in 1:length(strategy.pc1)){
        temp.strategy.pc1.av[[j]] <- tryCatch(par.averag(strategy.pc1[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.strategy.pc1.av <- bind_rows(temp.strategy.pc1.av)
    list.strategy.pc1.av[[i]] <- cbind(temp.strategy.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.strategy.pc1.av <- bind_rows(list.strategy.pc1.av)
post.strategy.pc1.av$tree <- rep(1:20, each = 50)

## Warts
list.warts.pc1.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.warts.pc1 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_warts_fit_residuals_pc1.RDS")))
    names(mcc.warts.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    warts.pc1 <- result.reorg("warts", pc = 1, nsimmap = 10)
    temp.warts.pc1.av <- vector(mode = "list", length = length(warts.pc1))
    for(j in 1:length(warts.pc1)){
        temp.warts.pc1.av[[j]] <- tryCatch(par.averag(warts.pc1[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.warts.pc1.av <- bind_rows(temp.warts.pc1.av)
    list.warts.pc1.av[[i]] <- cbind(temp.warts.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.warts.pc1.av <- bind_rows(list.warts.pc1.av)
post.warts.pc1.av$tree <- rep(1:20, each = 50)

## Hole Diameter - Discrete
list.holediam.disc.pc1.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.holediam.disc.pc1 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_holediam.disc_fit_residuals_pc1.RDS")))
    names(mcc.holediam.disc.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    holediam.disc.pc1 <- result.reorg("holediam.disc", pc = 1, nsimmap = 10)
    temp.holediam.disc.pc1.av <- vector(mode = "list", length = length(holediam.disc.pc1))
    for(j in 1:length(holediam.disc.pc1)){
        temp.holediam.disc.pc1.av[[j]] <- tryCatch(par.averag(holediam.disc.pc1[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.holediam.disc.pc1.av <- bind_rows(temp.holediam.disc.pc1.av)
    list.holediam.disc.pc1.av[[i]] <- cbind(temp.holediam.disc.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.holediam.disc.pc1.av <- bind_rows(list.holediam.disc.pc1.av)
post.holediam.disc.pc1.av$tree <- rep(1:20, each = 50)


fullpar.table.pc1 <- bind_rows(post.app.pc1.av, post.arch.pc1.av, post.domgrow.pc1.av, post.leafstruc.pc1.av, post.matsys.pc1.av, post.reward.pc1.av, post.strategy.pc1.av, post.warts.pc1.av, post.holediam.disc.pc1.av)
fullpar.table.pc1$disc.trait <- rep(c("app", "arch", "domgrow", "leafstruc", "matsys", "reward", "strategy", "warts", "holediam.disc"), each = 100)

write.table(fullpar.table.pc1, here::here("output/fullpar_table_climpc1_posterior.csv"), sep = ",", quote = FALSE, row.names = FALSE)

beepr::beep("fanfare")

rm(list = ls(pattern = "pc1"))
#rm(fullpar.table.pc1)





#### PC2
## Appendages
list.app.pc2.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.app.pc2 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_appendages_fit_residuals_pc2.RDS")))
    names(mcc.app.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    app.pc2 <- result.reorg("app", pc = 2, nsimmap = 10)
    temp.app.pc2.av <- vector(mode = "list", length = length(app.pc2))
    for(j in 1:length(app.pc2)){
        temp.app.pc2.av[[j]] <- tryCatch(par.averag(app.pc2[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.app.pc2.av <- bind_rows(temp.app.pc2.av)
    list.app.pc2.av[[i]] <- cbind(temp.app.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.app.pc2.av <- bind_rows(list.app.pc2.av)
post.app.pc2.av$tree <- rep(1:20, each = 50)

## Architecture
list.arch.pc2.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.arch.pc2 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_architecture_fit_residuals_pc2.RDS")))
    names(mcc.arch.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    arch.pc2 <- result.reorg("arch", pc = 2, nsimmap = 10)
    temp.arch.pc2.av <- vector(mode = "list", length = length(arch.pc2))
    for(j in 1:length(arch.pc2)){
        temp.arch.pc2.av[[j]] <- tryCatch(par.averag(arch.pc2[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.arch.pc2.av <- bind_rows(temp.arch.pc2.av)
    list.arch.pc2.av[[i]] <- cbind(temp.arch.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.arch.pc2.av <- bind_rows(list.arch.pc2.av)
post.arch.pc2.av$tree <- rep(1:20, each = 50)

## Domatium Growth
list.domgrow.pc2.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.domgrow.pc2 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_dom.growth_fit_residuals_pc2.RDS")))
    names(mcc.domgrow.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    domgrow.pc2 <- result.reorg("domgrow", pc = 2, nsimmap = 10)
    temp.domgrow.pc2.av <- vector(mode = "list", length = length(domgrow.pc2))
    for(j in 1:length(domgrow.pc2)){
        temp.domgrow.pc2.av[[j]] <- tryCatch(par.averag(domgrow.pc2[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.domgrow.pc2.av <- bind_rows(temp.domgrow.pc2.av)
    list.domgrow.pc2.av[[i]] <- cbind(temp.domgrow.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.domgrow.pc2.av <- bind_rows(list.domgrow.pc2.av)
post.domgrow.pc2.av$tree <- rep(1:20, each = 50)

## Leaf Structure
list.leafstruc.pc2.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.leafstruc.pc2 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_leaf.struct_fit_residuals_pc2.RDS")))
    names(mcc.leafstruc.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    leafstruc.pc2 <- result.reorg("leafstruc", pc = 2, nsimmap = 10)
    temp.leafstruc.pc2.av <- vector(mode = "list", length = length(leafstruc.pc2))
    for(j in 1:length(leafstruc.pc2)){
        temp.leafstruc.pc2.av[[j]] <- tryCatch(par.averag(leafstruc.pc2[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.leafstruc.pc2.av <- bind_rows(temp.leafstruc.pc2.av)
    list.leafstruc.pc2.av[[i]] <- cbind(temp.leafstruc.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.leafstruc.pc2.av <- bind_rows(list.leafstruc.pc2.av)
post.leafstruc.pc2.av$tree <- rep(1:20, each = 50)


## Mating System
list.matsys.pc2.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.matsys.pc2 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_mating.system_fit_residuals_pc2.RDS")))
    names(mcc.matsys.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    matsys.pc2 <- result.reorg("matsys", pc = 2, nsimmap = 10)
    temp.matsys.pc2.av <- vector(mode = "list", length = length(matsys.pc2))
    for(j in 1:length(matsys.pc2)){
        temp.matsys.pc2.av[[j]] <- tryCatch(par.averag(matsys.pc2[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.matsys.pc2.av <- bind_rows(temp.matsys.pc2.av)
    list.matsys.pc2.av[[i]] <- cbind(temp.matsys.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.matsys.pc2.av <- bind_rows(list.matsys.pc2.av)
post.matsys.pc2.av$tree <- rep(1:20, each = 50)

## Reward
list.reward.pc2.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.reward.pc2 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_reward_fit_residuals_pc2.RDS")))
    names(mcc.reward.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    reward.pc2 <- result.reorg("reward", pc = 2, nsimmap = 10)
    temp.reward.pc2.av <- vector(mode = "list", length = length(reward.pc2))
    for(j in 1:length(reward.pc2)){
        temp.reward.pc2.av[[j]] <- tryCatch(par.averag(reward.pc2[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.reward.pc2.av <- bind_rows(temp.reward.pc2.av)
    list.reward.pc2.av[[i]] <- cbind(temp.reward.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.reward.pc2.av <- bind_rows(list.reward.pc2.av)
post.reward.pc2.av$tree <- rep(1:20, each = 50)

## Strategy
list.strategy.pc2.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.strategy.pc2 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_strategy_fit_residuals_pc2.RDS")))
    names(mcc.strategy.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    strategy.pc2 <- result.reorg("strategy", pc = 2, nsimmap = 10)
    temp.strategy.pc2.av <- vector(mode = "list", length = length(strategy.pc2))
    for(j in 1:length(strategy.pc2)){
        temp.strategy.pc2.av[[j]] <- tryCatch(par.averag(strategy.pc2[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.strategy.pc2.av <- bind_rows(temp.strategy.pc2.av)
    list.strategy.pc2.av[[i]] <- cbind(temp.strategy.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.strategy.pc2.av <- bind_rows(list.strategy.pc2.av)
post.strategy.pc2.av$tree <- rep(1:20, each = 50)

## Warts
list.warts.pc2.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.warts.pc2 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_warts_fit_residuals_pc2.RDS")))
    names(mcc.warts.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    warts.pc2 <- result.reorg("warts", pc = 2, nsimmap = 10)
    temp.warts.pc2.av <- vector(mode = "list", length = length(warts.pc2))
    for(j in 1:length(warts.pc2)){
        temp.warts.pc2.av[[j]] <- tryCatch(par.averag(warts.pc2[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.warts.pc2.av <- bind_rows(temp.warts.pc2.av)
    list.warts.pc2.av[[i]] <- cbind(temp.warts.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.warts.pc2.av <- bind_rows(list.warts.pc2.av)
post.warts.pc2.av$tree <- rep(1:20, each = 50)

## Hole Diameter - Discrete
list.holediam.disc.pc2.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.holediam.disc.pc2 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_holediam.disc_fit_residuals_pc2.RDS")))
    names(mcc.holediam.disc.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    holediam.disc.pc2 <- result.reorg("holediam.disc", pc = 2, nsimmap = 10)
    temp.holediam.disc.pc2.av <- vector(mode = "list", length = length(holediam.disc.pc2))
    for(j in 1:length(holediam.disc.pc2)){
        temp.holediam.disc.pc2.av[[j]] <- tryCatch(par.averag(holediam.disc.pc2[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.holediam.disc.pc2.av <- bind_rows(temp.holediam.disc.pc2.av)
    list.holediam.disc.pc2.av[[i]] <- cbind(temp.holediam.disc.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.holediam.disc.pc2.av <- bind_rows(list.holediam.disc.pc2.av)
post.holediam.disc.pc2.av$tree <- rep(1:20, each = 50)


fullpar.table.pc2 <- bind_rows(post.app.pc2.av, post.arch.pc2.av, post.domgrow.pc2.av, post.leafstruc.pc2.av, post.matsys.pc2.av, post.reward.pc2.av, post.strategy.pc2.av, post.warts.pc2.av, post.holediam.disc.pc2.av)
fullpar.table.pc2$disc.trait <- rep(c("app", "arch", "domgrow", "leafstruc", "matsys", "reward", "strategy", "warts", "holediam.disc"), each = 100)

write.table(fullpar.table.pc2, here::here("output/fullpar_table_climpc2_posterior.csv"), sep = ",", quote = FALSE, row.names = FALSE)

beepr::beep("fanfare")

rm(list = ls(pattern = "pc2"))
#rm(fullpar.table.pc2)




#### PC3
## Appendages
list.app.pc3.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.app.pc3 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_appendages_fit_residuals_pc3.RDS")))
    names(mcc.app.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    app.pc3 <- result.reorg("app", pc = 3, nsimmap = 10)
    temp.app.pc3.av <- vector(mode = "list", length = length(app.pc3))
    for(j in 1:length(app.pc3)){
        temp.app.pc3.av[[j]] <- tryCatch(par.averag(app.pc3[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.app.pc3.av <- bind_rows(temp.app.pc3.av)
    list.app.pc3.av[[i]] <- cbind(temp.app.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.app.pc3.av <- bind_rows(list.app.pc3.av)
post.app.pc3.av$tree <- rep(1:20, each = 50)

## Architecture
list.arch.pc3.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.arch.pc3 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_architecture_fit_residuals_pc3.RDS")))
    names(mcc.arch.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    arch.pc3 <- result.reorg("arch", pc = 3, nsimmap = 10)
    temp.arch.pc3.av <- vector(mode = "list", length = length(arch.pc3))
    for(j in 1:length(arch.pc3)){
        temp.arch.pc3.av[[j]] <- tryCatch(par.averag(arch.pc3[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.arch.pc3.av <- bind_rows(temp.arch.pc3.av)
    list.arch.pc3.av[[i]] <- cbind(temp.arch.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.arch.pc3.av <- bind_rows(list.arch.pc3.av)
post.arch.pc3.av$tree <- rep(1:20, each = 50)

## Domatium Growth
list.domgrow.pc3.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.domgrow.pc3 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_dom.growth_fit_residuals_pc3.RDS")))
    names(mcc.domgrow.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    domgrow.pc3 <- result.reorg("domgrow", pc = 3, nsimmap = 10)
    temp.domgrow.pc3.av <- vector(mode = "list", length = length(domgrow.pc3))
    for(j in 1:length(domgrow.pc3)){
        temp.domgrow.pc3.av[[j]] <- tryCatch(par.averag(domgrow.pc3[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.domgrow.pc3.av <- bind_rows(temp.domgrow.pc3.av)
    list.domgrow.pc3.av[[i]] <- cbind(temp.domgrow.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.domgrow.pc3.av <- bind_rows(list.domgrow.pc3.av)
post.domgrow.pc3.av$tree <- rep(1:20, each = 50)

## Leaf Structure
list.leafstruc.pc3.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.leafstruc.pc3 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_leaf.struct_fit_residuals_pc3.RDS")))
    names(mcc.leafstruc.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    leafstruc.pc3 <- result.reorg("leafstruc", pc = 3, nsimmap = 10)
    temp.leafstruc.pc3.av <- vector(mode = "list", length = length(leafstruc.pc3))
    for(j in 1:length(leafstruc.pc3)){
        temp.leafstruc.pc3.av[[j]] <- tryCatch(par.averag(leafstruc.pc3[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.leafstruc.pc3.av <- bind_rows(temp.leafstruc.pc3.av)
    list.leafstruc.pc3.av[[i]] <- cbind(temp.leafstruc.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.leafstruc.pc3.av <- bind_rows(list.leafstruc.pc3.av)
post.leafstruc.pc3.av$tree <- rep(1:20, each = 50)


## Mating System
list.matsys.pc3.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.matsys.pc3 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_mating.system_fit_residuals_pc3.RDS")))
    names(mcc.matsys.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    matsys.pc3 <- result.reorg("matsys", pc = 3, nsimmap = 10)
    temp.matsys.pc3.av <- vector(mode = "list", length = length(matsys.pc3))
    for(j in 1:length(matsys.pc3)){
        temp.matsys.pc3.av[[j]] <- tryCatch(par.averag(matsys.pc3[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.matsys.pc3.av <- bind_rows(temp.matsys.pc3.av)
    list.matsys.pc3.av[[i]] <- cbind(temp.matsys.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.matsys.pc3.av <- bind_rows(list.matsys.pc3.av)
post.matsys.pc3.av$tree <- rep(1:20, each = 50)

## Reward
list.reward.pc3.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.reward.pc3 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_reward_fit_residuals_pc3.RDS")))
    names(mcc.reward.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    reward.pc3 <- result.reorg("reward", pc = 3, nsimmap = 10)
    temp.reward.pc3.av <- vector(mode = "list", length = length(reward.pc3))
    for(j in 1:length(reward.pc3)){
        temp.reward.pc3.av[[j]] <- tryCatch(par.averag(reward.pc3[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.reward.pc3.av <- bind_rows(temp.reward.pc3.av)
    list.reward.pc3.av[[i]] <- cbind(temp.reward.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.reward.pc3.av <- bind_rows(list.reward.pc3.av)
post.reward.pc3.av$tree <- rep(1:20, each = 50)

## Strategy
list.strategy.pc3.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.strategy.pc3 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_strategy_fit_residuals_pc3.RDS")))
    names(mcc.strategy.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    strategy.pc3 <- result.reorg("strategy", pc = 3, nsimmap = 10)
    temp.strategy.pc3.av <- vector(mode = "list", length = length(strategy.pc3))
    for(j in 1:length(strategy.pc3)){
        temp.strategy.pc3.av[[j]] <- tryCatch(par.averag(strategy.pc3[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.strategy.pc3.av <- bind_rows(temp.strategy.pc3.av)
    list.strategy.pc3.av[[i]] <- cbind(temp.strategy.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.strategy.pc3.av <- bind_rows(list.strategy.pc3.av)
post.strategy.pc3.av$tree <- rep(1:20, each = 50)

## Warts
list.warts.pc3.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.warts.pc3 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_warts_fit_residuals_pc3.RDS")))
    names(mcc.warts.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    warts.pc3 <- result.reorg("warts", pc = 3, nsimmap = 10)
    temp.warts.pc3.av <- vector(mode = "list", length = length(warts.pc3))
    for(j in 1:length(warts.pc3)){
        temp.warts.pc3.av[[j]] <- tryCatch(par.averag(warts.pc3[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.warts.pc3.av <- bind_rows(temp.warts.pc3.av)
    list.warts.pc3.av[[i]] <- cbind(temp.warts.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.warts.pc3.av <- bind_rows(list.warts.pc3.av)
post.warts.pc3.av$tree <- rep(1:20, each = 50)

## Hole Diameter - Discrete
list.holediam.disc.pc3.av <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mcc.holediam.disc.pc3 <- readRDS(here::here(paste0("output/posterior_trees/ouwie_mcc_tree", i, "_holediam.disc_fit_residuals_pc3.RDS")))
    names(mcc.holediam.disc.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
    holediam.disc.pc3 <- result.reorg("holediam.disc", pc = 3, nsimmap = 10)
    temp.holediam.disc.pc3.av <- vector(mode = "list", length = length(holediam.disc.pc3))
    for(j in 1:length(holediam.disc.pc3)){
        temp.holediam.disc.pc3.av[[j]] <- tryCatch(par.averag(holediam.disc.pc3[[j]], nsimmap = 10), error = function(y){return(setNames(as.data.frame(matrix(NA, nrow = 10, ncol = 12)), c(paste0("alpha", 1:4), paste0("sigma", 1:4), paste0("theta", 1:4))))})
    }
    temp.holediam.disc.pc3.av <- bind_rows(temp.holediam.disc.pc3.av)
    list.holediam.disc.pc3.av[[i]] <- cbind(temp.holediam.disc.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 10), simmap = rep(1:10, 5)))
}

## Using the replica with the greatest number of states to build the table
post.holediam.disc.pc3.av <- bind_rows(list.holediam.disc.pc3.av)
post.holediam.disc.pc3.av$tree <- rep(1:20, each = 50)


fullpar.table.pc3 <- bind_rows(post.app.pc3.av, post.arch.pc3.av, post.domgrow.pc3.av, post.leafstruc.pc3.av, post.matsys.pc3.av, post.reward.pc3.av, post.strategy.pc3.av, post.warts.pc3.av, post.holediam.disc.pc3.av)
fullpar.table.pc3$disc.trait <- rep(c("app", "arch", "domgrow", "leafstruc", "matsys", "reward", "strategy", "warts", "holediam.disc"), each = 100)

write.table(fullpar.table.pc3, here::here("output/fullpar_table_climpc3_posterior.csv"), sep = ",", quote = FALSE, row.names = FALSE)

beepr::beep("fanfare")

rm(list = ls(pattern = "pc3"))
#rm(fullpar.table.pc3)
