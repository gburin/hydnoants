library("phytools")
library("RColorBrewer")
library("plyr")
library("patchwork")
library("tidyverse")

here::i_am("R/ouwie_summary_mcc_climpc.R")

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

result.reorg <- function(x, pc){
    traits <- c("stemarea", "leafarea", "corleng", "petleng", "holediam")
    data <- get(paste0("mcc.", x, ".pc", pc))
    res <- vector(mode = "list", length = 100)
    for(k in 1:length(traits)){
        for(i in 1:100){
            res[[i]] <- vector(mode = "list", length = 6)
            for(j in 1:6){
                res[[i]][[j]] <- data[[(k + 2)]][[j]][[i]]
            }
        }
        assign(paste0(traits[k], ".res"), res)
    }
    return(list(stemarea.res, leafarea.res, corleng.res, petleng.res, holediam.res))
}

par.averag <- function(res){
    pars <- llply(res, function(y){suppressMessages(bind_rows(llply(y, par.extr)))})
    aicc <- ldply(res, function(y){t(ldply(y, "[[", "AICc"))})
    daicc <- apply(aicc, 1, function(x){x - min(x)})
    aicc.w <- exp(-0.5 * daicc)
    mod.av <- as.data.frame(matrix(NA, ncol = 12, nrow = 100))
    names(mod.av) <- paste0(rep(c("alpha", "sigma", "theta"), each = 4), 1:4)
    for(i in 1:100){
        pars[[i]]$aic.w <- aicc.w[, i]/sum(aicc.w[, i])
        temp <- suppressMessages(bind_cols(llply(pars[[i]][, 1:(ncol(pars[[i]]) - 1)], function(x){x * pars[[i]][, ncol(pars[[i]])]})))
        mod.av[i, match(names(pars[[i]])[-ncol(pars[[i]])], names(mod.av))] <- apply(temp, 2, sum, na.rm = TRUE)
    }
    return(mod.av)
}


#### PC1

mcc.app.pc1 <- readRDS(here::here("output/ouwie_mcc_appendages_fit_residuals_pc1.RDS"))
names(mcc.app.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
app.pc1 <- result.reorg("app", pc = 1)
app.pc1.av <- bind_rows(lapply(app.pc1, par.averag))
app.pc1.av <- cbind(app.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.arch.pc1 <- readRDS(here::here("output/ouwie_mcc_architecture_fit_residuals_pc1.RDS"))
names(mcc.arch.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
arch.pc1 <- result.reorg("arch", pc = 1)
arch.pc1.av <- bind_rows(lapply(arch.pc1, par.averag))
arch.pc1.av <- cbind(arch.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.domgrow.pc1 <- readRDS(here::here("output/ouwie_mcc_dom.growth_fit_residuals_pc1.RDS"))
names(mcc.domgrow.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
domgrow.pc1 <- result.reorg("domgrow", pc = 1)
domgrow.pc1.av <- bind_rows(lapply(domgrow.pc1, par.averag))
domgrow.pc1.av <- cbind(domgrow.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.leafstruc.pc1 <- readRDS(here::here("output/ouwie_mcc_leaf.struct_fit_residuals_pc1.RDS"))
names(mcc.leafstruc.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
leafstruc.pc1 <- result.reorg("leafstruc", pc = 1)
leafstruc.pc1.av <- bind_rows(lapply(leafstruc.pc1, par.averag))
leafstruc.pc1.av <- cbind(leafstruc.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.matsys.pc1 <- readRDS(here::here("output/ouwie_mcc_mating.system_fit_residuals_pc1.RDS"))
names(mcc.matsys.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
matsys.pc1 <- result.reorg("matsys", pc = 1)
matsys.pc1.av <- bind_rows(lapply(matsys.pc1, par.averag))
matsys.pc1.av <- cbind(matsys.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.reward.pc1 <- readRDS(here::here("output/ouwie_mcc_reward_fit_residuals_pc1.RDS"))
names(mcc.reward.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
reward.pc1 <- result.reorg("reward", pc = 1)
reward.pc1.av <- bind_rows(lapply(reward.pc1, par.averag))
reward.pc1.av <- cbind(reward.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.strategy.pc1 <- readRDS(here::here("output/ouwie_mcc_strategy_fit_residuals_pc1.RDS"))
names(mcc.strategy.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
strategy.pc1 <- result.reorg("strategy", pc = 1)
strategy.pc1.av <- bind_rows(lapply(strategy.pc1, par.averag))
strategy.pc1.av <- cbind(strategy.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.warts.pc1 <- readRDS(here::here("output/ouwie_mcc_warts_fit_residuals_pc1.RDS"))
names(mcc.warts.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
warts.pc1 <- result.reorg("warts", pc = 1)
warts.pc1.av <- bind_rows(lapply(warts.pc1, par.averag))
warts.pc1.av <- cbind(warts.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.holediam.disc.pc1 <- readRDS(here::here("output/ouwie_mcc_holediam.disc_fit_residuals_pc1.RDS"))
names(mcc.holediam.disc.pc1) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
holediam.disc.pc1 <- result.reorg("holediam.disc", pc = 1)
holediam.disc.pc1.av <- bind_rows(lapply(holediam.disc.pc1, par.averag))
holediam.disc.pc1.av <- cbind(holediam.disc.pc1.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))





fullpar.table.pc1 <- bind_rows(app.pc1.av, arch.pc1.av, domgrow.pc1.av, leafstruc.pc1.av, matsys.pc1.av, reward.pc1.av, strategy.pc1.av, warts.pc1.av, holediam.disc.pc1.av)
fullpar.table.pc1$disc.trait <- rep(c("app", "arch", "domgrow", "leafstruc", "matsys", "reward", "strategy", "warts", "holediam.disc"), each = 500)
fullpar.table.pc1 <- fullpar.table.pc1[, c(13:15,
                                   grep("alpha", names(fullpar.table.pc1)),
                                   grep("sigma", names(fullpar.table.pc1)),
                                   grep("theta", names(fullpar.table.pc1)))]

beepr::beep("fanfare")

write.table(fullpar.table.pc1, here::here("output/fullpar_table_climpc1_new.csv"), sep = ",", quote = FALSE, row.names = FALSE)

rm(list = ls(pattern = "pc1"))
#rm(fullpar.table.pc1)
   
#### PC2

## Organizing into lists
mcc.app.pc2 <- readRDS(here::here("output/ouwie_mcc_appendages_fit_residuals_pc2.RDS"))
names(mcc.app.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
app.pc2 <- result.reorg("app", pc = 2)
app.pc2.av <- bind_rows(lapply(app.pc2, par.averag))
app.pc2.av <- cbind(app.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.arch.pc2 <- readRDS(here::here("output/ouwie_mcc_architecture_fit_residuals_pc2.RDS"))
names(mcc.arch.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
arch.pc2 <- result.reorg("arch", pc = 2)
arch.pc2.av <- bind_rows(lapply(arch.pc2, par.averag))
arch.pc2.av <- cbind(arch.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.domgrow.pc2 <- readRDS(here::here("output/ouwie_mcc_dom.growth_fit_residuals_pc2.RDS"))
names(mcc.domgrow.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
domgrow.pc2 <- result.reorg("domgrow", pc = 2)
domgrow.pc2.av <- bind_rows(lapply(domgrow.pc2, par.averag))
domgrow.pc2.av <- cbind(domgrow.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.leafstruc.pc2 <- readRDS(here::here("output/ouwie_mcc_leaf.struct_fit_residuals_pc2.RDS"))
names(mcc.leafstruc.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
leafstruc.pc2 <- result.reorg("leafstruc", pc = 2)
leafstruc.pc2.av <- bind_rows(lapply(leafstruc.pc2, par.averag))
leafstruc.pc2.av <- cbind(leafstruc.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.matsys.pc2 <- readRDS(here::here("output/ouwie_mcc_mating.system_fit_residuals_pc2.RDS"))
names(mcc.matsys.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
matsys.pc2 <- result.reorg("matsys", pc = 2)
matsys.pc2.av <- bind_rows(lapply(matsys.pc2, par.averag))
matsys.pc2.av <- cbind(matsys.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.reward.pc2 <- readRDS(here::here("output/ouwie_mcc_reward_fit_residuals_pc2.RDS"))
names(mcc.reward.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
reward.pc2 <- result.reorg("reward", pc = 2)
reward.pc2.av <- bind_rows(lapply(reward.pc2, par.averag))
reward.pc2.av <- cbind(reward.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.strategy.pc2 <- readRDS(here::here("output/ouwie_mcc_strategy_fit_residuals_pc2.RDS"))
names(mcc.strategy.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
strategy.pc2 <- result.reorg("strategy", pc = 2)
strategy.pc2.av <- bind_rows(lapply(strategy.pc2, par.averag))
strategy.pc2.av <- cbind(strategy.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.warts.pc2 <- readRDS(here::here("output/ouwie_mcc_warts_fit_residuals_pc2.RDS"))
names(mcc.warts.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
warts.pc2 <- result.reorg("warts", pc = 2)
warts.pc2.av <- bind_rows(lapply(warts.pc2, par.averag))
warts.pc2.av <- cbind(warts.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.holediam.disc.pc2 <- readRDS(here::here("output/ouwie_mcc_holediam.disc_fit_residuals_pc2.RDS"))
names(mcc.holediam.disc.pc2) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
holediam.disc.pc2 <- result.reorg("holediam.disc", pc = 2)
holediam.disc.pc2.av <- bind_rows(lapply(holediam.disc.pc2, par.averag))
holediam.disc.pc2.av <- cbind(holediam.disc.pc2.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

fullpar.table.pc2 <- bind_rows(app.pc2.av, arch.pc2.av, domgrow.pc2.av, leafstruc.pc2.av, matsys.pc2.av, reward.pc2.av, strategy.pc2.av, warts.pc2.av, holediam.disc.pc2.av)
fullpar.table.pc2$disc.trait <- rep(c("app", "arch", "domgrow", "leafstruc", "matsys", "reward", "strategy", "warts", "holediam.disc"), each = 500)
fullpar.table.pc2 <- fullpar.table.pc2[, c(13:15,
                                   grep("alpha", names(fullpar.table.pc2)),
                                   grep("sigma", names(fullpar.table.pc2)),
                                   grep("theta", names(fullpar.table.pc2)))]

beepr::beep("fanfare")

write.table(fullpar.table.pc2, here::here("output/fullpar_table_climpc2_new.csv"), sep = ",", quote = FALSE, row.names = FALSE)

rm(list = ls(pattern = "pc2"))
#rm(fullpar.table)


#### PC3

## Organizing into lists
mcc.app.pc3 <- readRDS(here::here("output/ouwie_mcc_appendages_fit_residuals_pc3.RDS"))
names(mcc.app.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
app.pc3 <- result.reorg("app", pc = 3)
app.pc3.av <- bind_rows(lapply(app.pc3, par.averag))
app.pc3.av <- cbind(app.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.arch.pc3 <- readRDS(here::here("output/ouwie_mcc_architecture_fit_residuals_pc3.RDS"))
names(mcc.arch.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
arch.pc3 <- result.reorg("arch", pc = 3)
arch.pc3.av <- bind_rows(lapply(arch.pc3, par.averag))
arch.pc3.av <- cbind(arch.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.domgrow.pc3 <- readRDS(here::here("output/ouwie_mcc_dom.growth_fit_residuals_pc3.RDS"))
names(mcc.domgrow.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
domgrow.pc3 <- result.reorg("domgrow", pc = 3)
domgrow.pc3.av <- bind_rows(lapply(domgrow.pc3, par.averag))
domgrow.pc3.av <- cbind(domgrow.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.leafstruc.pc3 <- readRDS(here::here("output/ouwie_mcc_leaf.struct_fit_residuals_pc3.RDS"))
names(mcc.leafstruc.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
leafstruc.pc3 <- result.reorg("leafstruc", pc = 3)
leafstruc.pc3.av <- bind_rows(lapply(leafstruc.pc3, par.averag))
leafstruc.pc3.av <- cbind(leafstruc.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.matsys.pc3 <- readRDS(here::here("output/ouwie_mcc_mating.system_fit_residuals_pc3.RDS"))
names(mcc.matsys.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
matsys.pc3 <- result.reorg("matsys", pc = 3)
matsys.pc3.av <- bind_rows(lapply(matsys.pc3, par.averag))
matsys.pc3.av <- cbind(matsys.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.reward.pc3 <- readRDS(here::here("output/ouwie_mcc_reward_fit_residuals_pc3.RDS"))
names(mcc.reward.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
reward.pc3 <- result.reorg("reward", pc = 3)
reward.pc3.av <- bind_rows(lapply(reward.pc3, par.averag))
reward.pc3.av <- cbind(reward.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.strategy.pc3 <- readRDS(here::here("output/ouwie_mcc_strategy_fit_residuals_pc3.RDS"))
names(mcc.strategy.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
strategy.pc3 <- result.reorg("strategy", pc = 3)
strategy.pc3.av <- bind_rows(lapply(strategy.pc3, par.averag))
strategy.pc3.av <- cbind(strategy.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.warts.pc3 <- readRDS(here::here("output/ouwie_mcc_warts_fit_residuals_pc3.RDS"))
names(mcc.warts.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
warts.pc3 <- result.reorg("warts", pc = 3)
warts.pc3.av <- bind_rows(lapply(warts.pc3, par.averag))
warts.pc3.av <- cbind(warts.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

mcc.holediam.disc.pc3 <- readRDS(here::here("output/ouwie_mcc_holediam.disc_fit_residuals_pc3.RDS"))
names(mcc.holediam.disc.pc3) <- c("simmap.sample", "simmap.list", "stemarea", "leafarea", "corleng", "petleng", "holediam")
holediam.disc.pc3 <- result.reorg("holediam.disc", pc = 3)
holediam.disc.pc3.av <- bind_rows(lapply(holediam.disc.pc3, par.averag))
holediam.disc.pc3.av <- cbind(holediam.disc.pc3.av, data.frame(trait = rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 100), simmap = rep(1:100, 5)))

fullpar.table.pc3 <- bind_rows(app.pc3.av, arch.pc3.av, domgrow.pc3.av, leafstruc.pc3.av, matsys.pc3.av, reward.pc3.av, strategy.pc3.av, warts.pc3.av, holediam.disc.pc3.av)
fullpar.table.pc3$disc.trait <- rep(c("app", "arch", "domgrow", "leafstruc", "matsys", "reward", "strategy", "warts", "holediam.disc"), each = 500)
fullpar.table.pc3 <- fullpar.table.pc3[, c(13:15,
                                   grep("alpha", names(fullpar.table.pc3)),
                                   grep("sigma", names(fullpar.table.pc3)),
                                   grep("theta", names(fullpar.table.pc3)))]

beepr::beep("fanfare")

write.table(fullpar.table.pc3, here::here("output/fullpar_table_climpc3_new.csv"), sep = ",", quote = FALSE, row.names = FALSE)
