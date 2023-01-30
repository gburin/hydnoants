library("phytools")
library("RColorBrewer")
library("plyr")


load("../output/ouwie_mcc_full_fit.RData")

sum.simmap.app <- summary(simmap.list.appendages)
sum.simmap.arch <- summary(simmap.list.architecture)
sum.simmap.domgrow <- summary(simmap.list.dom.growth)
sum.simmap.leafstr <- summary(simmap.list.leaf.struct)
sum.simmap.matsys <- summary(simmap.list.mating.system)
sum.simmap.reward <- summary(simmap.list.reward)
sum.simmap.strategy <- summary(simmap.list.strategy)
sum.simmap.warts <- summary(simmap.list.warts)

## Organizing into lists
mcc.app <- lapply(ls(pattern = "fit.mcc.appendages"), get)
names(mcc.app) <- gsub("fit.mcc.appendages.", "", ls(pattern = "fit.mcc.appendages"))
mcc.arch <- lapply(ls(pattern = "fit.mcc.architecture"), get)
names(mcc.arch) <- gsub("fit.mcc.architecture.", "", ls(pattern = "fit.mcc.architecture"))
mcc.domgrow <- lapply(ls(pattern = "fit.mcc.dom.growth"), get)
names(mcc.domgrow) <- gsub("fit.mcc.dom.growth.", "", ls(pattern = "fit.mcc.dom.growth"))
mcc.leafstr <- lapply(ls(pattern = "fit.mcc.leaf.struct"), get)
names(mcc.leafstr) <- gsub("fit.mcc.leaf.struct.", "", ls(pattern = "fit.mcc.leaf.struct"))
mcc.matsys <- lapply(ls(pattern = "fit.mcc.mating.system"), get)
names(mcc.matsys) <- gsub("fit.mcc.mating.system.", "", ls(pattern = "fit.mcc.mating.system"))
mcc.reward <- lapply(ls(pattern = "fit.mcc.reward"), get)
names(mcc.reward) <- gsub("fit.mcc.reward.", "", ls(pattern = "fit.mcc.reward"))
mcc.strategy <- lapply(ls(pattern = "fit.mcc.strategy"), get)
names(mcc.strategy) <- gsub("fit.mcc.strategy.", "", ls(pattern = "fit.mcc.strategy"))
mcc.warts <- lapply(ls(pattern = "fit.mcc.warts"), get)
names(mcc.warts) <- gsub("fit.mcc.warts.", "", ls(pattern = "fit.mcc.warts"))


aicc.app <- dplyr::bind_rows(ldply(mcc.app[c("stemarea", "leafarea", "corleng", "petleng", "holediam")], function(z){sapply(z, function(x){sapply(x, "[[", "AICc")})}))
names(aicc.app) <- c("variable", "BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")
aicc.app$best <- apply(aicc.app[, 2:8], 1, function(x){names(x)[which.min(x)]})

aicc.arch <- dplyr::bind_rows(ldply(mcc.arch[c("stemarea", "leafarea", "corleng", "petleng", "holediam")], function(z){sapply(z, function(x){sapply(x, function(y){y$AICc})})}))
names(aicc.arch) <- c("variable", "BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")
aicc.arch$best <- apply(aicc.arch[, 2:8], 1, function(x){names(x)[which.min(x)]})

aicc.domgrow <- dplyr::bind_rows(ldply(mcc.domgrow[c("stemarea", "leafarea", "corleng", "petleng", "holediam")], function(z){sapply(z, function(x){sapply(x, function(y){y$AICc})})}))
names(aicc.domgrow) <- c("variable", "BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")
aicc.domgrow$best <- apply(aicc.domgrow[, 2:8], 1, function(x){names(x)[which.min(x)]})

aicc.leafstr <- dplyr::bind_rows(ldply(mcc.leafstr[c("stemarea", "leafarea", "corleng", "petleng", "holediam")], function(z){sapply(z, function(x){sapply(x, function(y){y$AICc})})}))
names(aicc.leafstr) <- c("variable", "BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")
aicc.leafstr$best <- apply(aicc.leafstr[, 2:8], 1, function(x){names(x)[which.min(x)]})

aicc.matsys <- dplyr::bind_rows(ldply(mcc.matsys[c("stemarea", "leafarea", "corleng", "petleng", "holediam")], function(z){sapply(z, function(x){sapply(x, function(y){y$AICc})})}))
names(aicc.matsys) <- c("variable", "BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")
aicc.matsys$best <- apply(aicc.matsys[, 2:8], 1, function(x){names(x)[which.min(x)]})

aicc.reward <- dplyr::bind_rows(ldply(mcc.reward[c("stemarea", "leafarea", "corleng", "petleng", "holediam")], function(z){sapply(z, function(x){sapply(x, function(y){y$AICc})})}))
names(aicc.reward) <- c("variable", "BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")
aicc.reward$best <- apply(aicc.reward[, 2:8], 1, function(x){names(x)[which.min(x)]})

aicc.strategy <- dplyr::bind_rows(ldply(mcc.strategy[c("stemarea", "leafarea", "corleng", "petleng", "holediam")], function(z){sapply(z, function(x){sapply(x, function(y){y$AICc})})}))
names(aicc.strategy) <- c("variable", "BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")
aicc.strategy$best <- apply(aicc.strategy[, 2:8], 1, function(x){names(x)[which.min(x)]})

aicc.warts <- dplyr::bind_rows(ldply(mcc.warts[c("stemarea", "leafarea", "corleng", "petleng", "holediam")], function(z){sapply(z, function(x){sapply(x, function(y){y$AICc})})}))
names(aicc.warts) <- c("variable", "BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")
aicc.warts$best <- apply(aicc.warts[, 2:8], 1, function(x){names(x)[which.min(x)]})

fullaicc.table <- bind_rows(aicc.leafstr, aicc.app, aicc.arch, aicc.domgrow, aicc.matsys, aicc.reward, aicc.strategy, aicc.warts)
fullaicc.table$trait <- rep(c("Leaf.Structure", "Appendages", "Architecture", "Domatium.Growth", "Mating.System", "Reward", "Strategy", "Warts"), each = 500)

write.table(fullaicc.table, "../output/fullaicc_table.csv", sep = ",", quote = FALSE, row.names = FALSE)



### Extracting parameter values

aicc.app.pars <- ldply(1:5, function(y){ldply(1:length(mcc.app[[y]]), function(z){ldply(1:length(mcc.app[[y]][[z]]), function(x){setNames(c(mcc.app[[y]][[z]][[x]]$solution[1,], mcc.app[[y]][[z]][[x]]$solution[2,], mcc.app[[y]][[z]][[x]]$theta[,1]), c(paste0("alpha", 1:ncol(mcc.app[[y]][[z]][[x]]$solution)), paste0("sigma", 1:ncol(mcc.app[[y]][[z]][[x]]$solution)), paste0("theta", 1:ncol(mcc.app[[y]][[z]][[x]]$solution))))})})})
aicc.app.pars$replica <- rep(rep(1:100, times = 7), 5)
aicc.app.pars$model <- rep(rep(names(aicc.app)[2:8], each = 100), 5)
aicc.app.pars$variable <- rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 700)

aicc.arch.pars <- ldply(1:5, function(y){ldply(1:length(mcc.arch[[y]]), function(z){ldply(1:length(mcc.arch[[y]][[z]]), function(x){setNames(c(mcc.arch[[y]][[z]][[x]]$solution[1,], mcc.arch[[y]][[z]][[x]]$solution[2,], mcc.arch[[y]][[z]][[x]]$theta[,1]), c(paste0("alpha", 1:ncol(mcc.arch[[y]][[z]][[x]]$solution)), paste0("sigma", 1:ncol(mcc.arch[[y]][[z]][[x]]$solution)), paste0("theta", 1:ncol(mcc.arch[[y]][[z]][[x]]$solution))))})})})
aicc.arch.pars$replica <- rep(rep(1:100, times = 7), 5)
aicc.arch.pars$model <- rep(rep(names(aicc.arch)[2:8], each = 100), 5)
aicc.arch.pars$variable <- rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 700)

aicc.domgrow.pars <- ldply(1:5, function(y){ldply(1:length(mcc.domgrow[[y]]), function(z){ldply(1:length(mcc.domgrow[[y]][[z]]), function(x){setNames(c(mcc.domgrow[[y]][[z]][[x]]$solution[1,], mcc.domgrow[[y]][[z]][[x]]$solution[2,], mcc.domgrow[[y]][[z]][[x]]$theta[,1]), c(paste0("alpha", 1:ncol(mcc.domgrow[[y]][[z]][[x]]$solution)), paste0("sigma", 1:ncol(mcc.domgrow[[y]][[z]][[x]]$solution)), paste0("theta", 1:ncol(mcc.domgrow[[y]][[z]][[x]]$solution))))})})})
aicc.domgrow.pars$replica <- rep(rep(1:100, times = 7), 5)
aicc.domgrow.pars$model <- rep(rep(names(aicc.domgrow)[2:8], each = 100), 5)
aicc.domgrow.pars$variable <- rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 700)

aicc.leafstr.pars <- ldply(1:5, function(y){ldply(1:length(mcc.leafstr[[y]]), function(z){ldply(1:length(mcc.leafstr[[y]][[z]]), function(x){setNames(c(mcc.leafstr[[y]][[z]][[x]]$solution[1,], mcc.leafstr[[y]][[z]][[x]]$solution[2,], mcc.leafstr[[y]][[z]][[x]]$theta[,1]), c(paste0("alpha", 1:ncol(mcc.leafstr[[y]][[z]][[x]]$solution)), paste0("sigma", 1:ncol(mcc.leafstr[[y]][[z]][[x]]$solution)), paste0("theta", 1:ncol(mcc.leafstr[[y]][[z]][[x]]$solution))))})})})
aicc.leafstr.pars$replica <- rep(rep(1:100, times = 7), 5)
aicc.leafstr.pars$model <- rep(rep(names(aicc.leafstr)[2:8], each = 100), 5)
aicc.leafstr.pars$variable <- rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 700)

aicc.matsys.pars <- ldply(1:5, function(y){ldply(1:length(mcc.matsys[[y]]), function(z){ldply(1:length(mcc.matsys[[y]][[z]]), function(x){setNames(c(mcc.matsys[[y]][[z]][[x]]$solution[1,], mcc.matsys[[y]][[z]][[x]]$solution[2,], mcc.matsys[[y]][[z]][[x]]$theta[,1]), c(paste0("alpha", 1:ncol(mcc.matsys[[y]][[z]][[x]]$solution)), paste0("sigma", 1:ncol(mcc.matsys[[y]][[z]][[x]]$solution)), paste0("theta", 1:ncol(mcc.matsys[[y]][[z]][[x]]$solution))))})})})
aicc.matsys.pars$replica <- rep(rep(1:100, times = 7), 5)
aicc.matsys.pars$model <- rep(rep(names(aicc.matsys)[2:8], each = 100), 5)
aicc.matsys.pars$variable <- rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 700)

aicc.reward.pars <- ldply(1:5, function(y){ldply(1:length(mcc.reward[[y]]), function(z){ldply(1:length(mcc.reward[[y]][[z]]), function(x){setNames(c(mcc.reward[[y]][[z]][[x]]$solution[1,], mcc.reward[[y]][[z]][[x]]$solution[2,], mcc.reward[[y]][[z]][[x]]$theta[,1]), c(paste0("alpha", 1:ncol(mcc.reward[[y]][[z]][[x]]$solution)), paste0("sigma", 1:ncol(mcc.reward[[y]][[z]][[x]]$solution)), paste0("theta", 1:ncol(mcc.reward[[y]][[z]][[x]]$solution))))})})})
aicc.reward.pars$replica <- rep(rep(1:100, times = 7), 5)
aicc.reward.pars$model <- rep(rep(names(aicc.reward)[2:8], each = 100), 5)
aicc.reward.pars$variable <- rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 700)

aicc.strategy.pars <- ldply(1:5, function(y){ldply(1:length(mcc.strategy[[y]]), function(z){ldply(1:length(mcc.strategy[[y]][[z]]), function(x){setNames(c(mcc.strategy[[y]][[z]][[x]]$solution[1,], mcc.strategy[[y]][[z]][[x]]$solution[2,], mcc.strategy[[y]][[z]][[x]]$theta[,1]), c(paste0("alpha", 1:ncol(mcc.strategy[[y]][[z]][[x]]$solution)), paste0("sigma", 1:ncol(mcc.strategy[[y]][[z]][[x]]$solution)), paste0("theta", 1:ncol(mcc.strategy[[y]][[z]][[x]]$solution))))})})})
aicc.strategy.pars$replica <- rep(rep(1:100, times = 7), 5)
aicc.strategy.pars$model <- rep(rep(names(aicc.strategy)[2:8], each = 100), 5)
aicc.strategy.pars$variable <- rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 700)

aicc.warts.pars <- ldply(1:5, function(y){ldply(1:length(mcc.warts[[y]]), function(z){ldply(1:length(mcc.warts[[y]][[z]]), function(x){setNames(c(mcc.warts[[y]][[z]][[x]]$solution[1,], mcc.warts[[y]][[z]][[x]]$solution[2,], mcc.warts[[y]][[z]][[x]]$theta[,1]), c(paste0("alpha", 1:ncol(mcc.warts[[y]][[z]][[x]]$solution)), paste0("sigma", 1:ncol(mcc.warts[[y]][[z]][[x]]$solution)), paste0("theta", 1:ncol(mcc.warts[[y]][[z]][[x]]$solution))))})})})
aicc.warts.pars$replica <- rep(rep(1:100, times = 7), 5)
aicc.warts.pars$model <- rep(rep(names(aicc.warts)[2:8], each = 100), 5)
aicc.warts.pars$variable <- rep(c("stemarea", "leafarea", "corleng", "petleng", "holediam"), each = 700)


fullpar.table <- bind_rows(aicc.leafstr.pars, aicc.app.pars, aicc.arch.pars, aicc.domgrow.pars, aicc.matsys.pars, aicc.reward.pars, aicc.strategy.pars, aicc.warts.pars)
fullpar.table$trait <- rep(c("Leaf.Structure", "Appendages", "Architecture", "Domatium.Growth", "Mating.System", "Reward", "Strategy", "Warts"), each = 3500)

write.table(fullpar.table, "../output/fullpar_table.csv", sep = ",", quote = FALSE, row.names = FALSE)
