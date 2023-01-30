here::i_am("R/data_prep_ouwie_climpc_clean.R")

#source(here::here("R/tree_plots.R"))
library("knitr")
library("RColorBrewer")
library("viridis")
library("kableExtra")
## load("../output/first_analysis.RData")
knitr::opts_chunk$set(fig.pos = "h", out.extra = "")
## source("../R/ouwie_summary_mcc.R")
pars.table.pc1 <- read.csv(here::here("output/fullpar_table_climpc1.csv"))
pars.table.pc1$trait <- c("Stem Area", "Leaf Area", "Corola Length", "Petiole Length", "Hole Diameter")[match(pars.table.pc1$trait, c("stemarea", "leafarea", "corleng", "petleng", "holediam"))]
pars.table.pc2 <- read.csv(here::here("output/fullpar_table_climpc2.csv"))
pars.table.pc2$trait <- c("Stem Area", "Leaf Area", "Corola Length", "Petiole Length", "Hole Diameter")[match(pars.table.pc2$trait, c("stemarea", "leafarea", "corleng", "petleng", "holediam"))]
pars.table.pc3 <- read.csv(here::here("output/fullpar_table_climpc3.csv"))
pars.table.pc3$trait <- c("Stem Area", "Leaf Area", "Corola Length", "Petiole Length", "Hole Diameter")[match(pars.table.pc3$trait, c("stemarea", "leafarea", "corleng", "petleng", "holediam"))]
colors <- setNames(brewer.pal(4, "Set1"), 1:4)

states.app <- c(#"outgroup",
                "none",
                "variable",
                "spines")
states.arch <- c(#"shrub",
                 "multiple",
                 "single")
states.domgrow <- c(#"outgroup",
                    "diffuse",
                    "apical")
states.leafstruc <- c("thick",
                      "var",
                      "thin",
                      "succ")
states.matsys <- c("heterostyl",
                   "nonHeterostyl",
                   "functUnisex")
states.reward <- c(#"outgroup",
                   "absent",
                   "present")
states.strategy <- c(#"outgroup",
                     "facultative",
                     "obligate",
                     "lost")
states.warts <- c(#"outgroup",
                  "variable",
                  "diff",
                  "lost")

## data.plot.pc1 <- reshape2::melt(pars.table.pc1, id.vars = c("trait", "disc.trait", "simmap"))
## theta.plot.pc1 <- data.plot.pc1[grep("theta", data.plot.pc1$variable),]
## alpha.plot.pc1 <- data.plot.pc1[grep("alpha", data.plot.pc1$variable),]
## sigma.plot.pc1 <- data.plot.pc1[grep("sigma", data.plot.pc1$variable),]

## data.plot.pc2 <- reshape2::melt(pars.table.pc2, id.vars = c("trait", "disc.trait", "simmap"))
## theta.plot.pc2 <- data.plot.pc2[grep("theta", data.plot.pc2$variable),]
## alpha.plot.pc2 <- data.plot.pc2[grep("alpha", data.plot.pc2$variable),]
## sigma.plot.pc2 <- data.plot.pc2[grep("sigma", data.plot.pc2$variable),]

## data.plot.pc3 <- reshape2::melt(pars.table.pc3, id.vars = c("trait", "disc.trait", "simmap"))
## theta.plot.pc3 <- data.plot.pc3[grep("theta", data.plot.pc3$variable),]
## alpha.plot.pc3 <- data.plot.pc3[grep("alpha", data.plot.pc3$variable),]
## sigma.plot.pc3 <- data.plot.pc3[grep("sigma", data.plot.pc3$variable),]


## Appendages
### PC1

theta.app.plot.pc1 <- na.omit(subset(theta.plot.pc1, disc.trait == "app"))
theta.app.plot.pc1$variable <- as.character(theta.app.plot.pc1$variable)
theta.app.plot.pc1 <- theta.app.plot.pc1[theta.app.plot.pc1$value >= -10 & theta.app.plot.pc1$value <= 10, ]

alpha.app.plot.pc1 <- na.omit(subset(alpha.plot.pc1, disc.trait == "app"))
alpha.app.plot.pc1$variable <- as.character(alpha.app.plot.pc1$variable)
alpha.app.plot.pc1 <- alpha.app.plot.pc1[alpha.app.plot.pc1$value >= -10 & alpha.app.plot.pc1$value <= 10, ]

sigma.app.plot.pc1 <- na.omit(subset(sigma.plot.pc1, disc.trait == "app"))
sigma.app.plot.pc1$variable <- as.character(sigma.app.plot.pc1$variable)
sigma.app.plot.pc1 <- sigma.app.plot.pc1[sigma.app.plot.pc1$value >= -10 & sigma.app.plot.pc1$value <= 10, ]


#### Differences
app.plot.diff.pc1 <- pars.table.pc1[pars.table.pc1$disc.trait == "app", -c(6, 10, 14, 15)]

for(i in 1:length(states.app)){
    for(j in 1:length(states.app)){
        app.plot.diff.pc1[, paste0("theta", i, j)] <- app.plot.diff.pc1[, paste0("theta", i)] > app.plot.diff.pc1[, paste0("theta", j)]
        app.plot.diff.pc1[, paste0("alpha", i, j)] <- app.plot.diff.pc1[, paste0("alpha", i)] > app.plot.diff.pc1[, paste0("alpha", j)]
        app.plot.diff.pc1[, paste0("sigma", i, j)] <- app.plot.diff.pc1[, paste0("sigma", i)] > app.plot.diff.pc1[, paste0("sigma", j)]
    }
}

agg.app.diff.pc1 <- aggregate(app.plot.diff.pc1[app.plot.diff.pc1$trait == "Hole Diameter", 12:38], by = list(app.plot.diff.pc1$trait[app.plot.diff.pc1$trait == "Hole Diameter"], app.plot.diff.pc1$disc.trait[app.plot.diff.pc1$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.app.diff.pc1$Group.1)){
    assign(paste0("app.diff.theta.pc1.", gsub(" ", ".", agg.app.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.app.diff.pc1[i, grep("theta", names(agg.app.diff.pc1))])), ncol = length(states.app), nrow = length(states.app), byrow = TRUE), row.names = states.app))
    assign(paste0("app.diff.alpha.pc1.", gsub(" ", ".", agg.app.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.app.diff.pc1[i, grep("alpha", names(agg.app.diff.pc1))])), ncol = length(states.app), nrow = length(states.app), byrow = TRUE), row.names = states.app))
    assign(paste0("app.diff.sigma.pc1.", gsub(" ", ".", agg.app.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.app.diff.pc1[i, grep("sigma", names(agg.app.diff.pc1))])), ncol = length(states.app), nrow = length(states.app), byrow = TRUE), row.names = states.app))
}



### PC2

theta.app.plot.pc2 <- subset(theta.plot.pc2, disc.trait == "app")
theta.app.plot.pc2$variable <- as.character(theta.app.plot.pc2$variable)
theta.app.plot.pc2 <- theta.app.plot.pc2[theta.app.plot.pc2$value >= -10 & theta.app.plot.pc2$value <= 10, ]

alpha.app.plot.pc2 <- subset(alpha.plot.pc2, disc.trait == "app")
alpha.app.plot.pc2$variable <- as.character(alpha.app.plot.pc2$variable)
alpha.app.plot.pc2 <- alpha.app.plot.pc2[alpha.app.plot.pc2$value >= -10 & alpha.app.plot.pc2$value <= 10, ]

sigma.app.plot.pc2 <- subset(sigma.plot.pc2, disc.trait == "app")
sigma.app.plot.pc2$variable <- as.character(sigma.app.plot.pc2$variable)
sigma.app.plot.pc2 <- sigma.app.plot.pc2[sigma.app.plot.pc2$value >= -10 & sigma.app.plot.pc2$value <= 10, ]

app.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "app", -c(6, 10, 14, 15)]
app.diff.pc2 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(app.plot.diff.pc2[, grep("theta", names(app.plot.diff.pc2))] >= app.plot.diff.pc2$theta1))


#### Differences
app.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "app", -c(6, 10, 14, 15)]

for(i in 1:length(states.app)){
    for(j in 1:length(states.app)){
        app.plot.diff.pc2[, paste0("theta", i, j)] <- app.plot.diff.pc2[, paste0("theta", i)] > app.plot.diff.pc2[, paste0("theta", j)]
        app.plot.diff.pc2[, paste0("alpha", i, j)] <- app.plot.diff.pc2[, paste0("alpha", i)] > app.plot.diff.pc2[, paste0("alpha", j)]
        app.plot.diff.pc2[, paste0("sigma", i, j)] <- app.plot.diff.pc2[, paste0("sigma", i)] > app.plot.diff.pc2[, paste0("sigma", j)]
    }
}

agg.app.diff.pc2 <- aggregate(app.plot.diff.pc2[app.plot.diff.pc2$trait == "Hole Diameter", 12:38], by = list(app.plot.diff.pc2$trait[app.plot.diff.pc2$trait == "Hole Diameter"], app.plot.diff.pc2$disc.trait[app.plot.diff.pc2$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.app.diff.pc2$Group.1)){
    assign(paste0("app.diff.theta.pc2.", gsub(" ", ".", agg.app.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.app.diff.pc2[i, grep("theta", names(agg.app.diff.pc2))])), ncol = length(states.app), nrow = length(states.app), byrow = TRUE), row.names = states.app))
    assign(paste0("app.diff.alpha.pc2.", gsub(" ", ".", agg.app.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.app.diff.pc2[i, grep("alpha", names(agg.app.diff.pc2))])), ncol = length(states.app), nrow = length(states.app), byrow = TRUE), row.names = states.app))
    assign(paste0("app.diff.sigma.pc2.", gsub(" ", ".", agg.app.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.app.diff.pc2[i, grep("sigma", names(agg.app.diff.pc2))])), ncol = length(states.app), nrow = length(states.app), byrow = TRUE), row.names = states.app))
}


### PC3

theta.app.plot.pc3 <- subset(theta.plot.pc3, disc.trait == "app")
theta.app.plot.pc3$variable <- as.character(theta.app.plot.pc3$variable)
theta.app.plot.pc3 <- theta.app.plot.pc3[theta.app.plot.pc3$value >= -10 & theta.app.plot.pc3$value <= 10, ]

alpha.app.plot.pc3 <- subset(alpha.plot.pc3, disc.trait == "app")
alpha.app.plot.pc3$variable <- as.character(alpha.app.plot.pc3$variable)
alpha.app.plot.pc3 <- alpha.app.plot.pc3[alpha.app.plot.pc3$value >= -10 & alpha.app.plot.pc3$value <= 10, ]

sigma.app.plot.pc3 <- subset(sigma.plot.pc3, disc.trait == "app")
sigma.app.plot.pc3$variable <- as.character(sigma.app.plot.pc3$variable)
sigma.app.plot.pc3 <- sigma.app.plot.pc3[sigma.app.plot.pc3$value >= -10 & sigma.app.plot.pc3$value <= 10, ]

app.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "app", -c(6, 10, 14, 15)]
app.diff.pc3 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(app.plot.diff.pc3[, grep("theta", names(app.plot.diff.pc3))] >= app.plot.diff.pc3$theta1))


#### Differences
app.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "app", -c(6, 10, 14, 15)]

for(i in 1:length(states.app)){
    for(j in 1:length(states.app)){
        app.plot.diff.pc3[, paste0("theta", i, j)] <- app.plot.diff.pc3[, paste0("theta", i)] > app.plot.diff.pc3[, paste0("theta", j)]
        app.plot.diff.pc3[, paste0("alpha", i, j)] <- app.plot.diff.pc3[, paste0("alpha", i)] > app.plot.diff.pc3[, paste0("alpha", j)]
        app.plot.diff.pc3[, paste0("sigma", i, j)] <- app.plot.diff.pc3[, paste0("sigma", i)] > app.plot.diff.pc3[, paste0("sigma", j)]
    }
}

agg.app.diff.pc3 <- aggregate(app.plot.diff.pc3[app.plot.diff.pc3$trait == "Hole Diameter", 12:38], by = list(app.plot.diff.pc3$trait[app.plot.diff.pc3$trait == "Hole Diameter"], app.plot.diff.pc3$disc.trait[app.plot.diff.pc3$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.app.diff.pc3$Group.1)){
    assign(paste0("app.diff.theta.pc3.", gsub(" ", ".", agg.app.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.app.diff.pc3[i, grep("theta", names(agg.app.diff.pc3))])), ncol = length(states.app), nrow = length(states.app), byrow = TRUE), row.names = states.app))
    assign(paste0("app.diff.alpha.pc3.", gsub(" ", ".", agg.app.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.app.diff.pc3[i, grep("alpha", names(agg.app.diff.pc3))])), ncol = length(states.app), nrow = length(states.app), byrow = TRUE), row.names = states.app))
    assign(paste0("app.diff.sigma.pc3.", gsub(" ", ".", agg.app.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.app.diff.pc3[i, grep("sigma", names(agg.app.diff.pc3))])), ncol = length(states.app), nrow = length(states.app), byrow = TRUE), row.names = states.app))
}







## Architecture
### PC1

theta.arch.plot.pc1 <- subset(theta.plot.pc1, disc.trait == "arch")
theta.arch.plot.pc1$variable <- as.character(theta.arch.plot.pc1$variable)
theta.arch.plot.pc1 <- theta.arch.plot.pc1[theta.arch.plot.pc1$value >= -10 & theta.arch.plot.pc1$value <= 10, ]

alpha.arch.plot.pc1 <- subset(alpha.plot.pc1, disc.trait == "arch")
alpha.arch.plot.pc1$variable <- as.character(alpha.arch.plot.pc1$variable)
alpha.arch.plot.pc1 <- alpha.arch.plot.pc1[alpha.arch.plot.pc1$value >= -10 & alpha.arch.plot.pc1$value <= 10, ]

sigma.arch.plot.pc1 <- subset(sigma.plot.pc1, disc.trait == "arch")
sigma.arch.plot.pc1$variable <- as.character(sigma.arch.plot.pc1$variable)
sigma.arch.plot.pc1 <- sigma.arch.plot.pc1[sigma.arch.plot.pc1$value >= -10 & sigma.arch.plot.pc1$value <= 10, ]


#### Differences
arch.plot.diff.pc1 <- pars.table.pc1[pars.table.pc1$disc.trait == "arch", -c(6, 10, 14, 15)]

for(i in 1:length(states.arch)){
    for(j in 1:length(states.arch)){
        arch.plot.diff.pc1[, paste0("theta", i, j)] <- arch.plot.diff.pc1[, paste0("theta", i)] > arch.plot.diff.pc1[, paste0("theta", j)]
        arch.plot.diff.pc1[, paste0("alpha", i, j)] <- arch.plot.diff.pc1[, paste0("alpha", i)] > arch.plot.diff.pc1[, paste0("alpha", j)]
        arch.plot.diff.pc1[, paste0("sigma", i, j)] <- arch.plot.diff.pc1[, paste0("sigma", i)] > arch.plot.diff.pc1[, paste0("sigma", j)]
    }
}

agg.arch.diff.pc1 <- aggregate(arch.plot.diff.pc1[arch.plot.diff.pc1$trait == "Hole Diameter", 12:ncol(arch.plot.diff.pc1)], by = list(arch.plot.diff.pc1$trait[arch.plot.diff.pc1$trait == "Hole Diameter"], arch.plot.diff.pc1$disc.trait[arch.plot.diff.pc1$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.arch.diff.pc1$Group.1)){
    assign(paste0("arch.diff.theta.pc1.", gsub(" ", ".", agg.arch.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.arch.diff.pc1[i, grep("theta", names(agg.arch.diff.pc1))])), ncol = length(states.arch), nrow = length(states.arch), byrow = TRUE), row.names = states.arch))
    assign(paste0("arch.diff.alpha.pc1.", gsub(" ", ".", agg.arch.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.arch.diff.pc1[i, grep("alpha", names(agg.arch.diff.pc1))])), ncol = length(states.arch), nrow = length(states.arch), byrow = TRUE), row.names = states.arch))
    assign(paste0("arch.diff.sigma.pc1.", gsub(" ", ".", agg.arch.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.arch.diff.pc1[i, grep("sigma", names(agg.arch.diff.pc1))])), ncol = length(states.arch), nrow = length(states.arch), byrow = TRUE), row.names = states.arch))
}

### PC2

theta.arch.plot.pc2 <- subset(theta.plot.pc2, disc.trait == "arch")
theta.arch.plot.pc2$variable <- as.character(theta.arch.plot.pc2$variable)
theta.arch.plot.pc2 <- theta.arch.plot.pc2[theta.arch.plot.pc2$value >= -10 & theta.arch.plot.pc2$value <= 10, ]

alpha.arch.plot.pc2 <- subset(alpha.plot.pc2, disc.trait == "arch")
alpha.arch.plot.pc2$variable <- as.character(alpha.arch.plot.pc2$variable)
alpha.arch.plot.pc2 <- alpha.arch.plot.pc2[alpha.arch.plot.pc2$value >= -10 & alpha.arch.plot.pc2$value <= 10, ]

sigma.arch.plot.pc2 <- subset(sigma.plot.pc2, disc.trait == "arch")
sigma.arch.plot.pc2$variable <- as.character(sigma.arch.plot.pc2$variable)
sigma.arch.plot.pc2 <- sigma.arch.plot.pc2[sigma.arch.plot.pc2$value >= -10 & sigma.arch.plot.pc2$value <= 10, ]

arch.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "arch", -c(6, 10, 14, 15)]
arch.diff.pc2 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(arch.plot.diff.pc2[, grep("theta", names(arch.plot.diff.pc2))] >= arch.plot.diff.pc2$theta1))


#### Differences
arch.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "arch", -c(6, 10, 14, 15)]

for(i in 1:length(states.arch)){
    for(j in 1:length(states.arch)){
        arch.plot.diff.pc2[, paste0("theta", i, j)] <- arch.plot.diff.pc2[, paste0("theta", i)] > arch.plot.diff.pc2[, paste0("theta", j)]
        arch.plot.diff.pc2[, paste0("alpha", i, j)] <- arch.plot.diff.pc2[, paste0("alpha", i)] > arch.plot.diff.pc2[, paste0("alpha", j)]
        arch.plot.diff.pc2[, paste0("sigma", i, j)] <- arch.plot.diff.pc2[, paste0("sigma", i)] > arch.plot.diff.pc2[, paste0("sigma", j)]
    }
}

agg.arch.diff.pc2 <- aggregate(arch.plot.diff.pc2[arch.plot.diff.pc2$trait == "Hole Diameter", 12:ncol(arch.plot.diff.pc2)], by = list(arch.plot.diff.pc2$trait[arch.plot.diff.pc2$trait == "Hole Diameter"], arch.plot.diff.pc2$disc.trait[arch.plot.diff.pc2$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.arch.diff.pc2$Group.1)){
    assign(paste0("arch.diff.theta.pc2.", gsub(" ", ".", agg.arch.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.arch.diff.pc2[i, grep("theta", names(agg.arch.diff.pc2))])), ncol = length(states.arch), nrow = length(states.arch), byrow = TRUE), row.names = states.arch))
    assign(paste0("arch.diff.alpha.pc2.", gsub(" ", ".", agg.arch.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.arch.diff.pc2[i, grep("alpha", names(agg.arch.diff.pc2))])), ncol = length(states.arch), nrow = length(states.arch), byrow = TRUE), row.names = states.arch))
    assign(paste0("arch.diff.sigma.pc2.", gsub(" ", ".", agg.arch.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.arch.diff.pc2[i, grep("sigma", names(agg.arch.diff.pc2))])), ncol = length(states.arch), nrow = length(states.arch), byrow = TRUE), row.names = states.arch))
}


### PC3

theta.arch.plot.pc3 <- subset(theta.plot.pc3, disc.trait == "arch")
theta.arch.plot.pc3$variable <- as.character(theta.arch.plot.pc3$variable)
theta.arch.plot.pc3 <- theta.arch.plot.pc3[theta.arch.plot.pc3$value >= -10 & theta.arch.plot.pc3$value <= 10, ]

alpha.arch.plot.pc3 <- subset(alpha.plot.pc3, disc.trait == "arch")
alpha.arch.plot.pc3$variable <- as.character(alpha.arch.plot.pc3$variable)
alpha.arch.plot.pc3 <- alpha.arch.plot.pc3[alpha.arch.plot.pc3$value >= -10 & alpha.arch.plot.pc3$value <= 10, ]

sigma.arch.plot.pc3 <- subset(sigma.plot.pc3, disc.trait == "arch")
sigma.arch.plot.pc3$variable <- as.character(sigma.arch.plot.pc3$variable)
sigma.arch.plot.pc3 <- sigma.arch.plot.pc3[sigma.arch.plot.pc3$value >= -10 & sigma.arch.plot.pc3$value <= 10, ]

arch.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "arch", -c(6, 10, 14, 15)]
arch.diff.pc3 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(arch.plot.diff.pc3[, grep("theta", names(arch.plot.diff.pc3))] >= arch.plot.diff.pc3$theta1))


#### Differences
arch.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "arch", -c(6, 10, 14, 15)]

for(i in 1:length(states.arch)){
    for(j in 1:length(states.arch)){
        arch.plot.diff.pc3[, paste0("theta", i, j)] <- arch.plot.diff.pc3[, paste0("theta", i)] > arch.plot.diff.pc3[, paste0("theta", j)]
        arch.plot.diff.pc3[, paste0("alpha", i, j)] <- arch.plot.diff.pc3[, paste0("alpha", i)] > arch.plot.diff.pc3[, paste0("alpha", j)]
        arch.plot.diff.pc3[, paste0("sigma", i, j)] <- arch.plot.diff.pc3[, paste0("sigma", i)] > arch.plot.diff.pc3[, paste0("sigma", j)]
    }
}

agg.arch.diff.pc3 <- aggregate(arch.plot.diff.pc3[arch.plot.diff.pc3$trait == "Hole Diameter", 12:ncol(arch.plot.diff.pc3)], by = list(arch.plot.diff.pc3$trait[arch.plot.diff.pc3$trait == "Hole Diameter"], arch.plot.diff.pc3$disc.trait[arch.plot.diff.pc3$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.arch.diff.pc3$Group.1)){
    assign(paste0("arch.diff.theta.pc3.", gsub(" ", ".", agg.arch.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.arch.diff.pc3[i, grep("theta", names(agg.arch.diff.pc3))])), ncol = length(states.arch), nrow = length(states.arch), byrow = TRUE), row.names = states.arch))
    assign(paste0("arch.diff.alpha.pc3.", gsub(" ", ".", agg.arch.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.arch.diff.pc3[i, grep("alpha", names(agg.arch.diff.pc3))])), ncol = length(states.arch), nrow = length(states.arch), byrow = TRUE), row.names = states.arch))
    assign(paste0("arch.diff.sigma.pc3.", gsub(" ", ".", agg.arch.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.arch.diff.pc3[i, grep("sigma", names(agg.arch.diff.pc3))])), ncol = length(states.arch), nrow = length(states.arch), byrow = TRUE), row.names = states.arch))
}






## Domatium Growth
### PC1

theta.domgrow.plot.pc1 <- subset(theta.plot.pc1, disc.trait == "domgrow")
theta.domgrow.plot.pc1$variable <- as.character(theta.domgrow.plot.pc1$variable)
theta.domgrow.plot.pc1 <- theta.domgrow.plot.pc1[theta.domgrow.plot.pc1$value >= -10 & theta.domgrow.plot.pc1$value <= 10, ]

alpha.domgrow.plot.pc1 <- subset(alpha.plot.pc1, disc.trait == "domgrow")
alpha.domgrow.plot.pc1$variable <- as.character(alpha.domgrow.plot.pc1$variable)
alpha.domgrow.plot.pc1 <- alpha.domgrow.plot.pc1[alpha.domgrow.plot.pc1$value >= -10 & alpha.domgrow.plot.pc1$value <= 10, ]

sigma.domgrow.plot.pc1 <- subset(sigma.plot.pc1, disc.trait == "domgrow")
sigma.domgrow.plot.pc1$variable <- as.character(sigma.domgrow.plot.pc1$variable)
sigma.domgrow.plot.pc1 <- sigma.domgrow.plot.pc1[sigma.domgrow.plot.pc1$value >= -10 & sigma.domgrow.plot.pc1$value <= 10, ]


#### Differences
domgrow.plot.diff.pc1 <- pars.table.pc1[pars.table.pc1$disc.trait == "domgrow", -c(6, 10, 14, 15)]

for(i in 1:length(states.domgrow)){
    for(j in 1:length(states.domgrow)){
        domgrow.plot.diff.pc1[, paste0("theta", i, j)] <- domgrow.plot.diff.pc1[, paste0("theta", i)] > domgrow.plot.diff.pc1[, paste0("theta", j)]
        domgrow.plot.diff.pc1[, paste0("alpha", i, j)] <- domgrow.plot.diff.pc1[, paste0("alpha", i)] > domgrow.plot.diff.pc1[, paste0("alpha", j)]
        domgrow.plot.diff.pc1[, paste0("sigma", i, j)] <- domgrow.plot.diff.pc1[, paste0("sigma", i)] > domgrow.plot.diff.pc1[, paste0("sigma", j)]
    }
}

agg.domgrow.diff.pc1 <- aggregate(domgrow.plot.diff.pc1[domgrow.plot.diff.pc1$trait == "Hole Diameter", 12:ncol(domgrow.plot.diff.pc1)], by = list(domgrow.plot.diff.pc1$trait[domgrow.plot.diff.pc1$trait == "Hole Diameter"], domgrow.plot.diff.pc1$disc.trait[domgrow.plot.diff.pc1$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.domgrow.diff.pc1$Group.1)){
    assign(paste0("domgrow.diff.theta.pc1.", gsub(" ", ".", agg.domgrow.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.domgrow.diff.pc1[i, grep("theta", names(agg.domgrow.diff.pc1))])), ncol = length(states.domgrow), nrow = length(states.domgrow), byrow = TRUE), row.names = states.domgrow))
    assign(paste0("domgrow.diff.alpha.pc1.", gsub(" ", ".", agg.domgrow.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.domgrow.diff.pc1[i, grep("alpha", names(agg.domgrow.diff.pc1))])), ncol = length(states.domgrow), nrow = length(states.domgrow), byrow = TRUE), row.names = states.domgrow))
    assign(paste0("domgrow.diff.sigma.pc1.", gsub(" ", ".", agg.domgrow.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.domgrow.diff.pc1[i, grep("sigma", names(agg.domgrow.diff.pc1))])), ncol = length(states.domgrow), nrow = length(states.domgrow), byrow = TRUE), row.names = states.domgrow))
}

### PC2

theta.domgrow.plot.pc2 <- subset(theta.plot.pc2, disc.trait == "domgrow")
theta.domgrow.plot.pc2$variable <- as.character(theta.domgrow.plot.pc2$variable)
theta.domgrow.plot.pc2 <- theta.domgrow.plot.pc2[theta.domgrow.plot.pc2$value >= -10 & theta.domgrow.plot.pc2$value <= 10, ]

alpha.domgrow.plot.pc2 <- subset(alpha.plot.pc2, disc.trait == "domgrow")
alpha.domgrow.plot.pc2$variable <- as.character(alpha.domgrow.plot.pc2$variable)
alpha.domgrow.plot.pc2 <- alpha.domgrow.plot.pc2[alpha.domgrow.plot.pc2$value >= -10 & alpha.domgrow.plot.pc2$value <= 10, ]

sigma.domgrow.plot.pc2 <- subset(sigma.plot.pc2, disc.trait == "domgrow")
sigma.domgrow.plot.pc2$variable <- as.character(sigma.domgrow.plot.pc2$variable)
sigma.domgrow.plot.pc2 <- sigma.domgrow.plot.pc2[sigma.domgrow.plot.pc2$value >= -10 & sigma.domgrow.plot.pc2$value <= 10, ]

domgrow.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "domgrow", -c(6, 10, 14, 15)]
domgrow.diff.pc2 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(domgrow.plot.diff.pc2[, grep("theta", names(domgrow.plot.diff.pc2))] >= domgrow.plot.diff.pc2$theta1))


#### Differences
domgrow.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "domgrow", -c(6, 10, 14, 15)]

for(i in 1:length(states.domgrow)){
    for(j in 1:length(states.domgrow)){
        domgrow.plot.diff.pc2[, paste0("theta", i, j)] <- domgrow.plot.diff.pc2[, paste0("theta", i)] > domgrow.plot.diff.pc2[, paste0("theta", j)]
        domgrow.plot.diff.pc2[, paste0("alpha", i, j)] <- domgrow.plot.diff.pc2[, paste0("alpha", i)] > domgrow.plot.diff.pc2[, paste0("alpha", j)]
        domgrow.plot.diff.pc2[, paste0("sigma", i, j)] <- domgrow.plot.diff.pc2[, paste0("sigma", i)] > domgrow.plot.diff.pc2[, paste0("sigma", j)]
    }
}

agg.domgrow.diff.pc2 <- aggregate(domgrow.plot.diff.pc2[domgrow.plot.diff.pc2$trait == "Hole Diameter", 12:ncol(domgrow.plot.diff.pc2)], by = list(domgrow.plot.diff.pc2$trait[domgrow.plot.diff.pc2$trait == "Hole Diameter"], domgrow.plot.diff.pc2$disc.trait[domgrow.plot.diff.pc2$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.domgrow.diff.pc2$Group.1)){
    assign(paste0("domgrow.diff.theta.pc2.", gsub(" ", ".", agg.domgrow.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.domgrow.diff.pc2[i, grep("theta", names(agg.domgrow.diff.pc2))])), ncol = length(states.domgrow), nrow = length(states.domgrow), byrow = TRUE), row.names = states.domgrow))
    assign(paste0("domgrow.diff.alpha.pc2.", gsub(" ", ".", agg.domgrow.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.domgrow.diff.pc2[i, grep("alpha", names(agg.domgrow.diff.pc2))])), ncol = length(states.domgrow), nrow = length(states.domgrow), byrow = TRUE), row.names = states.domgrow))
    assign(paste0("domgrow.diff.sigma.pc2.", gsub(" ", ".", agg.domgrow.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.domgrow.diff.pc2[i, grep("sigma", names(agg.domgrow.diff.pc2))])), ncol = length(states.domgrow), nrow = length(states.domgrow), byrow = TRUE), row.names = states.domgrow))
}


### PC3

theta.domgrow.plot.pc3 <- subset(theta.plot.pc3, disc.trait == "domgrow")
theta.domgrow.plot.pc3$variable <- as.character(theta.domgrow.plot.pc3$variable)
theta.domgrow.plot.pc3 <- theta.domgrow.plot.pc3[theta.domgrow.plot.pc3$value >= -10 & theta.domgrow.plot.pc3$value <= 10, ]

alpha.domgrow.plot.pc3 <- subset(alpha.plot.pc3, disc.trait == "domgrow")
alpha.domgrow.plot.pc3$variable <- as.character(alpha.domgrow.plot.pc3$variable)
alpha.domgrow.plot.pc3 <- alpha.domgrow.plot.pc3[alpha.domgrow.plot.pc3$value >= -10 & alpha.domgrow.plot.pc3$value <= 10, ]

sigma.domgrow.plot.pc3 <- subset(sigma.plot.pc3, disc.trait == "domgrow")
sigma.domgrow.plot.pc3$variable <- as.character(sigma.domgrow.plot.pc3$variable)
sigma.domgrow.plot.pc3 <- sigma.domgrow.plot.pc3[sigma.domgrow.plot.pc3$value >= -10 & sigma.domgrow.plot.pc3$value <= 10, ]

domgrow.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "domgrow", -c(6, 10, 14, 15)]
domgrow.diff.pc3 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(domgrow.plot.diff.pc3[, grep("theta", names(domgrow.plot.diff.pc3))] >= domgrow.plot.diff.pc3$theta1))


#### Differences
domgrow.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "domgrow", -c(6, 10, 14, 15)]

for(i in 1:length(states.domgrow)){
    for(j in 1:length(states.domgrow)){
        domgrow.plot.diff.pc3[, paste0("theta", i, j)] <- domgrow.plot.diff.pc3[, paste0("theta", i)] > domgrow.plot.diff.pc3[, paste0("theta", j)]
        domgrow.plot.diff.pc3[, paste0("alpha", i, j)] <- domgrow.plot.diff.pc3[, paste0("alpha", i)] > domgrow.plot.diff.pc3[, paste0("alpha", j)]
        domgrow.plot.diff.pc3[, paste0("sigma", i, j)] <- domgrow.plot.diff.pc3[, paste0("sigma", i)] > domgrow.plot.diff.pc3[, paste0("sigma", j)]
    }
}

agg.domgrow.diff.pc3 <- aggregate(domgrow.plot.diff.pc3[domgrow.plot.diff.pc3$trait == "Hole Diameter", 12:ncol(domgrow.plot.diff.pc3)], by = list(domgrow.plot.diff.pc3$trait[domgrow.plot.diff.pc3$trait == "Hole Diameter"], domgrow.plot.diff.pc3$disc.trait[domgrow.plot.diff.pc3$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.domgrow.diff.pc3$Group.1)){
    assign(paste0("domgrow.diff.theta.pc3.", gsub(" ", ".", agg.domgrow.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.domgrow.diff.pc3[i, grep("theta", names(agg.domgrow.diff.pc3))])), ncol = length(states.domgrow), nrow = length(states.domgrow), byrow = TRUE), row.names = states.domgrow))
    assign(paste0("domgrow.diff.alpha.pc3.", gsub(" ", ".", agg.domgrow.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.domgrow.diff.pc3[i, grep("alpha", names(agg.domgrow.diff.pc3))])), ncol = length(states.domgrow), nrow = length(states.domgrow), byrow = TRUE), row.names = states.domgrow))
    assign(paste0("domgrow.diff.sigma.pc3.", gsub(" ", ".", agg.domgrow.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.domgrow.diff.pc3[i, grep("sigma", names(agg.domgrow.diff.pc3))])), ncol = length(states.domgrow), nrow = length(states.domgrow), byrow = TRUE), row.names = states.domgrow))
}






## Leaf Structure
### PC1

theta.leafstruc.plot.pc1 <- subset(theta.plot.pc1, disc.trait == "leafstruc")
theta.leafstruc.plot.pc1$variable <- as.character(theta.leafstruc.plot.pc1$variable)
theta.leafstruc.plot.pc1 <- theta.leafstruc.plot.pc1[theta.leafstruc.plot.pc1$value >= -10 & theta.leafstruc.plot.pc1$value <= 10, ]

alpha.leafstruc.plot.pc1 <- subset(alpha.plot.pc1, disc.trait == "leafstruc")
alpha.leafstruc.plot.pc1$variable <- as.character(alpha.leafstruc.plot.pc1$variable)
alpha.leafstruc.plot.pc1 <- alpha.leafstruc.plot.pc1[alpha.leafstruc.plot.pc1$value >= -10 & alpha.leafstruc.plot.pc1$value <= 10, ]

sigma.leafstruc.plot.pc1 <- subset(sigma.plot.pc1, disc.trait == "leafstruc")
sigma.leafstruc.plot.pc1$variable <- as.character(sigma.leafstruc.plot.pc1$variable)
sigma.leafstruc.plot.pc1 <- sigma.leafstruc.plot.pc1[sigma.leafstruc.plot.pc1$value >= -10 & sigma.leafstruc.plot.pc1$value <= 10, ]


#### Differences
leafstruc.plot.diff.pc1 <- pars.table.pc1[pars.table.pc1$disc.trait == "leafstruc",]

for(i in 1:length(states.leafstruc)){
    for(j in 1:length(states.leafstruc)){
        leafstruc.plot.diff.pc1[, paste0("theta", i, j)] <- leafstruc.plot.diff.pc1[, paste0("theta", i)] > leafstruc.plot.diff.pc1[, paste0("theta", j)]
        leafstruc.plot.diff.pc1[, paste0("alpha", i, j)] <- leafstruc.plot.diff.pc1[, paste0("alpha", i)] > leafstruc.plot.diff.pc1[, paste0("alpha", j)]
        leafstruc.plot.diff.pc1[, paste0("sigma", i, j)] <- leafstruc.plot.diff.pc1[, paste0("sigma", i)] > leafstruc.plot.diff.pc1[, paste0("sigma", j)]
    }
}

agg.leafstruc.diff.pc1 <- aggregate(leafstruc.plot.diff.pc1[leafstruc.plot.diff.pc1$trait == "Hole Diameter", 16:ncol(leafstruc.plot.diff.pc1)], by = list(leafstruc.plot.diff.pc1$trait[leafstruc.plot.diff.pc1$trait == "Hole Diameter"], leafstruc.plot.diff.pc1$disc.trait[leafstruc.plot.diff.pc1$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.leafstruc.diff.pc1$Group.1)){
    assign(paste0("leafstruc.diff.theta.pc1.", gsub(" ", ".", agg.leafstruc.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.leafstruc.diff.pc1[i, grep("theta", names(agg.leafstruc.diff.pc1))])), ncol = length(states.leafstruc), nrow = length(states.leafstruc), byrow = TRUE), row.names = states.leafstruc))
    assign(paste0("leafstruc.diff.alpha.pc1.", gsub(" ", ".", agg.leafstruc.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.leafstruc.diff.pc1[i, grep("alpha", names(agg.leafstruc.diff.pc1))])), ncol = length(states.leafstruc), nrow = length(states.leafstruc), byrow = TRUE), row.names = states.leafstruc))
    assign(paste0("leafstruc.diff.sigma.pc1.", gsub(" ", ".", agg.leafstruc.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.leafstruc.diff.pc1[i, grep("sigma", names(agg.leafstruc.diff.pc1))])), ncol = length(states.leafstruc), nrow = length(states.leafstruc), byrow = TRUE), row.names = states.leafstruc))
}

### PC2

theta.leafstruc.plot.pc2 <- subset(theta.plot.pc2, disc.trait == "leafstruc")
theta.leafstruc.plot.pc2$variable <- as.character(theta.leafstruc.plot.pc2$variable)
theta.leafstruc.plot.pc2 <- theta.leafstruc.plot.pc2[theta.leafstruc.plot.pc2$value >= -10 & theta.leafstruc.plot.pc2$value <= 10, ]

alpha.leafstruc.plot.pc2 <- subset(alpha.plot.pc2, disc.trait == "leafstruc")
alpha.leafstruc.plot.pc2$variable <- as.character(alpha.leafstruc.plot.pc2$variable)
alpha.leafstruc.plot.pc2 <- alpha.leafstruc.plot.pc2[alpha.leafstruc.plot.pc2$value >= -10 & alpha.leafstruc.plot.pc2$value <= 10, ]

sigma.leafstruc.plot.pc2 <- subset(sigma.plot.pc2, disc.trait == "leafstruc")
sigma.leafstruc.plot.pc2$variable <- as.character(sigma.leafstruc.plot.pc2$variable)
sigma.leafstruc.plot.pc2 <- sigma.leafstruc.plot.pc2[sigma.leafstruc.plot.pc2$value >= -10 & sigma.leafstruc.plot.pc2$value <= 10, ]

leafstruc.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "leafstruc", ]
leafstruc.diff.pc2 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(leafstruc.plot.diff.pc2[, grep("theta", names(leafstruc.plot.diff.pc2))] >= leafstruc.plot.diff.pc2$theta1))


#### Differences
leafstruc.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "leafstruc", ]

for(i in 1:length(states.leafstruc)){
    for(j in 1:length(states.leafstruc)){
        leafstruc.plot.diff.pc2[, paste0("theta", i, j)] <- leafstruc.plot.diff.pc2[, paste0("theta", i)] > leafstruc.plot.diff.pc2[, paste0("theta", j)]
        leafstruc.plot.diff.pc2[, paste0("alpha", i, j)] <- leafstruc.plot.diff.pc2[, paste0("alpha", i)] > leafstruc.plot.diff.pc2[, paste0("alpha", j)]
        leafstruc.plot.diff.pc2[, paste0("sigma", i, j)] <- leafstruc.plot.diff.pc2[, paste0("sigma", i)] > leafstruc.plot.diff.pc2[, paste0("sigma", j)]
    }
}

agg.leafstruc.diff.pc2 <- aggregate(leafstruc.plot.diff.pc2[leafstruc.plot.diff.pc2$trait == "Hole Diameter", 16:ncol(leafstruc.plot.diff.pc2)], by = list(leafstruc.plot.diff.pc2$trait[leafstruc.plot.diff.pc2$trait == "Hole Diameter"], leafstruc.plot.diff.pc2$disc.trait[leafstruc.plot.diff.pc2$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.leafstruc.diff.pc2$Group.1)){
    assign(paste0("leafstruc.diff.theta.pc2.", gsub(" ", ".", agg.leafstruc.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.leafstruc.diff.pc2[i, grep("theta", names(agg.leafstruc.diff.pc2))])), ncol = length(states.leafstruc), nrow = length(states.leafstruc), byrow = TRUE), row.names = states.leafstruc))
    assign(paste0("leafstruc.diff.alpha.pc2.", gsub(" ", ".", agg.leafstruc.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.leafstruc.diff.pc2[i, grep("alpha", names(agg.leafstruc.diff.pc2))])), ncol = length(states.leafstruc), nrow = length(states.leafstruc), byrow = TRUE), row.names = states.leafstruc))
    assign(paste0("leafstruc.diff.sigma.pc2.", gsub(" ", ".", agg.leafstruc.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.leafstruc.diff.pc2[i, grep("sigma", names(agg.leafstruc.diff.pc2))])), ncol = length(states.leafstruc), nrow = length(states.leafstruc), byrow = TRUE), row.names = states.leafstruc))
}


### PC3

theta.leafstruc.plot.pc3 <- subset(theta.plot.pc3, disc.trait == "leafstruc")
theta.leafstruc.plot.pc3$variable <- as.character(theta.leafstruc.plot.pc3$variable)
theta.leafstruc.plot.pc3 <- theta.leafstruc.plot.pc3[theta.leafstruc.plot.pc3$value >= -10 & theta.leafstruc.plot.pc3$value <= 10, ]

alpha.leafstruc.plot.pc3 <- subset(alpha.plot.pc3, disc.trait == "leafstruc")
alpha.leafstruc.plot.pc3$variable <- as.character(alpha.leafstruc.plot.pc3$variable)
alpha.leafstruc.plot.pc3 <- alpha.leafstruc.plot.pc3[alpha.leafstruc.plot.pc3$value >= -10 & alpha.leafstruc.plot.pc3$value <= 10, ]

sigma.leafstruc.plot.pc3 <- subset(sigma.plot.pc3, disc.trait == "leafstruc")
sigma.leafstruc.plot.pc3$variable <- as.character(sigma.leafstruc.plot.pc3$variable)
sigma.leafstruc.plot.pc3 <- sigma.leafstruc.plot.pc3[sigma.leafstruc.plot.pc3$value >= -10 & sigma.leafstruc.plot.pc3$value <= 10, ]

leafstruc.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "leafstruc", ]
leafstruc.diff.pc3 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(leafstruc.plot.diff.pc3[, grep("theta", names(leafstruc.plot.diff.pc3))] >= leafstruc.plot.diff.pc3$theta1))


#### Differences
leafstruc.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "leafstruc", ]

for(i in 1:length(states.leafstruc)){
    for(j in 1:length(states.leafstruc)){
        leafstruc.plot.diff.pc3[, paste0("theta", i, j)] <- leafstruc.plot.diff.pc3[, paste0("theta", i)] > leafstruc.plot.diff.pc3[, paste0("theta", j)]
        leafstruc.plot.diff.pc3[, paste0("alpha", i, j)] <- leafstruc.plot.diff.pc3[, paste0("alpha", i)] > leafstruc.plot.diff.pc3[, paste0("alpha", j)]
        leafstruc.plot.diff.pc3[, paste0("sigma", i, j)] <- leafstruc.plot.diff.pc3[, paste0("sigma", i)] > leafstruc.plot.diff.pc3[, paste0("sigma", j)]
    }
}

agg.leafstruc.diff.pc3 <- aggregate(leafstruc.plot.diff.pc3[leafstruc.plot.diff.pc3$trait == "Hole Diameter", 16:ncol(leafstruc.plot.diff.pc3)], by = list(leafstruc.plot.diff.pc3$trait[leafstruc.plot.diff.pc3$trait == "Hole Diameter"], leafstruc.plot.diff.pc3$disc.trait[leafstruc.plot.diff.pc3$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.leafstruc.diff.pc3$Group.1)){
    assign(paste0("leafstruc.diff.theta.pc3.", gsub(" ", ".", agg.leafstruc.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.leafstruc.diff.pc3[i, grep("theta", names(agg.leafstruc.diff.pc3))])), ncol = length(states.leafstruc), nrow = length(states.leafstruc), byrow = TRUE), row.names = states.leafstruc))
    assign(paste0("leafstruc.diff.alpha.pc3.", gsub(" ", ".", agg.leafstruc.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.leafstruc.diff.pc3[i, grep("alpha", names(agg.leafstruc.diff.pc3))])), ncol = length(states.leafstruc), nrow = length(states.leafstruc), byrow = TRUE), row.names = states.leafstruc))
    assign(paste0("leafstruc.diff.sigma.pc3.", gsub(" ", ".", agg.leafstruc.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.leafstruc.diff.pc3[i, grep("sigma", names(agg.leafstruc.diff.pc3))])), ncol = length(states.leafstruc), nrow = length(states.leafstruc), byrow = TRUE), row.names = states.leafstruc))
}






## Mating System
### PC1

theta.matsys.plot.pc1 <- subset(theta.plot.pc1, disc.trait == "matsys")
theta.matsys.plot.pc1$variable <- as.character(theta.matsys.plot.pc1$variable)
theta.matsys.plot.pc1 <- theta.matsys.plot.pc1[theta.matsys.plot.pc1$value >= -10 & theta.matsys.plot.pc1$value <= 10, ]

alpha.matsys.plot.pc1 <- subset(alpha.plot.pc1, disc.trait == "matsys")
alpha.matsys.plot.pc1$variable <- as.character(alpha.matsys.plot.pc1$variable)
alpha.matsys.plot.pc1 <- alpha.matsys.plot.pc1[alpha.matsys.plot.pc1$value >= -10 & alpha.matsys.plot.pc1$value <= 10, ]

sigma.matsys.plot.pc1 <- subset(sigma.plot.pc1, disc.trait == "matsys")
sigma.matsys.plot.pc1$variable <- as.character(sigma.matsys.plot.pc1$variable)
sigma.matsys.plot.pc1 <- sigma.matsys.plot.pc1[sigma.matsys.plot.pc1$value >= -10 & sigma.matsys.plot.pc1$value <= 10, ]


#### Differences
matsys.plot.diff.pc1 <- pars.table.pc1[pars.table.pc1$disc.trait == "matsys", -c(6, 10, 14, 15)]

for(i in 1:length(states.matsys)){
    for(j in 1:length(states.matsys)){
        matsys.plot.diff.pc1[, paste0("theta", i, j)] <- matsys.plot.diff.pc1[, paste0("theta", i)] > matsys.plot.diff.pc1[, paste0("theta", j)]
        matsys.plot.diff.pc1[, paste0("alpha", i, j)] <- matsys.plot.diff.pc1[, paste0("alpha", i)] > matsys.plot.diff.pc1[, paste0("alpha", j)]
        matsys.plot.diff.pc1[, paste0("sigma", i, j)] <- matsys.plot.diff.pc1[, paste0("sigma", i)] > matsys.plot.diff.pc1[, paste0("sigma", j)]
    }
}

agg.matsys.diff.pc1 <- aggregate(matsys.plot.diff.pc1[matsys.plot.diff.pc1$trait == "Hole Diameter", 12:ncol(matsys.plot.diff.pc1)], by = list(matsys.plot.diff.pc1$trait[matsys.plot.diff.pc1$trait == "Hole Diameter"], matsys.plot.diff.pc1$disc.trait[matsys.plot.diff.pc1$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.matsys.diff.pc1$Group.1)){
    assign(paste0("matsys.diff.theta.pc1.", gsub(" ", ".", agg.matsys.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.matsys.diff.pc1[i, grep("theta", names(agg.matsys.diff.pc1))])), ncol = length(states.matsys), nrow = length(states.matsys), byrow = TRUE), row.names = states.matsys))
    assign(paste0("matsys.diff.alpha.pc1.", gsub(" ", ".", agg.matsys.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.matsys.diff.pc1[i, grep("alpha", names(agg.matsys.diff.pc1))])), ncol = length(states.matsys), nrow = length(states.matsys), byrow = TRUE), row.names = states.matsys))
    assign(paste0("matsys.diff.sigma.pc1.", gsub(" ", ".", agg.matsys.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.matsys.diff.pc1[i, grep("sigma", names(agg.matsys.diff.pc1))])), ncol = length(states.matsys), nrow = length(states.matsys), byrow = TRUE), row.names = states.matsys))
}

### PC2

theta.matsys.plot.pc2 <- subset(theta.plot.pc2, disc.trait == "matsys")
theta.matsys.plot.pc2$variable <- as.character(theta.matsys.plot.pc2$variable)
theta.matsys.plot.pc2 <- theta.matsys.plot.pc2[theta.matsys.plot.pc2$value >= -10 & theta.matsys.plot.pc2$value <= 10, ]

alpha.matsys.plot.pc2 <- subset(alpha.plot.pc2, disc.trait == "matsys")
alpha.matsys.plot.pc2$variable <- as.character(alpha.matsys.plot.pc2$variable)
alpha.matsys.plot.pc2 <- alpha.matsys.plot.pc2[alpha.matsys.plot.pc2$value >= -10 & alpha.matsys.plot.pc2$value <= 10, ]

sigma.matsys.plot.pc2 <- subset(sigma.plot.pc2, disc.trait == "matsys")
sigma.matsys.plot.pc2$variable <- as.character(sigma.matsys.plot.pc2$variable)
sigma.matsys.plot.pc2 <- sigma.matsys.plot.pc2[sigma.matsys.plot.pc2$value >= -10 & sigma.matsys.plot.pc2$value <= 10, ]

matsys.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "matsys", -c(6, 10, 14, 15)]
matsys.diff.pc2 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(matsys.plot.diff.pc2[, grep("theta", names(matsys.plot.diff.pc2))] >= matsys.plot.diff.pc2$theta1))


#### Differences
matsys.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "matsys", -c(6, 10, 14, 15)]

for(i in 1:length(states.matsys)){
    for(j in 1:length(states.matsys)){
        matsys.plot.diff.pc2[, paste0("theta", i, j)] <- matsys.plot.diff.pc2[, paste0("theta", i)] > matsys.plot.diff.pc2[, paste0("theta", j)]
        matsys.plot.diff.pc2[, paste0("alpha", i, j)] <- matsys.plot.diff.pc2[, paste0("alpha", i)] > matsys.plot.diff.pc2[, paste0("alpha", j)]
        matsys.plot.diff.pc2[, paste0("sigma", i, j)] <- matsys.plot.diff.pc2[, paste0("sigma", i)] > matsys.plot.diff.pc2[, paste0("sigma", j)]
    }
}

agg.matsys.diff.pc2 <- aggregate(matsys.plot.diff.pc2[matsys.plot.diff.pc2$trait == "Hole Diameter", 12:ncol(matsys.plot.diff.pc2)], by = list(matsys.plot.diff.pc2$trait[matsys.plot.diff.pc2$trait == "Hole Diameter"], matsys.plot.diff.pc2$disc.trait[matsys.plot.diff.pc2$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.matsys.diff.pc2$Group.1)){
    assign(paste0("matsys.diff.theta.pc2.", gsub(" ", ".", agg.matsys.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.matsys.diff.pc2[i, grep("theta", names(agg.matsys.diff.pc2))])), ncol = length(states.matsys), nrow = length(states.matsys), byrow = TRUE), row.names = states.matsys))
    assign(paste0("matsys.diff.alpha.pc2.", gsub(" ", ".", agg.matsys.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.matsys.diff.pc2[i, grep("alpha", names(agg.matsys.diff.pc2))])), ncol = length(states.matsys), nrow = length(states.matsys), byrow = TRUE), row.names = states.matsys))
    assign(paste0("matsys.diff.sigma.pc2.", gsub(" ", ".", agg.matsys.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.matsys.diff.pc2[i, grep("sigma", names(agg.matsys.diff.pc2))])), ncol = length(states.matsys), nrow = length(states.matsys), byrow = TRUE), row.names = states.matsys))
}


### PC3

theta.matsys.plot.pc3 <- subset(theta.plot.pc3, disc.trait == "matsys")
theta.matsys.plot.pc3$variable <- as.character(theta.matsys.plot.pc3$variable)
theta.matsys.plot.pc3 <- theta.matsys.plot.pc3[theta.matsys.plot.pc3$value >= -10 & theta.matsys.plot.pc3$value <= 10, ]

alpha.matsys.plot.pc3 <- subset(alpha.plot.pc3, disc.trait == "matsys")
alpha.matsys.plot.pc3$variable <- as.character(alpha.matsys.plot.pc3$variable)
alpha.matsys.plot.pc3 <- alpha.matsys.plot.pc3[alpha.matsys.plot.pc3$value >= -10 & alpha.matsys.plot.pc3$value <= 10, ]

sigma.matsys.plot.pc3 <- subset(sigma.plot.pc3, disc.trait == "matsys")
sigma.matsys.plot.pc3$variable <- as.character(sigma.matsys.plot.pc3$variable)
sigma.matsys.plot.pc3 <- sigma.matsys.plot.pc3[sigma.matsys.plot.pc3$value >= -10 & sigma.matsys.plot.pc3$value <= 10, ]

matsys.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "matsys", -c(6, 10, 14, 15)]
matsys.diff.pc3 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(matsys.plot.diff.pc3[, grep("theta", names(matsys.plot.diff.pc3))] >= matsys.plot.diff.pc3$theta1))


#### Differences
matsys.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "matsys", -c(6, 10, 14, 15)]

for(i in 1:length(states.matsys)){
    for(j in 1:length(states.matsys)){
        matsys.plot.diff.pc3[, paste0("theta", i, j)] <- matsys.plot.diff.pc3[, paste0("theta", i)] > matsys.plot.diff.pc3[, paste0("theta", j)]
        matsys.plot.diff.pc3[, paste0("alpha", i, j)] <- matsys.plot.diff.pc3[, paste0("alpha", i)] > matsys.plot.diff.pc3[, paste0("alpha", j)]
        matsys.plot.diff.pc3[, paste0("sigma", i, j)] <- matsys.plot.diff.pc3[, paste0("sigma", i)] > matsys.plot.diff.pc3[, paste0("sigma", j)]
    }
}

agg.matsys.diff.pc3 <- aggregate(matsys.plot.diff.pc3[matsys.plot.diff.pc3$trait == "Hole Diameter", 12:ncol(matsys.plot.diff.pc3)], by = list(matsys.plot.diff.pc3$trait[matsys.plot.diff.pc3$trait == "Hole Diameter"], matsys.plot.diff.pc3$disc.trait[matsys.plot.diff.pc3$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.matsys.diff.pc3$Group.1)){
    assign(paste0("matsys.diff.theta.pc3.", gsub(" ", ".", agg.matsys.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.matsys.diff.pc3[i, grep("theta", names(agg.matsys.diff.pc3))])), ncol = length(states.matsys), nrow = length(states.matsys), byrow = TRUE), row.names = states.matsys))
    assign(paste0("matsys.diff.alpha.pc3.", gsub(" ", ".", agg.matsys.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.matsys.diff.pc3[i, grep("alpha", names(agg.matsys.diff.pc3))])), ncol = length(states.matsys), nrow = length(states.matsys), byrow = TRUE), row.names = states.matsys))
    assign(paste0("matsys.diff.sigma.pc3.", gsub(" ", ".", agg.matsys.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.matsys.diff.pc3[i, grep("sigma", names(agg.matsys.diff.pc3))])), ncol = length(states.matsys), nrow = length(states.matsys), byrow = TRUE), row.names = states.matsys))
}








## Reward
### PC1

theta.reward.plot.pc1 <- subset(theta.plot.pc1, disc.trait == "reward")
theta.reward.plot.pc1$variable <- as.character(theta.reward.plot.pc1$variable)
theta.reward.plot.pc1 <- theta.reward.plot.pc1[theta.reward.plot.pc1$value >= -10 & theta.reward.plot.pc1$value <= 10, ]

alpha.reward.plot.pc1 <- subset(alpha.plot.pc1, disc.trait == "reward")
alpha.reward.plot.pc1$variable <- as.character(alpha.reward.plot.pc1$variable)
alpha.reward.plot.pc1 <- alpha.reward.plot.pc1[alpha.reward.plot.pc1$value >= -10 & alpha.reward.plot.pc1$value <= 10, ]

sigma.reward.plot.pc1 <- subset(sigma.plot.pc1, disc.trait == "reward")
sigma.reward.plot.pc1$variable <- as.character(sigma.reward.plot.pc1$variable)
sigma.reward.plot.pc1 <- sigma.reward.plot.pc1[sigma.reward.plot.pc1$value >= -10 & sigma.reward.plot.pc1$value <= 10, ]


#### Differences
reward.plot.diff.pc1 <- pars.table.pc1[pars.table.pc1$disc.trait == "reward", -c(6, 10, 14, 15)]

for(i in 1:length(states.reward)){
    for(j in 1:length(states.reward)){
        reward.plot.diff.pc1[, paste0("theta", i, j)] <- reward.plot.diff.pc1[, paste0("theta", i)] > reward.plot.diff.pc1[, paste0("theta", j)]
        reward.plot.diff.pc1[, paste0("alpha", i, j)] <- reward.plot.diff.pc1[, paste0("alpha", i)] > reward.plot.diff.pc1[, paste0("alpha", j)]
        reward.plot.diff.pc1[, paste0("sigma", i, j)] <- reward.plot.diff.pc1[, paste0("sigma", i)] > reward.plot.diff.pc1[, paste0("sigma", j)]
    }
}

agg.reward.diff.pc1 <- aggregate(reward.plot.diff.pc1[reward.plot.diff.pc1$trait == "Hole Diameter", 12:ncol(reward.plot.diff.pc1)], by = list(reward.plot.diff.pc1$trait[reward.plot.diff.pc1$trait == "Hole Diameter"], reward.plot.diff.pc1$disc.trait[reward.plot.diff.pc1$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.reward.diff.pc1$Group.1)){
    assign(paste0("reward.diff.theta.pc1.", gsub(" ", ".", agg.reward.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.reward.diff.pc1[i, grep("theta", names(agg.reward.diff.pc1))])), ncol = length(states.reward), nrow = length(states.reward), byrow = TRUE), row.names = states.reward))
    assign(paste0("reward.diff.alpha.pc1.", gsub(" ", ".", agg.reward.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.reward.diff.pc1[i, grep("alpha", names(agg.reward.diff.pc1))])), ncol = length(states.reward), nrow = length(states.reward), byrow = TRUE), row.names = states.reward))
    assign(paste0("reward.diff.sigma.pc1.", gsub(" ", ".", agg.reward.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.reward.diff.pc1[i, grep("sigma", names(agg.reward.diff.pc1))])), ncol = length(states.reward), nrow = length(states.reward), byrow = TRUE), row.names = states.reward))
}

### PC2

theta.reward.plot.pc2 <- subset(theta.plot.pc2, disc.trait == "reward")
theta.reward.plot.pc2$variable <- as.character(theta.reward.plot.pc2$variable)
theta.reward.plot.pc2 <- theta.reward.plot.pc2[theta.reward.plot.pc2$value >= -10 & theta.reward.plot.pc2$value <= 10, ]

alpha.reward.plot.pc2 <- subset(alpha.plot.pc2, disc.trait == "reward")
alpha.reward.plot.pc2$variable <- as.character(alpha.reward.plot.pc2$variable)
alpha.reward.plot.pc2 <- alpha.reward.plot.pc2[alpha.reward.plot.pc2$value >= -10 & alpha.reward.plot.pc2$value <= 10, ]

sigma.reward.plot.pc2 <- subset(sigma.plot.pc2, disc.trait == "reward")
sigma.reward.plot.pc2$variable <- as.character(sigma.reward.plot.pc2$variable)
sigma.reward.plot.pc2 <- sigma.reward.plot.pc2[sigma.reward.plot.pc2$value >= -10 & sigma.reward.plot.pc2$value <= 10, ]

reward.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "reward", -c(6, 10, 14, 15)]
reward.diff.pc2 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(reward.plot.diff.pc2[, grep("theta", names(reward.plot.diff.pc2))] >= reward.plot.diff.pc2$theta1))


#### Differences
reward.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "reward", -c(6, 10, 14, 15)]

for(i in 1:length(states.reward)){
    for(j in 1:length(states.reward)){
        reward.plot.diff.pc2[, paste0("theta", i, j)] <- reward.plot.diff.pc2[, paste0("theta", i)] > reward.plot.diff.pc2[, paste0("theta", j)]
        reward.plot.diff.pc2[, paste0("alpha", i, j)] <- reward.plot.diff.pc2[, paste0("alpha", i)] > reward.plot.diff.pc2[, paste0("alpha", j)]
        reward.plot.diff.pc2[, paste0("sigma", i, j)] <- reward.plot.diff.pc2[, paste0("sigma", i)] > reward.plot.diff.pc2[, paste0("sigma", j)]
    }
}

agg.reward.diff.pc2 <- aggregate(reward.plot.diff.pc2[reward.plot.diff.pc3$trait == "Hole Diameter", 12:ncol(reward.plot.diff.pc2)], by = list(reward.plot.diff.pc2$trait[reward.plot.diff.pc1$trait == "Hole Diameter"], reward.plot.diff.pc2$disc.trait[reward.plot.diff.pc1$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.reward.diff.pc2$Group.1)){
    assign(paste0("reward.diff.theta.pc2.", gsub(" ", ".", agg.reward.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.reward.diff.pc2[i, grep("theta", names(agg.reward.diff.pc2))])), ncol = length(states.reward), nrow = length(states.reward), byrow = TRUE), row.names = states.reward))
    assign(paste0("reward.diff.alpha.pc2.", gsub(" ", ".", agg.reward.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.reward.diff.pc2[i, grep("alpha", names(agg.reward.diff.pc2))])), ncol = length(states.reward), nrow = length(states.reward), byrow = TRUE), row.names = states.reward))
    assign(paste0("reward.diff.sigma.pc2.", gsub(" ", ".", agg.reward.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.reward.diff.pc2[i, grep("sigma", names(agg.reward.diff.pc2))])), ncol = length(states.reward), nrow = length(states.reward), byrow = TRUE), row.names = states.reward))
}


### PC3

theta.reward.plot.pc3 <- subset(theta.plot.pc3, disc.trait == "reward")
theta.reward.plot.pc3$variable <- as.character(theta.reward.plot.pc3$variable)
theta.reward.plot.pc3 <- theta.reward.plot.pc3[theta.reward.plot.pc3$value >= -10 & theta.reward.plot.pc3$value <= 10, ]

alpha.reward.plot.pc3 <- subset(alpha.plot.pc3, disc.trait == "reward")
alpha.reward.plot.pc3$variable <- as.character(alpha.reward.plot.pc3$variable)
alpha.reward.plot.pc3 <- alpha.reward.plot.pc3[alpha.reward.plot.pc3$value >= -10 & alpha.reward.plot.pc3$value <= 10, ]

sigma.reward.plot.pc3 <- subset(sigma.plot.pc3, disc.trait == "reward")
sigma.reward.plot.pc3$variable <- as.character(sigma.reward.plot.pc3$variable)
sigma.reward.plot.pc3 <- sigma.reward.plot.pc3[sigma.reward.plot.pc3$value >= -10 & sigma.reward.plot.pc3$value <= 10, ]

reward.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "reward", -c(6, 10, 14, 15)]
reward.diff.pc3 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(reward.plot.diff.pc3[, grep("theta", names(reward.plot.diff.pc3))] >= reward.plot.diff.pc3$theta1))


#### Differences
reward.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "reward", -c(6, 10, 14, 15)]

for(i in 1:length(states.reward)){
    for(j in 1:length(states.reward)){
        reward.plot.diff.pc3[, paste0("theta", i, j)] <- reward.plot.diff.pc3[, paste0("theta", i)] > reward.plot.diff.pc3[, paste0("theta", j)]
        reward.plot.diff.pc3[, paste0("alpha", i, j)] <- reward.plot.diff.pc3[, paste0("alpha", i)] > reward.plot.diff.pc3[, paste0("alpha", j)]
        reward.plot.diff.pc3[, paste0("sigma", i, j)] <- reward.plot.diff.pc3[, paste0("sigma", i)] > reward.plot.diff.pc3[, paste0("sigma", j)]
    }
}

agg.reward.diff.pc3 <- aggregate(reward.plot.diff.pc3[reward.plot.diff.pc3$trait == "Hole Diameter", 12:ncol(reward.plot.diff.pc3)], by = list(reward.plot.diff.pc3$trait[reward.plot.diff.pc3$trait == "Hole Diameter"], reward.plot.diff.pc3$disc.trait[reward.plot.diff.pc3$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.reward.diff.pc3$Group.1)){
    assign(paste0("reward.diff.theta.pc3.", gsub(" ", ".", agg.reward.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.reward.diff.pc3[i, grep("theta", names(agg.reward.diff.pc3))])), ncol = length(states.reward), nrow = length(states.reward), byrow = TRUE), row.names = states.reward))
    assign(paste0("reward.diff.alpha.pc3.", gsub(" ", ".", agg.reward.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.reward.diff.pc3[i, grep("alpha", names(agg.reward.diff.pc3))])), ncol = length(states.reward), nrow = length(states.reward), byrow = TRUE), row.names = states.reward))
    assign(paste0("reward.diff.sigma.pc3.", gsub(" ", ".", agg.reward.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.reward.diff.pc3[i, grep("sigma", names(agg.reward.diff.pc3))])), ncol = length(states.reward), nrow = length(states.reward), byrow = TRUE), row.names = states.reward))
}







## Strategy
### PC1

theta.strategy.plot.pc1 <- subset(theta.plot.pc1, disc.trait == "strategy")
theta.strategy.plot.pc1$variable <- as.character(theta.strategy.plot.pc1$variable)
theta.strategy.plot.pc1 <- theta.strategy.plot.pc1[theta.strategy.plot.pc1$value >= -10 & theta.strategy.plot.pc1$value <= 10, ]

alpha.strategy.plot.pc1 <- subset(alpha.plot.pc1, disc.trait == "strategy")
alpha.strategy.plot.pc1$variable <- as.character(alpha.strategy.plot.pc1$variable)
alpha.strategy.plot.pc1 <- alpha.strategy.plot.pc1[alpha.strategy.plot.pc1$value >= -10 & alpha.strategy.plot.pc1$value <= 10, ]

sigma.strategy.plot.pc1 <- subset(sigma.plot.pc1, disc.trait == "strategy")
sigma.strategy.plot.pc1$variable <- as.character(sigma.strategy.plot.pc1$variable)
sigma.strategy.plot.pc1 <- sigma.strategy.plot.pc1[sigma.strategy.plot.pc1$value >= -10 & sigma.strategy.plot.pc1$value <= 10, ]


#### Differences
strategy.plot.diff.pc1 <- pars.table.pc1[pars.table.pc1$disc.trait == "strategy", -c(6, 10, 14, 15)]

for(i in 1:length(states.strategy)){
    for(j in 1:length(states.strategy)){
        strategy.plot.diff.pc1[, paste0("theta", i, j)] <- strategy.plot.diff.pc1[, paste0("theta", i)] > strategy.plot.diff.pc1[, paste0("theta", j)]
        strategy.plot.diff.pc1[, paste0("alpha", i, j)] <- strategy.plot.diff.pc1[, paste0("alpha", i)] > strategy.plot.diff.pc1[, paste0("alpha", j)]
        strategy.plot.diff.pc1[, paste0("sigma", i, j)] <- strategy.plot.diff.pc1[, paste0("sigma", i)] > strategy.plot.diff.pc1[, paste0("sigma", j)]
    }
}

agg.strategy.diff.pc1 <- aggregate(strategy.plot.diff.pc1[strategy.plot.diff.pc1$trait == "Hole Diameter", 12:ncol(strategy.plot.diff.pc1)], by = list(strategy.plot.diff.pc1$trait[strategy.plot.diff.pc1$trait == "Hole Diameter"], strategy.plot.diff.pc1$disc.trait[strategy.plot.diff.pc1$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.strategy.diff.pc1$Group.1)){
    assign(paste0("strategy.diff.theta.pc1.", gsub(" ", ".", agg.strategy.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.strategy.diff.pc1[i, grep("theta", names(agg.strategy.diff.pc1))])), ncol = length(states.strategy), nrow = length(states.strategy), byrow = TRUE), row.names = states.strategy))
    assign(paste0("strategy.diff.alpha.pc1.", gsub(" ", ".", agg.strategy.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.strategy.diff.pc1[i, grep("alpha", names(agg.strategy.diff.pc1))])), ncol = length(states.strategy), nrow = length(states.strategy), byrow = TRUE), row.names = states.strategy))
    assign(paste0("strategy.diff.sigma.pc1.", gsub(" ", ".", agg.strategy.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.strategy.diff.pc1[i, grep("sigma", names(agg.strategy.diff.pc1))])), ncol = length(states.strategy), nrow = length(states.strategy), byrow = TRUE), row.names = states.strategy))
}

### PC2

theta.strategy.plot.pc2 <- subset(theta.plot.pc2, disc.trait == "strategy")
theta.strategy.plot.pc2$variable <- as.character(theta.strategy.plot.pc2$variable)
theta.strategy.plot.pc2 <- theta.strategy.plot.pc2[theta.strategy.plot.pc2$value >= -10 & theta.strategy.plot.pc2$value <= 10, ]

alpha.strategy.plot.pc2 <- subset(alpha.plot.pc2, disc.trait == "strategy")
alpha.strategy.plot.pc2$variable <- as.character(alpha.strategy.plot.pc2$variable)
alpha.strategy.plot.pc2 <- alpha.strategy.plot.pc2[alpha.strategy.plot.pc2$value >= -10 & alpha.strategy.plot.pc2$value <= 10, ]

sigma.strategy.plot.pc2 <- subset(sigma.plot.pc2, disc.trait == "strategy")
sigma.strategy.plot.pc2$variable <- as.character(sigma.strategy.plot.pc2$variable)
sigma.strategy.plot.pc2 <- sigma.strategy.plot.pc2[sigma.strategy.plot.pc2$value >= -10 & sigma.strategy.plot.pc2$value <= 10, ]

strategy.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "strategy", -c(6, 10, 14, 15)]
strategy.diff.pc2 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(strategy.plot.diff.pc2[, grep("theta", names(strategy.plot.diff.pc2))] >= strategy.plot.diff.pc2$theta1))


#### Differences
strategy.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "strategy", -c(6, 10, 14, 15)]

for(i in 1:length(states.strategy)){
    for(j in 1:length(states.strategy)){
        strategy.plot.diff.pc2[, paste0("theta", i, j)] <- strategy.plot.diff.pc2[, paste0("theta", i)] > strategy.plot.diff.pc2[, paste0("theta", j)]
        strategy.plot.diff.pc2[, paste0("alpha", i, j)] <- strategy.plot.diff.pc2[, paste0("alpha", i)] > strategy.plot.diff.pc2[, paste0("alpha", j)]
        strategy.plot.diff.pc2[, paste0("sigma", i, j)] <- strategy.plot.diff.pc2[, paste0("sigma", i)] > strategy.plot.diff.pc2[, paste0("sigma", j)]
    }
}

agg.strategy.diff.pc2 <- aggregate(strategy.plot.diff.pc2[strategy.plot.diff.pc2$trait == "Hole Diameter", 12:ncol(strategy.plot.diff.pc2)], by = list(strategy.plot.diff.pc2$trait[strategy.plot.diff.pc2$trait == "Hole Diameter"], strategy.plot.diff.pc2$disc.trait[strategy.plot.diff.pc2$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.strategy.diff.pc2$Group.1)){
    assign(paste0("strategy.diff.theta.pc2.", gsub(" ", ".", agg.strategy.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.strategy.diff.pc2[i, grep("theta", names(agg.strategy.diff.pc2))])), ncol = length(states.strategy), nrow = length(states.strategy), byrow = TRUE), row.names = states.strategy))
    assign(paste0("strategy.diff.alpha.pc2.", gsub(" ", ".", agg.strategy.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.strategy.diff.pc2[i, grep("alpha", names(agg.strategy.diff.pc2))])), ncol = length(states.strategy), nrow = length(states.strategy), byrow = TRUE), row.names = states.strategy))
    assign(paste0("strategy.diff.sigma.pc2.", gsub(" ", ".", agg.strategy.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.strategy.diff.pc2[i, grep("sigma", names(agg.strategy.diff.pc2))])), ncol = length(states.strategy), nrow = length(states.strategy), byrow = TRUE), row.names = states.strategy))
}


### PC3

theta.strategy.plot.pc3 <- subset(theta.plot.pc3, disc.trait == "strategy")
theta.strategy.plot.pc3$variable <- as.character(theta.strategy.plot.pc3$variable)
theta.strategy.plot.pc3 <- theta.strategy.plot.pc3[theta.strategy.plot.pc3$value >= -10 & theta.strategy.plot.pc3$value <= 10, ]

alpha.strategy.plot.pc3 <- subset(alpha.plot.pc3, disc.trait == "strategy")
alpha.strategy.plot.pc3$variable <- as.character(alpha.strategy.plot.pc3$variable)
alpha.strategy.plot.pc3 <- alpha.strategy.plot.pc3[alpha.strategy.plot.pc3$value >= -10 & alpha.strategy.plot.pc3$value <= 10, ]

sigma.strategy.plot.pc3 <- subset(sigma.plot.pc3, disc.trait == "strategy")
sigma.strategy.plot.pc3$variable <- as.character(sigma.strategy.plot.pc3$variable)
sigma.strategy.plot.pc3 <- sigma.strategy.plot.pc3[sigma.strategy.plot.pc3$value >= -10 & sigma.strategy.plot.pc3$value <= 10, ]

strategy.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "strategy", -c(6, 10, 14, 15)]
strategy.diff.pc3 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(strategy.plot.diff.pc3[, grep("theta", names(strategy.plot.diff.pc3))] >= strategy.plot.diff.pc3$theta1))


#### Differences
strategy.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "strategy", -c(6, 10, 14, 15)]

for(i in 1:length(states.strategy)){
    for(j in 1:length(states.strategy)){
        strategy.plot.diff.pc3[, paste0("theta", i, j)] <- strategy.plot.diff.pc3[, paste0("theta", i)] > strategy.plot.diff.pc3[, paste0("theta", j)]
        strategy.plot.diff.pc3[, paste0("alpha", i, j)] <- strategy.plot.diff.pc3[, paste0("alpha", i)] > strategy.plot.diff.pc3[, paste0("alpha", j)]
        strategy.plot.diff.pc3[, paste0("sigma", i, j)] <- strategy.plot.diff.pc3[, paste0("sigma", i)] > strategy.plot.diff.pc3[, paste0("sigma", j)]
    }
}

agg.strategy.diff.pc3 <- aggregate(strategy.plot.diff.pc3[strategy.plot.diff.pc3$trait == "Hole Diameter", 12:ncol(strategy.plot.diff.pc3)], by = list(strategy.plot.diff.pc3$trait[strategy.plot.diff.pc3$trait == "Hole Diameter"], strategy.plot.diff.pc3$disc.trait[strategy.plot.diff.pc3$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.strategy.diff.pc3$Group.1)){
    assign(paste0("strategy.diff.theta.pc3.", gsub(" ", ".", agg.strategy.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.strategy.diff.pc3[i, grep("theta", names(agg.strategy.diff.pc3))])), ncol = length(states.strategy), nrow = length(states.strategy), byrow = TRUE), row.names = states.strategy))
    assign(paste0("strategy.diff.alpha.pc3.", gsub(" ", ".", agg.strategy.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.strategy.diff.pc3[i, grep("alpha", names(agg.strategy.diff.pc3))])), ncol = length(states.strategy), nrow = length(states.strategy), byrow = TRUE), row.names = states.strategy))
    assign(paste0("strategy.diff.sigma.pc3.", gsub(" ", ".", agg.strategy.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.strategy.diff.pc3[i, grep("sigma", names(agg.strategy.diff.pc3))])), ncol = length(states.strategy), nrow = length(states.strategy), byrow = TRUE), row.names = states.strategy))
}






## Warts
### PC1

theta.warts.plot.pc1 <- subset(theta.plot.pc1, disc.trait == "warts")
theta.warts.plot.pc1$variable <- as.character(theta.warts.plot.pc1$variable)
theta.warts.plot.pc1 <- theta.warts.plot.pc1[theta.warts.plot.pc1$value >= -10 & theta.warts.plot.pc1$value <= 10, ]

alpha.warts.plot.pc1 <- subset(alpha.plot.pc1, disc.trait == "warts")
alpha.warts.plot.pc1$variable <- as.character(alpha.warts.plot.pc1$variable)
alpha.warts.plot.pc1 <- alpha.warts.plot.pc1[alpha.warts.plot.pc1$value >= -10 & alpha.warts.plot.pc1$value <= 10, ]

sigma.warts.plot.pc1 <- subset(sigma.plot.pc1, disc.trait == "warts")
sigma.warts.plot.pc1$variable <- as.character(sigma.warts.plot.pc1$variable)
sigma.warts.plot.pc1 <- sigma.warts.plot.pc1[sigma.warts.plot.pc1$value >= -10 & sigma.warts.plot.pc1$value <= 10, ]


#### Differences
warts.plot.diff.pc1 <- pars.table.pc1[pars.table.pc1$disc.trait == "warts", -c(6, 10, 14, 15)]

for(i in 1:length(states.warts)){
    for(j in 1:length(states.warts)){
        warts.plot.diff.pc1[, paste0("theta", i, j)] <- warts.plot.diff.pc1[, paste0("theta", i)] > warts.plot.diff.pc1[, paste0("theta", j)]
        warts.plot.diff.pc1[, paste0("alpha", i, j)] <- warts.plot.diff.pc1[, paste0("alpha", i)] > warts.plot.diff.pc1[, paste0("alpha", j)]
        warts.plot.diff.pc1[, paste0("sigma", i, j)] <- warts.plot.diff.pc1[, paste0("sigma", i)] > warts.plot.diff.pc1[, paste0("sigma", j)]
    }
}

agg.warts.diff.pc1 <- aggregate(warts.plot.diff.pc1[warts.plot.diff.pc1$trait == "Hole Diameter", 12:ncol(warts.plot.diff.pc1)], by = list(warts.plot.diff.pc1$trait[warts.plot.diff.pc1$trait == "Hole Diameter"], warts.plot.diff.pc1$disc.trait[warts.plot.diff.pc1$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.warts.diff.pc1$Group.1)){
    assign(paste0("warts.diff.theta.pc1.", gsub(" ", ".", agg.warts.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.warts.diff.pc1[i, grep("theta", names(agg.warts.diff.pc1))])), ncol = length(states.warts), nrow = length(states.warts), byrow = TRUE), row.names = states.warts))
    assign(paste0("warts.diff.alpha.pc1.", gsub(" ", ".", agg.warts.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.warts.diff.pc1[i, grep("alpha", names(agg.warts.diff.pc1))])), ncol = length(states.warts), nrow = length(states.warts), byrow = TRUE), row.names = states.warts))
    assign(paste0("warts.diff.sigma.pc1.", gsub(" ", ".", agg.warts.diff.pc1$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.warts.diff.pc1[i, grep("sigma", names(agg.warts.diff.pc1))])), ncol = length(states.warts), nrow = length(states.warts), byrow = TRUE), row.names = states.warts))
}

### PC2

theta.warts.plot.pc2 <- subset(theta.plot.pc2, disc.trait == "warts")
theta.warts.plot.pc2$variable <- as.character(theta.warts.plot.pc2$variable)
theta.warts.plot.pc2 <- theta.warts.plot.pc2[theta.warts.plot.pc2$value >= -10 & theta.warts.plot.pc2$value <= 10, ]

alpha.warts.plot.pc2 <- subset(alpha.plot.pc2, disc.trait == "warts")
alpha.warts.plot.pc2$variable <- as.character(alpha.warts.plot.pc2$variable)
alpha.warts.plot.pc2 <- alpha.warts.plot.pc2[alpha.warts.plot.pc2$value >= -10 & alpha.warts.plot.pc2$value <= 10, ]

sigma.warts.plot.pc2 <- subset(sigma.plot.pc2, disc.trait == "warts")
sigma.warts.plot.pc2$variable <- as.character(sigma.warts.plot.pc2$variable)
sigma.warts.plot.pc2 <- sigma.warts.plot.pc2[sigma.warts.plot.pc2$value >= -10 & sigma.warts.plot.pc2$value <= 10, ]

warts.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "warts", -c(6, 10, 14, 15)]
warts.diff.pc2 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(warts.plot.diff.pc2[, grep("theta", names(warts.plot.diff.pc2))] >= warts.plot.diff.pc2$theta1))


#### Differences
warts.plot.diff.pc2 <- pars.table.pc2[pars.table.pc2$disc.trait == "warts", -c(6, 10, 14, 15)]

for(i in 1:length(states.warts)){
    for(j in 1:length(states.warts)){
        warts.plot.diff.pc2[, paste0("theta", i, j)] <- warts.plot.diff.pc2[, paste0("theta", i)] > warts.plot.diff.pc2[, paste0("theta", j)]
        warts.plot.diff.pc2[, paste0("alpha", i, j)] <- warts.plot.diff.pc2[, paste0("alpha", i)] > warts.plot.diff.pc2[, paste0("alpha", j)]
        warts.plot.diff.pc2[, paste0("sigma", i, j)] <- warts.plot.diff.pc2[, paste0("sigma", i)] > warts.plot.diff.pc2[, paste0("sigma", j)]
    }
}

agg.warts.diff.pc2 <- aggregate(warts.plot.diff.pc2[warts.plot.diff.pc2$trait == "Hole Diameter", 12:ncol(warts.plot.diff.pc2)], by = list(warts.plot.diff.pc2$trait[warts.plot.diff.pc2$trait == "Hole Diameter"], warts.plot.diff.pc2$disc.trait[warts.plot.diff.pc2$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.warts.diff.pc2$Group.1)){
    assign(paste0("warts.diff.theta.pc2.", gsub(" ", ".", agg.warts.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.warts.diff.pc2[i, grep("theta", names(agg.warts.diff.pc2))])), ncol = length(states.warts), nrow = length(states.warts), byrow = TRUE), row.names = states.warts))
    assign(paste0("warts.diff.alpha.pc2.", gsub(" ", ".", agg.warts.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.warts.diff.pc2[i, grep("alpha", names(agg.warts.diff.pc2))])), ncol = length(states.warts), nrow = length(states.warts), byrow = TRUE), row.names = states.warts))
    assign(paste0("warts.diff.sigma.pc2.", gsub(" ", ".", agg.warts.diff.pc2$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.warts.diff.pc2[i, grep("sigma", names(agg.warts.diff.pc2))])), ncol = length(states.warts), nrow = length(states.warts), byrow = TRUE), row.names = states.warts))
}


### PC3

theta.warts.plot.pc3 <- subset(theta.plot.pc3, disc.trait == "warts")
theta.warts.plot.pc3$variable <- as.character(theta.warts.plot.pc3$variable)
theta.warts.plot.pc3 <- theta.warts.plot.pc3[theta.warts.plot.pc3$value >= -10 & theta.warts.plot.pc3$value <= 10, ]

alpha.warts.plot.pc3 <- subset(alpha.plot.pc3, disc.trait == "warts")
alpha.warts.plot.pc3$variable <- as.character(alpha.warts.plot.pc3$variable)
alpha.warts.plot.pc3 <- alpha.warts.plot.pc3[alpha.warts.plot.pc3$value >= -10 & alpha.warts.plot.pc3$value <= 10, ]

sigma.warts.plot.pc3 <- subset(sigma.plot.pc3, disc.trait == "warts")
sigma.warts.plot.pc3$variable <- as.character(sigma.warts.plot.pc3$variable)
sigma.warts.plot.pc3 <- sigma.warts.plot.pc3[sigma.warts.plot.pc3$value >= -10 & sigma.warts.plot.pc3$value <= 10, ]

warts.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "warts", -c(6, 10, 14, 15)]
warts.diff.pc3 <- data.frame(state = c("None", "Spines", "Variable"),
                           none = sum(warts.plot.diff.pc3[, grep("theta", names(warts.plot.diff.pc3))] >= warts.plot.diff.pc3$theta1))


#### Differences
warts.plot.diff.pc3 <- pars.table.pc3[pars.table.pc3$disc.trait == "warts", -c(6, 10, 14, 15)]

for(i in 1:length(states.warts)){
    for(j in 1:length(states.warts)){
        warts.plot.diff.pc3[, paste0("theta", i, j)] <- warts.plot.diff.pc3[, paste0("theta", i)] > warts.plot.diff.pc3[, paste0("theta", j)]
        warts.plot.diff.pc3[, paste0("alpha", i, j)] <- warts.plot.diff.pc3[, paste0("alpha", i)] > warts.plot.diff.pc3[, paste0("alpha", j)]
        warts.plot.diff.pc3[, paste0("sigma", i, j)] <- warts.plot.diff.pc3[, paste0("sigma", i)] > warts.plot.diff.pc3[, paste0("sigma", j)]
    }
}

agg.warts.diff.pc3 <- aggregate(warts.plot.diff.pc3[warts.plot.diff.pc3$trait == "Hole Diameter", 12:ncol(warts.plot.diff.pc3)], by = list(warts.plot.diff.pc3$trait[warts.plot.diff.pc3$trait == "Hole Diameter"], warts.plot.diff.pc3$disc.trait[warts.plot.diff.pc3$trait == "Hole Diameter"]), FUN = sum)

for(i in 1:length(agg.warts.diff.pc3$Group.1)){
    assign(paste0("warts.diff.theta.pc3.", gsub(" ", ".", agg.warts.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.warts.diff.pc3[i, grep("theta", names(agg.warts.diff.pc3))])), ncol = length(states.warts), nrow = length(states.warts), byrow = TRUE), row.names = states.warts))
    assign(paste0("warts.diff.alpha.pc3.", gsub(" ", ".", agg.warts.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.warts.diff.pc3[i, grep("alpha", names(agg.warts.diff.pc3))])), ncol = length(states.warts), nrow = length(states.warts), byrow = TRUE), row.names = states.warts))
    assign(paste0("warts.diff.sigma.pc3.", gsub(" ", ".", agg.warts.diff.pc3$Group.1[i])), as.data.frame(matrix(as.numeric(unlist(agg.warts.diff.pc3[i, grep("sigma", names(agg.warts.diff.pc3))])), ncol = length(states.warts), nrow = length(states.warts), byrow = TRUE), row.names = states.warts))
}

