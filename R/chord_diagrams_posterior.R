here::i_am("R/chord_diagrams_posterior.R")
source(here::here("R/data_prep_ouwie_climpc_posterior.R"))
library("circlize")
library("MetBrewer")
library("ggplot2")

## Calculating differences in relation to state 1

pars.table.pc1.clean <- pars.table.pc1[apply(pars.table.pc1[, grep("theta", names(pars.table.pc1))], 1, function(x){!any(x > 15, na.rm = TRUE)}),]
pars.table.pc1.clean <- pars.table.pc1.clean[apply(pars.table.pc1.clean[, grep("theta", names(pars.table.pc1.clean))], 1, function(x){!any(x < -15, na.rm = TRUE)}),]

reldiff.pc1 <- data.frame(trait = pars.table.pc1.clean$trait,
                          disc.trait = pars.table.pc1.clean$disc.trait,
                          alpha1 = round(pars.table.pc1.clean$alpha1 - pars.table.pc1.clean$alpha1, 4),
                          alpha2 = round(pars.table.pc1.clean$alpha2 - pars.table.pc1.clean$alpha1, 4),
                          alpha3 = round(pars.table.pc1.clean$alpha3 - pars.table.pc1.clean$alpha1, 4),
                          alpha4 = round(pars.table.pc1.clean$alpha4 - pars.table.pc1.clean$alpha1, 4),
                          sigma1 = round(pars.table.pc1.clean$sigma1 - pars.table.pc1.clean$sigma1, 4),
                          sigma2 = round(pars.table.pc1.clean$sigma2 - pars.table.pc1.clean$sigma1, 4),
                          sigma3 = round(pars.table.pc1.clean$sigma3 - pars.table.pc1.clean$sigma1, 4),
                          sigma4 = round(pars.table.pc1.clean$sigma4 - pars.table.pc1.clean$sigma1, 4),
                          theta1 = round(pars.table.pc1.clean$theta1 - pars.table.pc1.clean$theta1, 4),
                          theta2 = round(pars.table.pc1.clean$theta2 - pars.table.pc1.clean$theta1, 4),
                          theta3 = round(pars.table.pc1.clean$theta3 - pars.table.pc1.clean$theta1, 4),
                          theta4 = round(pars.table.pc1.clean$theta4 - pars.table.pc1.clean$theta1, 4),
                          simmap = pars.table.pc1.clean$simmap)

reldiff.pc1.app <- subset(reldiff.pc1, disc.trait == "app") ## 385 out of 500
ranks.pc1.app <- apply(reldiff.pc1.app[, grep("theta", names(reldiff.pc1.app))], 1, order)
## aggregate(ranks.pc1.app, 1, mean)

reldiff.pc1.arch <- subset(reldiff.pc1, disc.trait == "arch") ## 411 out of 500
reldiff.pc1.domgrow <- subset(reldiff.pc1, disc.trait == "domgrow") ## 450 out of 500
reldiff.pc1.leafstruc <- subset(reldiff.pc1, disc.trait == "leafstruc") ## 200 out of 500
reldiff.pc1.matsys <- subset(reldiff.pc1, disc.trait == "matsys") ## 251 out of 500
reldiff.pc1.reward <- subset(reldiff.pc1, disc.trait == "reward") ## 403 out of 500
reldiff.pc1.strategy <- subset(reldiff.pc1, disc.trait == "strategy") ## 437 out of 500
reldiff.pc1.warts <- subset(reldiff.pc1, disc.trait == "warts") ## 441 out of 500
reldiff.pc1.holediam.disc <- subset(reldiff.pc1, disc.trait == "holediam.disc") ## 400 out of 500


ouwie.chord <- function(trait, disc.trait, par, pc, pal){
    raw.data <- get(paste(disc.trait, "diff", par, paste0("pc", pc), trait, sep = "."))
    data <- as.matrix(raw.data)
    rownames(data) <- rownames(raw.data)
    colnames(data) <- rownames(raw.data)
    cols <- colorRampPalette(met.brewer(name = pal, n = 100, type = "continuous"))(length(rownames(data)))
    chordDiagram(data,
                 preAllocateTracks = list(grid.col = setNames(cols, rownames(data))),
                 row.col = cols,
                 grid.col = setNames(cols, rownames(data)),
                 annotationTrack = c("grid"),
                 annotationTrackHeight = mm_h(5),
                 order = rownames(data)
                 )
    circos.clear()
}

## disc.traits <- unique(theta.plot.pc1$disc.trait)
disc.traits <- c("domgrow", "reward", "warts", "holediam.disc")
pars <- c("theta", "alpha", "sigma")
traits <- gsub(" ", ".", unique(theta.plot.pc1$trait))
nstates <- c(3, 2, 2, 4, 3, 2, 3, 3, 3)
pals <- c(#"Egypt",
          #"Derain",
          #"Pissaro",
          "Greek",
          #"VanGogh2",
          #"Troy",
          #"Peru1",
          #"Redon",
          #"Monet",
          "Hiroshige",
          #"Hokusai3",
          "Isfahan2",
          "Veronese"
          #"Morgenstern"
          )


pdf(here::here("output/chord_diags_diffs_posterior_new.pdf"), height = 40, width = 40)
par(mfrow = c(length(disc.traits), length(traits) * length(pars)), mar = c(0, 0, 0, 0))
for(i in 1:length(disc.traits)){
    for(j in 1:length(traits)){
        for(k in 1:length(pars)){
            circos.clear()
            tryCatch(ouwie.chord(trait = traits[j], disc.trait = disc.traits[i], par = pars[k], pc = 1, pal = pals[i]), error = function(x){plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10))})
        }
    }
}
dev.off()
