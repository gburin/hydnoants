library("ape")
library("tidyverse")
library("cowplot")
library("ggtree")
library("ggforce")
library("ggnewscale")
library("patchwork")
library("ggsci")
library("grafify")
library("wesanderson")
here::i_am("R/tree_plots.R")
trees <- read.nexus(here::here("data/posterior_trimmed_135taxa.tre"))
tree <- ladderize(trees[[300]])

tiplabel.fix <- read.csv(here::here("data/tiplabel_fix.csv"))

tree$tip.label[match(gsub(" ", "_", tiplabel.fix[, "phylo"]), tree$tip.label)] <- gsub(" ", "_", tiplabel.fix[, "fix"])

traits <- read.csv(here::here("data/Dataset_traits_Hydnophytinae_CODED_2022.csv"), as.is = TRUE)
names(traits) <- c("species", paste0("bio", 1:19), "strategy", "warts", "stem.area", "leaf.area", "hole.diameter", "architecture", "dom.growth", "reward", "corola.length", "mating.system", "petiole.length", "leaf.struct", "appendages", "holediam.disc")

traits$species[grep("Wor2g", traits$species)] <- gsub("Wor2g", "Worthing", traits$species[grep("Wor2g", traits$species)])

traits$hole.diameter[which(traits$hole.diameter == "NA (outgroup has no domatia)")] <- NA
traits$hole.diameter <- as.numeric(traits$hole.diameter)

strategy <- c("outgroup", "facultative", "obligate", "lost")
warts <- c("absent", "variable", "differentiated", "lost")
reward <- c("absent", "present")
archit <- c("shrub", "multiple", "single")
dom.growth <- c("outgroup", "diffuse", "apical")
mating.system <- c("heterostylous", "non-heterostylous", "funct_unisexual")
leaf.struct <- c("thick", "variable", "thin", "succulent")
appendages <- c("outgroup", "none", "variable", "spines")
holediam.disc <- c("absent", "oneLarge", "sevLargeBase", "allLarge")

traits$strategy <- strategy[traits$strategy + 1]
traits$warts <- warts[traits$warts + 1]
traits$reward <- reward[traits$reward + 1]
traits$architecture <- archit[traits$architecture + 1]
traits$dom.growth <- dom.growth[traits$dom.growth + 1]
traits$mating.system <- mating.system[traits$mating.system + 1]
traits$leaf.struct <- leaf.struct[traits$leaf.struct + 1]
traits$appendages <- appendages[traits$appendages + 1]
traits$holediam.disc <- holediam.disc[traits$holediam.disc + 1]

## 1, 2, 6, 7

traits.plot <- traits[traits$strategy != "outgroup", -1]
rownames(traits.plot) <- traits[traits$strategy != "outgroup", 1]

tree.full <- tree
tree <- drop.tip(tree, traits$species[traits$strategy == "outgroup"])

## Fan layout

## p <- ggtree(tree, layout = "fan", open.angle = 30, xlim = c(0, 35)) %<+% traits

## ## p1 <-
## ##     p +
## ##     #geom_tiplab(offset = 1, hjust = 0) +
## ##     theme(legend.position = "bottom") +
## ##     scale_shape_manual(values = c(15, 16, 17, 18)) +
## ##     scale_colour_brewer(palette = "Set1")

## col1 <- setNames((grafify:::graf_palettes)$okabe_ito[c(1, 3, 4, 2)], NULL)

## p1 <-
##     gheatmap(p = p, data = traits.plot[, c("dom.growth", "reward")], width = 0.25, colnames_angle = 45, hjust = 1) +
##     geom_tiplab(size = 3, offset = 23) +
##     theme_tree2() +
##     scale_fill_manual(values = alpha(col1, 0.75), name = "DomGrowth")

## p1.2 <- p1 + new_scale_fill()

## ## col2 <- setNames((grafify:::graf_palettes)$okabe_ito[c(1, 3, 6, 7, 8, 4)], NULL)
## col2 <- setNames((grafify:::graf_palettes)$muted[c(1:6)], NULL)

## p1.3 <-
##     gheatmap(p = p1.2, data = traits.plot[, c("warts", "holediam.disc")], offset = 4.5, width = 0.25, colnames_angle = 45, hjust = 1) +
##     #scale_fill_brewer(palette = "Set1", name = "Reward")
##     scale_fill_manual(values = alpha(col2, 0.75))

## p1.4 <- p1.3 + new_scale_fill()

## ## p1.5 <-
## ##     gheatmap(p = p1.4, data = traits.plot[, c("warts")], offset = 12, width = 0.1, colnames_angle = 45, hjust = 1) +
## ##     scale_fill_brewer(palette = "Set3", direction = -1, name = "Warts")

## ## p1.6 <- p1.5 + new_scale_fill()

## ## p1.7 <-
## ##     gheatmap(p = p1.6, data = traits.plot[, c("holediam.disc")], offset = 18, width = 0.1, colnames_angle = 45, hjust = 1) +
## ##     scale_fill_brewer(palette = "Dark2", direction = 1, name = "HoleDiam(disc)")

## ## p1.8 <- p1.7 + new_scale_fill()

## p2 <-
##     gheatmap(p = p1.4, data = traits.plot[, c("stem.area", "leaf.area", "corola.length", "petiole.length")], offset = 12, width = 0.5, colnames_angle = 45, hjust = 1) +
##     scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral"), name = "Continuous Traits") +
##     theme(legend.position = "left")

## #p2

## ## p2.1 <- p2 + new_scale_fill()

## ## p3 <-
##     ## gheatmap(p = p2.1, data = traits.plot[, c("hole.diameter", "corola.length", "petiole.length")], offset = 30, width = 0.1, colnames_angle = 45, hjust = 1) +
##     ## theme(legend.position = "bottom") +
##     ## scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral"), name = "Hole Diameter, Corola\nand Petiole Length")


## ggsave(filename = here::here("output/preliminary_phylo_fan_4traits.pdf"), plot = p2, width = 30, height = 30)




## ## Normal layout

## p <- ggtree(tree) %<+% traits

## #p <- ggtree(tree, layout = "fan", open.angle = 30) %<+% traits

## ## p1 <-
## ##     p +
## ##     #geom_tiplab(offset = 1, hjust = 0) +
## ##     theme(legend.position = "bottom") +
## ##     scale_shape_manual(values = c(15, 16, 17, 18)) +
## ##     scale_colour_brewer(palette = "Set1")

## col1 <- setNames((grafify:::graf_palettes)$okabe_ito[c(1, 3, 4, 2)], NULL)

## p1 <-
##     gheatmap(p = p, data = traits.plot[, c("dom.growth", "reward")], width = 0.1, colnames_angle = 45, hjust = 1) +
##     geom_tiplab(size = 2, offset = 9) +
##     theme_tree2() +
##     scale_fill_manual(values = alpha(col1, 0.75), name = "DomGrowth") +
##     xlim(0, 35) +
##     ylim(-5, 90)

## p1.2 <- p1 + new_scale_fill()

## ## col2 <- setNames((grafify:::graf_palettes)$okabe_ito[c(1, 3, 6, 7, 8, 4)], NULL)
## col2 <- setNames((grafify:::graf_palettes)$muted[c(1:6)], NULL)

## p1.3 <-
##     gheatmap(p = p1.2, data = traits.plot[, c("warts", "holediam.disc")], offset = 1.85, width = 0.1, colnames_angle = 45, hjust = 1) +
##     #scale_fill_brewer(palette = "Set1", name = "Reward")
##     scale_fill_manual(values = alpha(col2, 0.75))

## p1.4 <- p1.3 + new_scale_fill()

## ## p1.5 <-
## ##     gheatmap(p = p1.4, data = traits.plot[, c("warts")], offset = 12, width = 0.1, colnames_angle = 45, hjust = 1) +
## ##     scale_fill_brewer(palette = "Set3", direction = -1, name = "Warts")

## ## p1.6 <- p1.5 + new_scale_fill()

## ## p1.7 <-
## ##     gheatmap(p = p1.6, data = traits.plot[, c("holediam.disc")], offset = 18, width = 0.1, colnames_angle = 45, hjust = 1) +
## ##     scale_fill_brewer(palette = "Dark2", direction = 1, name = "HoleDiam(disc)")

## ## p1.8 <- p1.7 + new_scale_fill()

## p2 <-
##     gheatmap(p = p1.4, data = traits.plot[, c("stem.area", "leaf.area", "corola.length", "petiole.length")], offset = 4.5, width = 0.2, colnames_angle = 45, hjust = 1) +
##     scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral"), name = "Continuous Traits") +
##     theme(legend.position = "left")

## #p2

## ## p2.1 <- p2 + new_scale_fill()

## ## p3 <-
##     ## gheatmap(p = p2.1, data = traits.plot[, c("hole.diameter", "corola.length", "petiole.length")], offset = 30, width = 0.1, colnames_angle = 45, hjust = 1) +
##     ## theme(legend.position = "bottom") +
##     ## scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral"), name = "Hole Diameter, Corola\nand Petiole Length")


## ggsave(filename = here::here("output/preliminary_phylo_regular_4traits.pdf"), plot = p2, width = 20, height = 20)







## Including outgroup

traits.plot <- traits[, -1]
rownames(traits.plot) <- traits[, 1]
rownames(traits.plot) <- gsub("P_", "Psychotria_", rownames(traits.plot))

tree.full <- trees[[300]]

tiplabel.fix <- read.csv(here::here("data/tiplabel_fix.csv"))
tiplabel.fix <- rbind(tiplabel.fix, data.frame(phylo = tree.full$tip.label[grep("hellwigii", tree.full$tip.label)], fix = c("Hydnophytum_hellwigii_Gay_487_(BM)", "Hydnophytum_hellwigii_Schlechter_13674_(BM)")))
indices.sub <- match(gsub(" ", "_", tiplabel.fix[, "phylo"]), tree.full$tip.label)

tree.full$tip.label <- setNames(sapply(tree.full$tip.label, function(x){paste(strsplit(x, split = "_")[[1]][1:2], collapse = "_")}), NULL)
tree.full$tip.label[indices.sub] <- gsub(" ", "_", tiplabel.fix[, "fix"])
tree.full$tip.label <- gsub("P_", "Psychotria_", tree.full$tip.label)
tree.full$tip.label <- gsub("_", " ", tree.full$tip.label)

indices.rownames <- match(gsub(" ", "_", tiplabel.fix[, "phylo"]), rownames(traits.plot))
names.temp <- rownames(traits.plot)
names.temp <- setNames(sapply(names.temp, function(x){paste(strsplit(x, split = "_")[[1]][1:2], collapse = "_")}), NULL)
names.temp[indices.rownames] <- gsub(" ", "_", tiplabel.fix[, "fix"])
names.temp <- gsub("P ", "Psychotria ", names.temp)
names.temp <- gsub("_", " ", names.temp)
rownames(traits.plot) <- names.temp

## Fan layout

p <- ggtree(tree.full, layout = "fan", open.angle = 30, xlim = c(0, 35)) %<+% traits

col1 <- setNames((grafify:::graf_palettes)$okabe_ito[c(1, 3, 4, 2, 7)], NULL)

p1 <-
    gheatmap(p = p, data = traits.plot[, c("dom.growth", "reward")], width = 0.1, colnames_angle = 45, hjust = 1) +
    geom_tiplab(size = 3, offset = 30) +
    theme_tree2() +
    scale_fill_manual(values = alpha(col1, 0.75), name = "DomGrowth")

p1.2 <- p1 + new_scale_fill()

## col2 <- setNames((grafify:::graf_palettes)$okabe_ito[c(1, 3, 6, 7, 8, 4)], NULL)
col2 <- setNames((grafify:::graf_palettes)$muted[c(2:8)], NULL)

p1.3 <-
    gheatmap(p = p1.2, data = traits.plot[, c("warts", "holediam.disc")], offset = 6.15, width = 0.1, colnames_angle = 45, hjust = 1) +
    #scale_fill_brewer(palette = "Set1", name = "Reward")
    scale_fill_manual(values = alpha(col2, 0.75))

p1.4 <- p1.3 + new_scale_fill()

p2 <-
    gheatmap(p = p1.4, data = traits.plot[, c("stem.area", "leaf.area", "corola.length", "petiole.length")], offset = 13, width = 0.25, colnames_angle = 45, hjust = 1) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral"), name = "Continuous Traits") +
    theme(legend.position = "left")

p2

ggsave(filename = here::here("output/preliminary_phylo_fan_4traits_withOutgroup.pdf"), plot = p2, width = 30, height = 30)




## Normal layout

p <- ggtree(tree.full) %<+% traits

#p <- ggtree(tree, layout = "fan", open.angle = 30) %<+% traits

## p1 <-
##     p +
##     #geom_tiplab(offset = 1, hjust = 0) +
##     theme(legend.position = "bottom") +
##     scale_shape_manual(values = c(15, 16, 17, 18)) +
##     scale_colour_brewer(palette = "Set1")

col1 <- setNames((grafify:::graf_palettes)$okabe_ito[c(1, 3, 4, 2, 7)], NULL)

p1 <-
    gheatmap(p = p, data = traits.plot[, c("dom.growth", "reward")], width = 0.1, colnames_angle = 45, hjust = 1) +
    geom_tiplab(size = 2, offset = 27.5) +
    theme_tree2() +
    scale_fill_manual(values = alpha(col1, 0.75), name = "DomGrowth") +
    xlim(0, 100) +
    ylim(-5, 140)

p1.2 <- p1 + new_scale_fill()

## col2 <- setNames((grafify:::graf_palettes)$okabe_ito[c(1, 3, 6, 7, 8, 4)], NULL)
col2 <- setNames((grafify:::graf_palettes)$muted[c(2:8)], NULL)

p1.3 <-
    gheatmap(p = p1.2, data = traits.plot[, c("warts", "holediam.disc")], offset = 5.9, width = 0.1, colnames_angle = 45, hjust = 1) +
    #scale_fill_brewer(palette = "Set1", name = "Reward")
    scale_fill_manual(values = alpha(col2, 0.75))

p1.4 <- p1.3 + new_scale_fill()

p2 <-
    gheatmap(p = p1.4, data = traits.plot[, c("stem.area", "leaf.area", "corola.length", "petiole.length")], offset = 13, width = 0.2, colnames_angle = 45, hjust = 1) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral"), name = "Continuous Traits") +
    theme(legend.position = "left")

p2


ggsave(filename = here::here("output/preliminary_phylo_regular_4traits_withOutgroup.pdf"), plot = p2, width = 20, height = 20)

































## Normal layout - tip labels close to phylogeny

p <- ggtree(tree) %<+% traits

## p1 <-
##     p +
##     #geom_tiplab(offset = 1, hjust = 0) +
##     theme(legend.position = "bottom") +
##     scale_shape_manual(values = c(15, 16, 17, 18)) +
##     scale_colour_brewer(palette = "Set1")

p1 <-
    gheatmap(p = p, data = traits.plot[, c("strategy", "architecture")], width = 0.1, colnames_angle = 45, hjust = 1, offset = 20) +
    #geom_tiplab(size = 2, offset = 40) +
    theme_tree2() +
    scale_fill_npg(alpha = 0.75, name = "Mutualistic Strategy \nand Plant Architecture")
p1.2 <- p1 + new_scale_fill()

p1.3 <-
    gheatmap(p = p1.2, data = traits.plot[, c("warts", "dom.growth")], offset = 26, width = 0.1, colnames_angle = 45, hjust = 1) +
    scale_fill_jco(alpha = 0.75, name = "Warts and\nDomatia Growth")

p1.4 <- p1.3 + new_scale_fill()

p1.5 <-
    gheatmap(p = p1.4, data = traits.plot[, c("reward", "mating.system")], offset = 32, width = 0.1, colnames_angle = 45, hjust = 1) +
    scale_fill_brewer(palette = "Set3", direction = -1, name = "Reward Type, Mating\nSystem and Leaf Structure")

p1.6 <- p1.5 + new_scale_fill()

p1.7 <-
    gheatmap(p = p1.6, data = traits.plot[, c("leaf.struct", "holediam.disc")], offset = 38, width = 0.1, colnames_angle = 45, hjust = 1) +
    scale_fill_brewer(palette = "Dark2", direction = 1, name = "Leaf Structure and\nHole Diameter (discrete)")

p1.8 <- p1.7 + new_scale_fill()

p2 <-
    gheatmap(p = p1.8, data = traits.plot[, c("stem.area", "leaf.area", "hole.diameter", "corola.length", "petiole.length")], offset = 46, width = 0.2, colnames_angle = 45, hjust = 1) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral"), name = "Continuous Traits") +
    theme(legend.position = "left")

#p2 + geom_tiplab(size = 2, offset = 0)

## p2.1 <- p2 + new_scale_fill()

## p3 <-
    ## gheatmap(p = p2.1, data = traits.plot[, c("hole.diameter", "corola.length", "petiole.length")], offset = 30, width = 0.1, colnames_angle = 45, hjust = 1) +
    ## theme(legend.position = "bottom") +
    ## scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral"), name = "Hole Diameter, Corola\nand Petiole Length")


ggsave(filename = here::here("output/preliminary_phylo_regular_labelsfirst.pdf"), plot = p2 + geom_tiplab(size = 2, offset = 0), width = 15, height = 20)





## Fan layout - tip labels close to phylogeny

p <- ggtree(tree, layout = "fan", open.angle = 30) %<+% traits

## p1 <-
##     p +
##     #geom_tiplab(offset = 1, hjust = 0) +
##     theme(legend.position = "bottom") +
##     scale_shape_manual(values = c(15, 16, 17, 18)) +
##     scale_colour_brewer(palette = "Set1")

p1 <-
    gheatmap(p = p, data = traits.plot[, c("strategy", "architecture")], width = 0.1, colnames_angle = 45, hjust = 1, offset = 20) +
    #geom_tiplab(size = 2, offset = 40) +
    theme_tree2() +
    scale_fill_npg(alpha = 0.75, name = "Mutualistic Strategy \nand Plant Architecture")
p1.2 <- p1 + new_scale_fill()

p1.3 <-
    gheatmap(p = p1.2, data = traits.plot[, c("warts", "dom.growth")], offset = 26, width = 0.1, colnames_angle = 45, hjust = 1) +
    scale_fill_jco(alpha = 0.75, name = "Warts and\nDomatia Growth")

p1.4 <- p1.3 + new_scale_fill()

p1.5 <-
    gheatmap(p = p1.4, data = traits.plot[, c("reward", "mating.system")], offset = 32, width = 0.1, colnames_angle = 45, hjust = 1) +
    scale_fill_brewer(palette = "Set3", direction = -1, name = "Reward Type, Mating\nSystem and Leaf Structure")

p1.6 <- p1.5 + new_scale_fill()

p1.7 <-
    gheatmap(p = p1.6, data = traits.plot[, c("leaf.struct", "holediam.disc")], offset = 38, width = 0.1, colnames_angle = 45, hjust = 1) +
    scale_fill_brewer(palette = "Dark2", direction = 1, name = "Leaf Structure and\nHole Diameter (discrete)")

p1.8 <- p1.7 + new_scale_fill()

p2 <-
    gheatmap(p = p1.8, data = traits.plot[, c("stem.area", "leaf.area", "hole.diameter", "corola.length", "petiole.length")], offset = 46, width = 0.2, colnames_angle = 45, hjust = 1) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral"), name = "Continuous Traits") +
    theme(legend.position = "left")

#p2 + geom_tiplab(size = 2, offset = 0)

## p2.1 <- p2 + new_scale_fill()

## p3 <-
    ## gheatmap(p = p2.1, data = traits.plot[, c("hole.diameter", "corola.length", "petiole.length")], offset = 30, width = 0.1, colnames_angle = 45, hjust = 1) +
    ## theme(legend.position = "bottom") +
    ## scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral"), name = "Hole Diameter, Corola\nand Petiole Length")


ggsave(filename = here::here("output/preliminary_phylo_fan_labelsfirst.pdf"), plot = p2 + geom_tiplab(size = 2, offset = 0), width = 15, height = 20)
