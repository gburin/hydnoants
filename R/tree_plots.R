library("ape")
library("tidyverse")
library("cowplot")
library("ggtree")
library("ggforce")
library("ggnewscale")
library("patchwork")
library("ggsci")

tree <- read.tree("../data/Tree_Hydno_multi_traits.tre")

traits <- read.csv("../data/Dataset_traits_Hydnophytinae_CODED.csv", as.is = TRUE)
names(traits) <- c("species", "strategy", "warts", "stem.area", "leaf.area", "hole.diameter", "architecture", "dom.growth", "reward", "corola.length", "mating.system", "petiole.length", "leaf.struct", "appendages")

traits$species[grep("Wor2g", traits$species)] <- gsub("Wor2g", "Worthing", traits$species[grep("Wor2g", traits$species)])

traits$hole.diameter[which(traits$hole.diameter == "NA (outgroup has no domatia)")] <- NA
traits$hole.diameter <- as.numeric(traits$hole.diameter)

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

## 1, 2, 6, 7

traits.plot <- traits[,-1]
rownames(traits.plot) <- traits[,1]

p <- ggtree(tree, layout = "fan", open.angle = 30) %<+% traits

## p1 <-
##     p +
##     #geom_tiplab(offset = 1, hjust = 0) +
##     theme(legend.position = "bottom") +
##     scale_shape_manual(values = c(15, 16, 17, 18)) +
##     scale_colour_brewer(palette = "Set1")

p1 <-
    gheatmap(p = p, data = traits.plot[, c("strategy", "architecture")], width = 0.2, colnames_angle = 45, hjust = 1) +
    geom_tiplab(size = 2, offset = 30) +
    scale_fill_npg(alpha = 0.75, name = "Mutualistic Strategy \nand Plant Architecture") +
    theme(legend.position = "bottom")

p1.2 <- p1 + new_scale_fill()

p1.3 <-
    gheatmap(p = p1.2, data = traits.plot[, c("warts", "dom.growth")], offset = 4, width = 0.2, colnames_angle = 45, hjust = 1) +
    scale_fill_jco(alpha = 0.75, name = "Warts and\n Domatia Growth") +
    theme(legend.position = "bottom")

p1.4 <- p1.3 + new_scale_fill()

p1.5 <-
    gheatmap(p = p1.4, data = traits.plot[, c("reward", "mating.system", "leaf.struct")], offset = 8, width = 0.3, colnames_angle = 45, hjust = 1) +
    scale_fill_brewer(palette = "Set3", direction = -1, name = "Reward Type, Mating\nSystem and Leaf Structure") +
    theme(legend.position = "bottom")

p1.6 <- p1.5 + new_scale_fill()

p2 <-
    gheatmap(p = p1.6, data = traits.plot[, c(3, 4)], offset = 18, width = 0.2, colnames_angle = 45, hjust = 1) +
    theme(legend.position = "bottom") +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral"), name = "Stem and Leaf Area")

p2.1 <- p2 + new_scale_fill()

p3 <-
    gheatmap(p = p2.1, data = traits.plot[, c(5, 9, 11)], offset = 22, width = 0.3, colnames_angle = 45, hjust = 1) +
    theme(legend.position = "bottom") +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral"), name = "Hole Diameter, Corola\nand Petiole Length")


ggsave(filename = "../output/exploratory_phylo.pdf", width = 22, height = 22)
