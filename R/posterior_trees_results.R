library("here")
library("ape")
library("phytools")
library("foreach")
library("doMC")
library("OUwie")
library("caper")
library("corHMM")
library("data.table")
library("expm")
library("ggplot2")
library("MetBrewer")
library("ggh4x")
library("patchwork")

i_am("R/posterior_trees_results.R")

for(i in 1:length(list.files(here("output/houwie/posterior_trees/"), pattern = ".RDS"))){
    assign(gsub("-", ".", gsub(".RDS", "", list.files(here("output/houwie/posterior_trees/"), pattern = ".RDS")[i])),
           readRDS(paste0("output/houwie/posterior_trees/", list.files(here("output/houwie/posterior_trees/"), pattern = ".RDS")[i]))
           )
}


## Domatium Growth
### Corola Length
model.pars.domgrow.corleng <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("domgrow.corleng.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("domgrow.corleng.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("domgrow.corleng.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.domgrow.corleng[[i]] <- avg.par
    model.pars.domgrow.corleng[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
domgrow.corleng.full <- plyr::ldply(model.pars.domgrow.corleng)
domgrow.corleng.full$tip_state <- factor(c("Apical", "Diffuse")[as.numeric(domgrow.corleng.full$tip_state)], levels = c("Apical", "Diffuse"))
domgrow.corleng.plot <- reshape2::melt(domgrow.corleng.full)
domgrow.corleng.plot$value <- round(domgrow.corleng.plot$value, 3)

dgclth <- ggplot(subset(domgrow.corleng.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Greek", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = element_blank(), y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
    theme(legend.position = "none",
          strip.text = element_text(size = 8))

dgclalp <- ggplot(subset(domgrow.corleng.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Greek", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = element_blank(), y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
    theme(legend.position = "none",
          strip.text = element_text(size = 8))

dgclsig <- ggplot(subset(domgrow.corleng.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Greek", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
    theme(legend.position = "none",
          strip.text = element_text(size = 8))

domgrow.corleng <-
    (dgclth / dgclalp / dgclsig) +
    plot_layout(guides = "collect", heights = c(0.33, 0.33, 0.33)) +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/domgrow_corleng_full.pdf"), plot = domgrow.corleng, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/domgrow_corleng_theta.pdf"), plot = dgclth, width = 10, height = 5)
## ggsave(here("output/houwie/posterior_trees/figs/domgrow_corleng_alpha.pdf"), plot = dgclalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/domgrow_corleng_sigma.pdf"), plot = dgclsig, width = 12, height = 6)


### Leaf Area
model.pars.domgrow.leafarea <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("domgrow.leafarea.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("domgrow.leafarea.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("domgrow.leafarea.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.domgrow.leafarea[[i]] <- avg.par
    model.pars.domgrow.leafarea[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
domgrow.leafarea.full <- plyr::ldply(model.pars.domgrow.leafarea)
domgrow.leafarea.full$tip_state <- factor(c("Apical", "Diffuse")[as.numeric(domgrow.leafarea.full$tip_state)], levels = c("Apical", "Diffuse"))
domgrow.leafarea.plot <- reshape2::melt(domgrow.leafarea.full)
domgrow.leafarea.plot$value <- round(domgrow.leafarea.plot$value, 3)

dglath <- ggplot(subset(domgrow.leafarea.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Greek", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

dglaalp <- ggplot(subset(domgrow.leafarea.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Greek", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

dglasig <- ggplot(subset(domgrow.leafarea.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Greek", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

domgrow.leafarea <-
    (dglath / dglaalp / dglasig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/domgrow_leafarea_full.pdf"), plot = domgrow.leafarea, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/domgrow_leafarea_theta.pdf"), plot = dglath, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/domgrow_leafarea_alpha.pdf"), plot = dglaalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/domgrow_leafarea_sigma.pdf"), plot = dglasig, width = 12, height = 6)



### Petiole Length
model.pars.domgrow.petleng <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("domgrow.petleng.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("domgrow.petleng.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("domgrow.petleng.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.domgrow.petleng[[i]] <- avg.par
    model.pars.domgrow.petleng[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
domgrow.petleng.full <- plyr::ldply(model.pars.domgrow.petleng)
domgrow.petleng.full$tip_state <- factor(c("Apical", "Diffuse")[as.numeric(domgrow.petleng.full$tip_state)], levels = c("Apical", "Diffuse"))
domgrow.petleng.plot <- reshape2::melt(domgrow.petleng.full)
domgrow.petleng.plot$value <- round(domgrow.petleng.plot$value, 3)

dgplth <- ggplot(subset(domgrow.petleng.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Greek", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

dgplalp <- ggplot(subset(domgrow.petleng.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Greek", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

dgplsig <- ggplot(subset(domgrow.petleng.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Greek", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

domgrow.petleng <-
    (dgplth / dgplalp / dgplsig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/domgrow_petleng_full.pdf"), plot = domgrow.petleng, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/domgrow_petleng_theta.pdf"), plot = dgplth, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/domgrow_petleng_alpha.pdf"), plot = dgplalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/domgrow_petleng_sigma.pdf"), plot = dgplsig, width = 12, height = 6)



### Stem Area
model.pars.domgrow.stemarea <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("domgrow.stemarea.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("domgrow.stemarea.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("domgrow.stemarea.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.domgrow.stemarea[[i]] <- avg.par
    model.pars.domgrow.stemarea[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
domgrow.stemarea.full <- plyr::ldply(model.pars.domgrow.stemarea)
domgrow.stemarea.full$tip_state <- factor(c("Apical", "Diffuse")[as.numeric(domgrow.stemarea.full$tip_state)], levels = c("Apical", "Diffuse"))
domgrow.stemarea.plot <- reshape2::melt(domgrow.stemarea.full)
domgrow.stemarea.plot$value <- round(domgrow.stemarea.plot$value, 3)

dgsath <- ggplot(subset(domgrow.stemarea.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Greek", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

dgsaalp <- ggplot(subset(domgrow.stemarea.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Greek", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

dgsasig <- ggplot(subset(domgrow.stemarea.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Greek", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

domgrow.stemarea <-
    (dgsath / dgsaalp / dgsasig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/domgrow_stemarea_full.pdf"), plot = domgrow.stemarea, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/domgrow_stemarea_theta.pdf"), plot = dgsath, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/domgrow_stemarea_alpha.pdf"), plot = dgsaalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/domgrow_stemarea_sigma.pdf"), plot = dgsasig, width = 12, height = 6)











## Holediam
### Corola Length
model.pars.holediam.corleng <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("holediam.corleng.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("holediam.corleng.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("holediam.corleng.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.holediam.corleng[[i]] <- avg.par
    model.pars.holediam.corleng[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
holediam.corleng.full <- plyr::ldply(model.pars.holediam.corleng)
holediam.corleng.full$tip_state <- factor(c("Several Large\nat Base", "One Large\nat Base", "All Large")[as.numeric(holediam.corleng.full$tip_state)], levels = c("Several Large\nat Base", "One Large\nat Base", "All Large"))
holediam.corleng.plot <- reshape2::melt(holediam.corleng.full)
holediam.corleng.plot$value <- round(holediam.corleng.plot$value, 3)

hdclth <- ggplot(subset(holediam.corleng.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Veronese", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

hdclalp <- ggplot(subset(holediam.corleng.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Veronese", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

hdclsig <- ggplot(subset(holediam.corleng.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Veronese", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

holediam.corleng <-
    (hdclth / hdclalp / hdclsig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/holediam_corleng_full.pdf"), plot = holediam.corleng, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/holediam_corleng_theta.pdf"), plot = hdclth, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/holediam_corleng_alpha.pdf"), plot = hdclalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/holediam_corleng_sigma.pdf"), plot = hdclsig, width = 12, height = 6)


### Leaf Area
model.pars.holediam.leafarea <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("holediam.leafarea.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("holediam.leafarea.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("holediam.leafarea.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.holediam.leafarea[[i]] <- avg.par
    model.pars.holediam.leafarea[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
holediam.leafarea.full <- plyr::ldply(model.pars.holediam.leafarea)
holediam.leafarea.full$tip_state <- factor(c("Several Large\nat Base", "One Large\nat Base", "All Large")[as.numeric(holediam.leafarea.full$tip_state)], levels = c("Several Large\nat Base", "One Large\nat Base", "All Large"))
holediam.leafarea.plot <- reshape2::melt(holediam.leafarea.full)
holediam.leafarea.plot$value <- round(holediam.leafarea.plot$value, 3)

hdlath <- ggplot(subset(holediam.leafarea.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Veronese", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

hdlaalp <- ggplot(subset(holediam.leafarea.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Veronese", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

hdlasig <- ggplot(subset(holediam.leafarea.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Veronese", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

holediam.leafarea <-
    (hdlath / hdlaalp / hdlasig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/holediam_leafarea_full.pdf"), plot = holediam.leafarea, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/holediam_leafarea_theta.pdf"), plot = hdlath, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/holediam_leafarea_alpha.pdf"), plot = hdlaalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/holediam_leafarea_sigma.pdf"), plot = hdlasig, width = 12, height = 6)



### Petiole Length
model.pars.holediam.petleng <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("holediam.petleng.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("holediam.petleng.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("holediam.petleng.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.holediam.petleng[[i]] <- avg.par
    model.pars.holediam.petleng[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
holediam.petleng.full <- plyr::ldply(model.pars.holediam.petleng)
holediam.petleng.full$tip_state <- factor(c("Several Large\nat Base", "One Large\nat Base", "All Large")[as.numeric(holediam.petleng.full$tip_state)], levels = c("Several Large\nat Base", "One Large\nat Base", "All Large"))
holediam.petleng.plot <- reshape2::melt(holediam.petleng.full)
holediam.petleng.plot$value <- round(holediam.petleng.plot$value, 3)

hdplth <- ggplot(subset(holediam.petleng.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Veronese", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

hdplalp <- ggplot(subset(holediam.petleng.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Veronese", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

hdplsig <- ggplot(subset(holediam.petleng.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Veronese", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

holediam.petleng <-
    (hdplth / hdplalp / hdplsig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/holediam_petleng_full.pdf"), plot = holediam.petleng, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/holediam_petleng_theta.pdf"), plot = hdplth, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/holediam_petleng_alpha.pdf"), plot = hdplalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/holediam_petleng_sigma.pdf"), plot = hdplsig, width = 12, height = 6)



### Stem Area
model.pars.holediam.stemarea <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("holediam.stemarea.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("holediam.stemarea.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("holediam.stemarea.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.holediam.stemarea[[i]] <- avg.par
    model.pars.holediam.stemarea[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
holediam.stemarea.full <- plyr::ldply(model.pars.holediam.stemarea)
holediam.stemarea.full$tip_state <- factor(c("Several Large\nat Base", "One Large\nat Base", "All Large")[as.numeric(holediam.stemarea.full$tip_state)], levels = c("Several Large\nat Base", "One Large\nat Base", "All Large"))
holediam.stemarea.plot <- reshape2::melt(holediam.stemarea.full)
holediam.stemarea.plot$value <- round(holediam.stemarea.plot$value, 3)

hdsath <- ggplot(subset(holediam.stemarea.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Veronese", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

hdsaalp <- ggplot(subset(holediam.stemarea.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Veronese", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

hdsasig <- ggplot(subset(holediam.stemarea.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Veronese", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

holediam.stemarea <-
    (hdsath / hdsaalp / hdsasig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/holediam_stemarea_full.pdf"), plot = holediam.stemarea, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/holediam_stemarea_theta.pdf"), plot = hdsath, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/holediam_stemarea_alpha.pdf"), plot = hdsaalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/holediam_stemarea_sigma.pdf"), plot = hdsasig, width = 12, height = 6)













## Reward
### Corola Length
model.pars.reward.corleng <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("reward.corleng.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("reward.corleng.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("reward.corleng.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.reward.corleng[[i]] <- avg.par
    model.pars.reward.corleng[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
reward.corleng.full <- plyr::ldply(model.pars.reward.corleng)
reward.corleng.full$tip_state <- factor(c("Absent", "Present")[as.numeric(reward.corleng.full$tip_state) + 1], levels = c("Absent", "Present"))
reward.corleng.plot <- reshape2::melt(reward.corleng.full)
reward.corleng.plot$value <- round(reward.corleng.plot$value, 3)

rewclth <- ggplot(subset(reward.corleng.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Hiroshige", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

rewclalp <- ggplot(subset(reward.corleng.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Hiroshige", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

rewclsig <- ggplot(subset(reward.corleng.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Hiroshige", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

reward.corleng <-
    (rewclth / rewclalp / rewclsig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/reward_corleng_full.pdf"), plot = reward.corleng, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/reward_corleng_theta.pdf"), plot = rewclth, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/reward_corleng_alpha.pdf"), plot = rewclalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/reward_corleng_sigma.pdf"), plot = rewclsig, width = 12, height = 6)


### Leaf Area
model.pars.reward.leafarea <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("reward.leafarea.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("reward.leafarea.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("reward.leafarea.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.reward.leafarea[[i]] <- avg.par
    model.pars.reward.leafarea[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
reward.leafarea.full <- plyr::ldply(model.pars.reward.leafarea)
reward.leafarea.full$tip_state <- factor(c("Absent", "Present")[as.numeric(reward.leafarea.full$tip_state) + 1], levels = c("Absent", "Present"))
reward.leafarea.plot <- reshape2::melt(reward.leafarea.full)
reward.leafarea.plot$value <- round(reward.leafarea.plot$value, 3)

rewlath <- ggplot(subset(reward.leafarea.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Hiroshige", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

rewlaalp <- ggplot(subset(reward.leafarea.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Hiroshige", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

rewlasig <- ggplot(subset(reward.leafarea.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Hiroshige", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

reward.leafarea <-
    (rewlath / rewlaalp / rewlasig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/reward_leafarea_full.pdf"), plot = reward.leafarea, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/reward_leafarea_theta.pdf"), plot = rewlath, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/reward_leafarea_alpha.pdf"), plot = rewlaalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/reward_leafarea_sigma.pdf"), plot = rewlasig, width = 12, height = 6)



### Petiole Length
model.pars.reward.petleng <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("reward.petleng.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("reward.petleng.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("reward.petleng.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.reward.petleng[[i]] <- avg.par
    model.pars.reward.petleng[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
reward.petleng.full <- plyr::ldply(model.pars.reward.petleng)
reward.petleng.full$tip_state <- factor(c("Absent", "Present")[as.numeric(reward.petleng.full$tip_state) + 1], levels = c("Absent", "Present"))
reward.petleng.plot <- reshape2::melt(reward.petleng.full)
reward.petleng.plot$value <- round(reward.petleng.plot$value, 3)

rewplth <- ggplot(subset(reward.petleng.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Hiroshige", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

rewplalp <- ggplot(subset(reward.petleng.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Hiroshige", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

rewplsig <- ggplot(subset(reward.petleng.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Hiroshige", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

reward.petleng <-
    (rewplth / rewplalp / rewplsig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/reward_petleng_full.pdf"), plot = reward.petleng, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/reward_petleng_theta.pdf"), plot = rewplth, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/reward_petleng_alpha.pdf"), plot = rewplalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/reward_petleng_sigma.pdf"), plot = rewplsig, width = 12, height = 6)



### Stem Area
model.pars.reward.stemarea <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("reward.stemarea.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("reward.stemarea.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("reward.stemarea.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.reward.stemarea[[i]] <- avg.par
    model.pars.reward.stemarea[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
reward.stemarea.full <- plyr::ldply(model.pars.reward.stemarea)
reward.stemarea.full$tip_state <- factor(c("Absent", "Present")[as.numeric(reward.stemarea.full$tip_state) + 1], levels = c("Absent", "Present"))
reward.stemarea.plot <- reshape2::melt(reward.stemarea.full)
reward.stemarea.plot$value <- round(reward.stemarea.plot$value, 3)

rewsath <- ggplot(subset(reward.stemarea.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Hiroshige", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

rewsaalp <- ggplot(subset(reward.stemarea.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Hiroshige", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

rewsasig <- ggplot(subset(reward.stemarea.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Hiroshige", n = 2, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

reward.stemarea <-
    (rewsath / rewsaalp / rewsasig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/reward_stemarea_full.pdf"), plot = reward.stemarea, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/reward_stemarea_theta.pdf"), plot = rewsath, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/reward_stemarea_alpha.pdf"), plot = rewsaalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/reward_stemarea_sigma.pdf"), plot = rewsasig, width = 12, height = 6)













## Warts
### Corola Length
model.pars.warts.corleng <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("warts.corleng.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("warts.corleng.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("warts.corleng.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.warts.corleng[[i]] <- avg.par
    model.pars.warts.corleng[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
warts.corleng.full <- plyr::ldply(model.pars.warts.corleng)
warts.corleng.full$tip_state <- factor(c("Variable", "Differentiated", "Lost")[as.numeric(warts.corleng.full$tip_state)], levels = c("Variable", "Differentiated", "Lost"))
warts.corleng.plot <- reshape2::melt(warts.corleng.full)
warts.corleng.plot$value <- round(warts.corleng.plot$value, 3)

warclth <- ggplot(subset(warts.corleng.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Isfahan2", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

warclalp <- ggplot(subset(warts.corleng.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Isfahan2", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

warclsig <- ggplot(subset(warts.corleng.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Isfahan2", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

warts.corleng <-
    (warclth / warclalp / warclsig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/warts_corleng_full.pdf"), plot = warts.corleng, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/warts_corleng_theta.pdf"), plot = warclth, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/warts_corleng_alpha.pdf"), plot = warclalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/warts_corleng_sigma.pdf"), plot = warclsig, width = 12, height = 6)


### Leaf Area
model.pars.warts.leafarea <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("warts.leafarea.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("warts.leafarea.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("warts.leafarea.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.warts.leafarea[[i]] <- avg.par
    model.pars.warts.leafarea[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
warts.leafarea.full <- plyr::ldply(model.pars.warts.leafarea)
warts.leafarea.full$tip_state <- factor(c("Variable", "Differentiated", "Lost")[as.numeric(warts.leafarea.full$tip_state)], levels = c("Variable", "Differentiated", "Lost"))
warts.leafarea.plot <- reshape2::melt(warts.leafarea.full)
warts.leafarea.plot$value <- round(warts.leafarea.plot$value, 3)

warlath <- ggplot(subset(warts.leafarea.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Isfahan2", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

warlaalp <- ggplot(subset(warts.leafarea.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Isfahan2", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

warlasig <- ggplot(subset(warts.leafarea.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Isfahan2", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

warts.leafarea <-
    (warlath / warlaalp / warlasig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/warts_leafarea_full.pdf"), plot = warts.leafarea, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/warts_leafarea_theta.pdf"), plot = warlath, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/warts_leafarea_alpha.pdf"), plot = warlaalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/warts_leafarea_sigma.pdf"), plot = warlasig, width = 12, height = 6)



### Petiole Length
model.pars.warts.petleng <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("warts.petleng.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("warts.petleng.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("warts.petleng.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.warts.petleng[[i]] <- avg.par
    model.pars.warts.petleng[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
warts.petleng.full <- plyr::ldply(model.pars.warts.petleng)
warts.petleng.full$tip_state <- factor(c("Variable", "Differentiated", "Lost")[as.numeric(warts.petleng.full$tip_state)], levels = c("Variable", "Differentiated", "Lost"))
warts.petleng.plot <- reshape2::melt(warts.petleng.full)
warts.petleng.plot$value <- round(warts.petleng.plot$value, 3)

warplth <- ggplot(subset(warts.petleng.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Isfahan2", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

warplalp <- ggplot(subset(warts.petleng.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Isfahan2", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

warplsig <- ggplot(subset(warts.petleng.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Isfahan2", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

warts.petleng <-
    (warplth / warplalp / warplsig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/warts_petleng_full.pdf"), plot = warts.petleng, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/warts_petleng_theta.pdf"), plot = warplth, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/warts_petleng_alpha.pdf"), plot = warplalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/warts_petleng_sigma.pdf"), plot = warplsig, width = 12, height = 6)



### Stem Area
model.pars.warts.stemarea <- vector(mode = "list", length = 20)
for(i in 1:20){
    print(i)
    mod.list <- c(paste0("warts.stemarea.er.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("warts.stemarea.sym.nohid.tree", i, "[[", 1:6, "]]"),
                  paste0("warts.stemarea.ard.nohid.tree", i, "[[", 1:6, "]]"))
    mod.tab <- getModelTable(lapply(mod.list, function(x){eval(parse(text = x))}))
    avg.par <- getModelAvgParams(lapply(mod.list, function(x){eval(parse(text = x))}))
    model.pars.warts.stemarea[[i]] <- avg.par
    model.pars.warts.stemarea[[i]]$tree <- paste0("tree", sprintf("%02d", i))
}
warts.stemarea.full <- plyr::ldply(model.pars.warts.stemarea)
warts.stemarea.full$tip_state <- factor(c("Variable", "Differentiated", "Lost")[as.numeric(warts.stemarea.full$tip_state)], levels = c("Variable", "Differentiated", "Lost"))
warts.stemarea.plot <- reshape2::melt(warts.stemarea.full)
warts.stemarea.plot$value <- round(warts.stemarea.plot$value, 3)

warsath <- ggplot(subset(warts.stemarea.plot, variable %in% "expected_mean"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Isfahan2", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~theta~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

warsaalp <- ggplot(subset(warts.stemarea.plot, variable %in% "alpha"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Isfahan2", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~alpha~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

warsasig <- ggplot(subset(warts.stemarea.plot, variable %in% "expected_var"), aes(x = tip_state, y = value, colour = tip_state, group = tree)) +
    geom_point(size = 1.5, shape = 21, alpha = 0.3) +
    #geom_line(aes(group = tree), colour = "darkgrey") +
    scale_colour_manual(values = met.brewer("Isfahan2", n = 3, type = "continuous")) +
    stat_summary(fun = mean, geom = "point", aes(group = 1), size = 1.5) +
    #stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group = 1), width = 0.15, colour = "black") +
    #ylim(-5, 5) +
    labs(x = "Tip State", y = expression('Parameter Value ('~sigma^2~')')) +
    facet_wrap( ~ tree, scales = "free_y", ncol = 4, strip.position = "right") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none",
          strip.text = element_text(size = 8))

warts.stemarea <-
    (warsath / warsaalp / warsasig) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(here("output/houwie/posterior_trees/figs/warts_stemarea_full.pdf"), plot = warts.stemarea, width = 9, height = 11)

## ggsave(here("output/houwie/posterior_trees/figs/warts_stemarea_theta.pdf"), plot = warsath, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/warts_stemarea_alpha.pdf"), plot = warsaalp, width = 12, height = 6)
## ggsave(here("output/houwie/posterior_trees/figs/warts_stemarea_sigma.pdf"), plot = warsasig, width = 12, height = 6)
