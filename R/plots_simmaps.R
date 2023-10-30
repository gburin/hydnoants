library("ape")
library("phytools")
library("foreach")
library("doMC")
library("OUwie")
library("caper")
library("MetBrewer")
library("ggplot2")

here::i_am("R/plot_simmaps.R")

simmaps.domgrow <- readRDS(here::here("output/simmaps_dom.growth.RDS"))
simmaps.reward <- readRDS(here::here("output/simmaps_reward.RDS"))
simmaps.warts <- readRDS(here::here("output/simmaps_warts.RDS"))
simmaps.holediam <- readRDS(here::here("output/simmaps_holediam.disc.RDS"))

summary.domgrow <- describe.simmap(simmaps.domgrow)
summary.reward <- describe.simmap(simmaps.reward)
summary.warts <- describe.simmap(simmaps.warts)
summary.holediam <- describe.simmap(simmaps.holediam)

plot(summary.domgrow, colors = setNames(colorRampPalette(met.brewer("Greek", n = 100, type = "continuous"))(3), c(1, 0, 2)))

plot(summary.reward, colors = setNames(colorRampPalette(met.brewer("Hiroshige", n = 100, type = "continuous"))(2), c(0, 1)))

plot(summary.warts, colors = setNames(colorRampPalette(met.brewer("Isfahan2", n = 100, type = "continuous"))(4), c(1, 0, 2, 3)))

plot(summary.holediam, colors = setNames(colorRampPalette(met.brewer("Veronese", n = 100, type = "continuous"))(4), c(1, 0, 2, 3)))
