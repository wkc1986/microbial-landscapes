
# setup -------------------------------------------------------------------


library(data.table)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
library(cowplot)

utils.dir <- "utils/"
output.dir <- "mapper-output/"
figs.dir <- "../figures/"
for (script in list.files(utils.dir, full.names = TRUE)) source(script)
source("load_cholera_data.R")
cholera.mapper <- read.mapper.graph(paste0(output.dir, "cholera"))
cholera.v2p <- merge.mpr.samples(cholera.mapper, gordon.samples)
# lay out Mapper graph without singletons
plotter <- function(v) {
  plot.mapper.graph(cholera.mapper$graph,
                    node = geom_node_point(aes_(size = ~size, fill = v),
                                           shape = 21),
                    seed = 1,
                    exclude.singletons = TRUE) +
    guides(size = FALSE) +
    coord_equal()
}
basin.palette <- "Set3"
na.color <- "grey50"

# fraction diarrhea -------------------------------------------------------

plotter(~f.state) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill = "fraction\ndiarrhea")
fstate <- last_plot()

# basins ------------------------------------------------------------------

plotter(~as.factor(basin)) +
  scale_fill_brewer(palette = basin.palette, na.value = na.color)
basin <- last_plot()

# basin time series -------------------------------------------------------

theme_set(theme_cowplot())
p2basin <- cholera.v2p[, .N, by = .(subject, diagnosis, id, hour, basin)]
p2basin[, time := mapply(function(diagnosis, id, hour) {
  if (diagnosis == "diarrhea") {
    hour
  } else {
    as.numeric(strsplit(id, "d")[[1]][[2]])
  }
}, diagnosis = diagnosis, id = id, hour = hour)]
p2basin[, time.unit := sapply(diagnosis, function(x) {
  if (x == "diarrhea") "hour" else "day"
})]
p2basin[, i := rank(basin), by = .(subject, id)]
p2basin[, basin := as.factor(basin)]
theme_set(theme_cowplot(font_size = 10))
setkey(p2basin, diagnosis)
p2basin["diarrhea"] %>%
  ggplot(aes(x = time, y = basin)) +
  geom_tile(aes(fill = basin)) +
  facet_grid(subject ~ ., scales = "free_y", switch = "y") +
  theme(axis.ticks.y = element_blank(), panel.spacing = unit(1, "points")) +
  labs(x = "hour", fill = "basin") +
  scale_fill_brewer(palette = basin.palette, na.value = na.color,
                    drop = FALSE) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), strip.placement = "outside")
basin.series <- last_plot()

# basin distribution ------------------------------------------------------


d2basin <- p2basin[, .(N = sum(N)), by = .(subject, diagnosis, basin)]
d2basin[, frac := N / sum(N), by = .(subject, diagnosis)]
setkey(d2basin, diagnosis)
d2basin["diarrhea"] %>%
  ggplot(aes(x = subject, y = frac)) +
  geom_col(aes(fill = basin)) +
  coord_flip() +
  scale_fill_brewer(palette = basin.palette, na.value = na.color, drop = FALSE) +
  facet_grid(subject ~ ., scales = "free_y") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), panel.spacing = unit(1, "points")) +
  labs(fill = "basin", y = "fraction samples") #+
basin.distribs <- last_plot()
basin.combined <- plot_grid(basin.series + theme(legend.position = "none"),
                            basin.distribs +
                              theme(strip.background = element_blank(),
                                    strip.text = element_blank(),
                                    legend.position = "none"),
                            align = "hv", axis = "t"
                            )

# basin correlation function ----------------------------------------------


cholera.persistence <- split(cholera.v2p, by = c("subject", "diagnosis")) %>%
  lapply(function(df) basin.persistence(df$basin, df$hour, scale = TRUE)) %>%
  rbindlist(idcol = "group") %>%
  .[, c("subject", "diagnosis") := tstrsplit(group, "\\.")] %>%
  .[, basin := as.factor(as.numeric(basin))]
setkey(cholera.persistence, diagnosis)
cholera.persistence["diarrhea"] %>%
  .[!is.na(basin)] %>%
  ggplot(aes(x = delta.t, y = f, color = basin)) +
  geom_smooth(aes(fill = basin)) +
  stat_summary(geom = "point", fun.y = mean, size = 0.1) +
  scale_color_brewer(palette = basin.palette, na.value = na.color, drop = FALSE) +
  scale_fill_brewer(palette = basin.palette, na.value = na.color,  drop = FALSE) +
  facet_wrap(~ subject, nrow = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "interval (hours)", y = "correlation")
correlation <- last_plot()


# combined figure ---------------------------------------------------------


plot_grid(plot_grid(fstate, basin,
                    labels = "AUTO", align = "h", nrow = 1),
          basin.combined,
          correlation + theme(legend.position = "none"),
          rel_heights = c(2, 1, 1),
          labels = c("", "C", "D"),
          label_y = c(1, 1.2, 1),
          ncol = 1)
save_plot(last_plot(), filename = paste0(figs.dir, "/fig2.pdf"),
          ncol = 1, base_width = 8, base_height = 8)

# event vertex distribution -----------------------------------------------

vd <- vertex.distribution(cholera.v2p$vertex, cholera.v2p[, .(subject, diagnosis)])
ggplot(vd, aes(x = vertex)) +
  geom_col(aes(y = f, fill = diagnosis)) +
  facet_wrap(~subject, scales = "free_y")

state <- paste(vd$diagnosis, vd$subject)
sd <- state.divergence(vd$f, vd$vertex, state)
ggplot(sd, aes(x = state.x, y = state.y)) +
  geom_tile(aes(fill = JSD)) +
  scale_fill_distiller(palette = "Greys")


# correlation between quasipotential and rate of change -------------------

consecutive.samples <- gordon.samples %>%
  split(by = "subject") %>%
  lapply(function(df) pair.consecutive(df$sample, df$hour)) %>%
  rbindlist(idcol = "subject")
consecutive.samples[, delta.hour := time.y - time.x]
consecutive.samples <- merge(consecutive.samples, jsd)
consecutive.samples[, roc := jsd / delta.hour]
ggplot(consecutive.samples[delta.hour <= 24], aes(x = delta.hour, y = jsd)) +
  geom_point(aes(color = subject)) +
  scale_y_log10()

# make distance matrix and knn
dm <- dlist2dm(jsd$sample.x, jsd$sample.y, jsd$jsd)
knn <- dist2knn(dm, round(nrow(dm) / 10))
knn <- data.frame(sample = names(knn), knn = knn)
# recalculate potential for each point from vertex potential
# potential <- cholera.v2p[, .(potential = mean(scaled.knn)), by = point]

consecutive.samples <- merge(consecutive.samples, knn, by.x = "sample.x",
                             by.y = "sample")
consecutive.samples <- merge(consecutive.samples, knn, by.x = "sample.y",
                             by.y = "sample")
ggplot(consecutive.samples[delta.hour <= 24], aes(x = knn.x, y = sqrt(jsd) / delta.hour)) +
  geom_point(aes(color = subject)) #+
  # scale_y_log10()


# frequency in jumps of decreasing or increasing knn ----------------------

consecutive.samples[, delta.knn := knn.y - knn.x]
ggplot(consecutive.samples, aes(x = delta.knn / delta.hour)) +
  stat_ecdf(aes(color = subject))
