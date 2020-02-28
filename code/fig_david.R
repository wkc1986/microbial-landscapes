# setup -------------------------------------------------------------------


library(data.table)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
library(cowplot)

mapper.dir <- "mapper-output/david/"
utils.dir <- "utils/"
figs.dir <- "../figures/"
for (script in list.files(utils.dir, full.names = TRUE)) source(script)
source("load_david_data.R")
david.mapper <- read.mapper.graph(mapper.dir)
david.v2p <- merge.mpr.samples(david.mapper, david.samples)
plotter <- function(v) {
  plot.mapper.graph(david.mapper$graph,
                    node = geom_node_point(aes_(size = ~size, fill = v),
                                           shape = 21),
                    seed = 13,
                    exclude.singletons = TRUE) +
    guides(size = FALSE) +
    scale_size_area(max_size = 4)
}

# frac subject ------------------------------------------------------------

plotter(~f.subject) +
  scale_fill_distiller(palette = "Spectral", breaks = c(0, 0.5, 1),
                       labels = c("all B", "", "all A")) +
  labs(fill = "subject")
fsubject <- last_plot()


# events ------------------------------------------------------------------


set.seed(13)
lo <- david.mapper$graph %>% # mapper graph with no outlier
  activate(nodes) %>%
  filter(!in.singleton) %>%
  create_layout("fr", niter = 1000)
title.size <- 10
base.size <- 8
plt.event <- function(subj, ev, color, label) {
  setkey(david.v2p, subject, event)
  vs <- david.v2p[.(subj, ev)]$vertex
  ggraph(lo) +
    geom_edge_link0() +
    geom_node_point(aes(size = size), color = "grey") +
    geom_node_point(aes(size = size),
                    data = function(df) filter(df, vertex %in% vs),
                    fill = color, shape = 21) +
    labs(title = label) +
    theme_graph(base_family = "Helvetica", title_size = title.size,
                      base_size = base.size) +
    theme(legend.position = "none", plot.margin = margin(1, 1, 1, 1)) +
    scale_size_area(max_size = 2) +
    coord_equal()
}
plot_grid(plotlist = mapply(plt.event,
                            subj = c("A", "A", "B", "B"),
                            ev = c("US (pre)", "US (post)",
                                      "pre-Salmonella", "post-Salmonella"
                                      ),
                            color = c("blue", "purple", "red", "orange"),
                            label = c("A, pre-travel", "A, post-travel",
                                      "B, pre-Salmonella", "B, post-Salmonella"
                                      ),
                            SIMPLIFY = FALSE)
          )
events <- last_plot()


# basins ------------------------------------------------------------------

plotter(~as.factor(basin)) + labs(fill = "basin")
basins <- last_plot()


# basin series -----------------------------------------------


david.v2p[, basin := as.factor(basin)]
david.sample.basins <- david.v2p[, .N, by = .(subject, day, point, event, basin)]
theme_set(theme_cowplot())
david.sample.basins %>%
  ggplot(aes(x = day, y = basin)) +
  geom_tile(data = function(dt) dt[is.na(basin)], fill = "grey50") +
  geom_tile(aes(fill = event), data = function(dt) dt[!is.na(basin)]) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_discrete(drop = FALSE) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(subject ~ .)
series <- last_plot()

# distribs ----------------------------------------------------------------

# david.sample.basins %>%
#   group_by(subject, event) %>%
#   mutate(frac = N / sum(N)) %>%
#   group_by(subject, event, basin) %>%
#   summarize(frac = sum(frac)) %>%
#   filter(!is.na(basin)) %>%
  # --- COLOR BAR ---
  # ggplot(aes(x = event)) +
  # geom_col(aes(y = frac, group = basin, fill = basin),
  #          position = "stack") +
  # coord_flip() +
  # labs(y = "fraction samples") +
  # facet_grid(subject ~ ., scales = "free")
  # --- HEATMAP-ISH ---
  # ggplot(aes(x = basin, y = event)) +
  # geom_tile(aes(fill = frac)) +
  # # scale_fill_distiller(palette = "Greys", direction = -1) +
  # facet_grid(subject ~ ., scales = "free")
  # --- REGULAR BAR ---
  # ggplot(aes(x = basin, y = frac)) +
  # geom_col() +
  # facet_wrap(~ subject + event)

david.sample.basins %>%
  group_by(subject, event, basin) %>%
  summarise(N = sum(N)) %>%
  group_by(subject, event) %>%
  mutate(frac = N / sum(N)) %>%
  mutate(label = factor(paste(subject, event),
                        levels = c("A US (pre)", "A US (post)", "A travel",
                                   "A travel + diarrhea 1",
                                   "A travel + diarrhea 2",
                                   "B pre-Salmonella", "B post-Salmonella",
                                   "B Salmonella"))) %>%
  filter(label %in% c("A US (pre)", "A US (post)",
                      "B pre-Salmonella", "B post-Salmonella")) %>%
  filter(!is.na(basin)) %>%
  ggplot(aes(x = basin, y = frac)) +
  geom_col() +
  facet_wrap(~ label, dir = "v", nrow = 2) +
  labs(x = "state", y = "frequency") +
  scale_y_continuous(breaks = c(0, 0.4)) +
  theme(axis.text.x = element_blank())

distribs <- last_plot()

# correlation function ----------------------------------------------------

david.persistence <- david.v2p %>%
  split(by = c("subject", "event", "healthy")) %>%
  lapply(function(df) basin.persistence(df$basin, df$day, scale = TRUE)) %>%
  rbindlist(idcol = "id") %>%
  .[, basin := as.factor(as.numeric(basin))] %>%
  .[, N := uniqueN(delta.t), by = .(id, basin)] %>%
  .[, c("subject", "event", "healthy") := tstrsplit(id, "\\.",
                                                    type.convert = TRUE)] %>%
  .[, id := NULL] %>%
  .[!is.na(basin)]

# plot correlation for only the 3 most frequent basins per  phenotype
rank.states <- david.v2p[, .(count = .N), by = .(subject, event, basin)]
rank.states <- rank.states[!is.na(basin)]
rank.states[, rank := frank(-count), by = .(subject, event)]
top <- rank.states[rank <= 3]
setkey(david.persistence, subject, event, basin)
setkey(top, subject, event, basin)
david.persistence <- david.persistence[top]

theme_set(theme_cowplot() + theme(#legend.position = "none",
                                  axis.text = element_text(size = title.size),
                                  axis.title = element_text(size = title.size),
                                  strip.text = element_text(size = title.size)))
graphic.size <- 1
plot.correlation <- function(df) {
  ggplot(df, aes(x = delta.t, y = f)) +
  geom_point(size = graphic.size / 5) +
  geom_smooth(size = graphic.size) +
  scale_color_brewer(palette = "Paired", drop = TRUE) +
  scale_fill_brewer(palette = "Paired", drop = TRUE) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "interval (days)", y = "correlation", color = "state", fill = "state") +
  # facet_grid(paste(subject, event, sep = ", ") ~ .) +
  facet_wrap(~ state, dir = "v") +
  theme(strip.text.y = element_text(angle = 0))
}
david.persistence[!is.na(basin) & healthy & event != "travel"] %>%
  filter(!is.na(basin)) %>%
  filter(healthy) %>%
  filter(event != "travel") %>%
  mutate(label = paste(subject, event)) %>%
  mutate(label = factor(label, levels = c("A US (pre)", "A US (post)",
                                          "B pre-Salmonella", "B post-Salmonella"))) %>%
  mutate(state = factor(paste("state", basin),
                        levels = paste("state", levels(basin)))) %>%
  # split(.$label) %>%
  # lapply(plot.correlation)
  ggplot(aes(x = delta.t, y = f)) +
  geom_point(aes(color = state), size = graphic.size / 5) +
  geom_smooth(aes(color = state, fill = state), size = graphic.size) +
  scale_color_brewer(palette = "Paired", drop = TRUE) +
  scale_fill_brewer(palette = "Paired", drop = TRUE) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "interval (days)", y = "correlation", color = "state", fill = "state") +
  # facet_grid(paste(subject, event, sep = ", ") ~ .) +
  facet_wrap(~ label, ncol = 2, dir = "v") +
  theme(strip.text.y = element_text(angle = 0))
# plist <- mapply(function(g, l) g + ggtitle(l), g = plist, l = names(plist),
#                 SIMPLIFY = FALSE)
# plot_grid(plist[[1]], plist[[3]], plist[[2]], plist[[4]], nrow = 2)
correlation <- last_plot()


# compiled ----------------------------------------------------------------

# plot_grid(fsubject, events,
#           basins, plot_grid(distribs, correlation, nrow = 2),
#           ncol = 2, rel_heights = c(2, 3), labels = "AUTO")
plot_grid(fsubject, events, ncol = 2, labels = "AUTO")
save_plot(paste0(figs.dir, "/fig3.pdf"), last_plot(), ncol = 2,
          base_width = 4, base_height = 3)
plot_grid(distribs, correlation, nrow = 2, labels = "AUTO", align = "v",
          axis = "lr")
save_plot(paste0(figs.dir, "/fig4.pdf"), last_plot(), nrow = 1,
          base_width = 8, base_height = 6)
basins
save_plot(paste0(figs.dir, "/sup_fig1.pdf"), last_plot(), base_width = 8,
          base_height = 6)


# correlation between starting and change in knn --------------------------

dm <- dlist2dm(jsds$sample.x, jsds$sample.y, sqrt(jsds$jsd))
knn <- dist2knn(dm, round(nrow(dm) / 10))
knn <- data.frame(sample = names(knn), knn = knn)
david.samples <- merge(david.samples, knn, by = "sample")
emp.cor <- david.samples %>%
  group_by(subject, healthy, event) %>%
  summarise(pearson = knn.test(sample, day, knn))
nnulls <- 100
null.cor <- lapply(seq(nnulls), function(i) {
  david.samples %>%
    group_by(subject, healthy, event) %>%
    mutate(day = sample(day)) %>%
    summarise(pearson = knn.test(sample, day, knn))
})
null.cor <- rbindlist(null.cor, idcol = "iteration")
lvls <- c("US (pre)", "travel", "travel + diarrhea 1", "travel + diarrhea 2",
  "US (post)", "pre-Salmonella", "Salmonella", "post-Salmonella")
emp.cor <- mutate(emp.cor, event = factor(event, levels = lvls))
null.cor <- mutate(null.cor, event = factor(event, levels = lvls))
ggplot(null.cor, aes(x = pearson)) +
  geom_histogram() +
  geom_vline(aes(xintercept = pearson), data = emp.cor, color = "blue") +
  facet_wrap(~ subject + event)
save_plot(paste0(figs.dir, "/sup_fig7.pdf"), last_plot(), base_width = 8,
          base_height = 6)

# fraction subject for each state -----------------------------------------

david.v2p %>%
  group_by(basin) %>%
  summarize(f.subject = mean(f.subject)) %>%
  ggplot(aes(x = basin, y = f.subject)) +
  geom_col() +
  labs(x = "state", y = "fraction subject A")
ggsave("../figures/sup_fig_fsubject.pdf", width = 7, height = 5)


# mean taxon composition for each state -----------------------------------

taxonomy <- read.taxonomy("../data/david/david.otu.taxonomy")
sample.counts <- merge(david, taxonomy, by = "otu") %>%
  group_by(sample, kingdom, phylum, class, order, family, genus) %>%
  summarize(count = sum(count))
state.profiles <- merge(david.v2p, sample.counts,
                        by.x = "point", by.y = "sample", allow.cartesian = TRUE) %>%
  group_by(basin, kingdom, phylum, class, order, family, genus) %>%
  summarize(count = sum(count)) %>%
  group_by(basin) %>%
  mutate(relative.abundance = count / sum(count)) %>%
  mutate(count = NULL) %>%
  as.data.table %>%
  setnames("basin", "state")
if (!dir.exists("../data/supplemental-data/")) {
  dir.create("../data/supplemental-data/")
}
fwrite(state.profiles, "../data/supplemental-data/david-state-compositions.txt",
       sep = "\t", na = "NA")
