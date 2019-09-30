
library(tidyverse)
library(GGally)
library(cowplot)


count_repeat_frame <- read_tsv("../count_repeat_frame.tsv")

subtelomere_levels <- c("1ptel", "1qtel", "2ptel", "2qtel", "3ptel",
                        "3qtel", "4ptel", "4qtel", "5ptel", "5qtel",
                        "6ptel", "6qtel", "7qtel", "7ptel", "8ptel", "8qtel",
                        "9ptel", "9qtel", "10ptel", "10qtel", "11ptel", "11qtel",
                        "12ptel", "12qtel", "13qtel", "14qtel", "15qtel", "16ptel",
                        "16qtel", "17ptel", "17qtel", "18ptel", "18qtel", "19ptel",
                        "19qtel", "20ptel", "20qtel", "21qtel", "22qtel",
                        "XpYptel", "Xqtel", "Yqtel")


sample_cell_type_sorted <- count_repeat_frame %>% select(sample, cell_type) %>% unique() %>% arrange(cell_type)

count_repeat_frame <- count_repeat_frame %>%
  mutate(counts_trn = counts / total_reads,
         chr = factor(chr,
                      levels = subtelomere_levels,
                      ordered = TRUE),
         sample = factor(sample,
                         levels = sample_cell_type_sorted$sample,
                         ordered = TRUE),
         has_repeat = tel29bp_repeats > 0,
         has_full_match = tel29bp_full_matches > 0,
         has_rajika = rajika_repeats > 0)


count_repeat_frame <- filter(count_repeat_frame,
                             cell_type != "GM847")



#lm_fit <- lm(rlog_trn ~ has_repeat + cell_type, data = filter(count_repeat_frame, cell_type != "GM847"))


p1 <- ggplot(count_repeat_frame %>% drop_na() %>% select(sample, cell_type),
             aes(x = sample, y = "Cell type", fill = cell_type)) +
  geom_tile() +
  scale_fill_viridis_d(option = "A") +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

p2 <- ggplot(count_repeat_frame %>% drop_na(),
             aes(y = chr, x = sample, fill = rlog_trn)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plot_grid(NULL,
          p1,
          p2,
          ncol = 1,
          nrow = 3,
          rel_heights = c(0.05, 0.05, 0.90),
          align = "v")


p3 <- ggplot(count_repeat_frame %>% select(chr, has_full_match) %>% unique(),
             aes(x = "Tel29bp full match", y = chr, fill = has_full_match)) +
  geom_tile() +
  scale_fill_viridis_d(option = "C") +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.y = element_blank())

p4 <- ggplot(count_repeat_frame %>% drop_na() %>% group_by(sample, cell_type) %>% summarise(total_counts = sum(counts)),
             aes(x = sample, y = total_counts), colour = "blue") +
  geom_col() +
  scale_y_continuous("Subtelomere\ncounts") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p5 <- ggplot(count_repeat_frame %>% drop_na() %>% select(sample, cell_type) %>%
               mutate(Alt = ifelse(cell_type %in% c("HEK293T", "HeLa"),
                                   "ALT(-)",
                                   "ALT(+)")),
                  aes(x = sample, y = "Alt", fill = Alt)) +
  geom_tile() +
  scale_fill_viridis_d(option = "C") +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())


plots <- align_plots(p4,
                     p1 + theme(legend.position = "none"),
                     p5 + theme(legend.position = "none"),
                     p2 + theme(legend.position = "none"),
                     align = "v",
                     axis = "l")

widths_for_grid <- c(0.9, 0.1)

grid_of_plots <- plot_grid(
  plot_grid(plots[[1]], NULL,
            ncol = 2,
            rel_widths = widths_for_grid),
  plot_grid(plots[[2]], NULL,
            ncol = 2,
            rel_widths = widths_for_grid),
  plot_grid(plots[[3]], NULL,
            ncol = 2,
            rel_widths = widths_for_grid),
  plot_grid(plots[[4]], p3 + theme(legend.position = "none"),
            ncol = 2,
            rel_widths = widths_for_grid,
            align = "h"),
  rel_heights = c(0.1, 0.04, 0.04, 0.92),
  ncol = 1,
  align = "v",
  axis = "lr"
)

legend_grid <- plot_grid(
  NULL,
  get_legend(p1),
  get_legend(p5),
  get_legend(p2),
  get_legend(p3),
  NULL,
  NULL,
  ncol = 1
)

main_visualisation <- plot_grid(
  grid_of_plots,
  legend_grid,
  ncol = 2,
  rel_widths = c(0.95, 0.5)
)



subtel_nuc_content <- read_tsv("../subtel_nuc_content.tsv") %>%
  mutate(subtelomere = factor(subtelomere,
                      levels = subtelomere_levels,
                      ordered = TRUE))

CpG_heatmap <- ggplot(subtel_nuc_content,
                      aes(x = pos,
                          y = subtelomere,
                          fill = CpG_mean_over_window)) +
  geom_raster() +
  scale_fill_viridis_c()

test_with_CpG <- plot_grid(p2, CpG_heatmap, align = "h")

pdf("paper_plot_exploration_1.pdf",
    paper = "a4")

print(main_visualisation)

print(test_with_CpG)

dev.off()



