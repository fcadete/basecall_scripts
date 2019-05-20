
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)

rhietman_mapont_files <- list.files(path = "alignment_outputs/",
                                    pattern = "*rhietman_mapont_primary.depth")

depth_frame <- data_frame()

for (file in rhietman_mapont_files) {
  
  if (grepl("Joana_TERRA_VNP", file)) {
  depth_frame <- rbind(depth_frame,
                       cbind(data_frame(file = file),
                             read_tsv(paste0("alignment_outputs/", file),
                                      col_names = c("chr", "pos", "depth"),
                                      col_types = "cii")))
  }
}

depth_frame <- as_data_frame(depth_frame)

depth_frame <- depth_frame %>%
                 separate(chr, into = c("chr"), extra = "drop") %>%
                 separate(file, sep = "_Joana_TERRA_VNP", into = c("file"), extra = "drop")

depth_frame$chr <- factor(depth_frame$chr,
                          levels = c("1ptel", "1qtel", "2ptel", "2qtel", "3ptel",
                                     "3qtel", "4ptel", "4qtel", "5ptel", "5qtel",
                                     "6ptel", "6qtel", "7qtel", "7ptel", "8ptel", "8qtel",
                                     "9ptel", "9qtel", "10ptel", "10qtel", "11ptel", "11qtel",
                                     "12ptel", "12qtel", "13qtel", "14qtel", "15qtel", "16ptel",
                                     "16qtel", "17ptel", "17qtel", "18ptel", "18qtel", "19ptel",
                                     "19qtel", "20ptel", "20qtel", "21qtel", "22qtel",
                                     "XpYptel", "Xqtel", "Yqtel"),
                          ordered = TRUE)
                                       
tel_29_bp_frame <- read_delim("tel_29_bp_blast_results.bed",
                              delim = " ",
                              col_names=c("chr", "start", "end", "strand"))


tel_29_bp_frame <- tel_29_bp_frame %>% separate(chr, into = c("chr"), extra = "drop")

tel_29_bp_frame$chr <- factor(tel_29_bp_frame$chr,
                          levels = c("1ptel", "1qtel", "2ptel", "2qtel", "3ptel",
                                     "3qtel", "4ptel", "4qtel", "5ptel", "5qtel",
                                     "6ptel", "6qtel", "7qtel", "7ptel", "8ptel", "8qtel",
                                     "9ptel", "9qtel", "10ptel", "10qtel", "11ptel", "11qtel",
                                     "12ptel", "12qtel", "13qtel", "14qtel", "15qtel", "16ptel",
                                     "16qtel", "17ptel", "17qtel", "18ptel", "18qtel", "19ptel",
                                     "19qtel", "20ptel", "20qtel", "21qtel", "22qtel",
                                     "XpYptel", "Xqtel", "Yqtel"),
                          ordered = TRUE)


for (this_file in unique(depth_frame$file)) {
  ggsave(paste0("coverage_plots/coverage_", this_file, ".png"),
         depth_frame %>% filter(file == this_file) %>%
          group_by(chr, pos) %>%
          summarise(depth = sum(depth, na.rm = TRUE)) %>%
          ggplot(aes(x = pos, y = depth)) +
          geom_line() +
          facet_wrap(~ chr) +
          labs(main = this_file) +
          scale_x_continuous(name = this_file) +
          theme_bw(),
      width = 20, height = 15, units = "cm", device = png())

ggsave(paste0("coverage_plots/coverage_", this_file, "_zoomed.png"),
       left_join(data_frame(chr = factor(rep(unique(depth_frame$chr), each = 10000),
                                          levels = c("1ptel", "1qtel", "2ptel", "2qtel", "3ptel",
                                                     "3qtel", "4ptel", "4qtel", "5ptel", "5qtel",
                                                     "6ptel", "6qtel", "7qtel", "7ptel", "8ptel", "8qtel",
                                                     "9ptel", "9qtel", "10ptel", "10qtel", "11ptel",
                                                     "11qtel", "12ptel", "12qtel", "13qtel", "14qtel",
                                                     "15qtel", "16ptel", "16qtel", "17ptel", "17qtel",
                                                     "18ptel", "18qtel",  "19ptel", "19qtel",
                                                     "20ptel", "20qtel", "21qtel",  "22qtel",
                                                     "XpYptel", "Xqtel", "Yqtel"),
                                           ordered = TRUE),
                               pos = rep(1:10000, length(unique(depth_frame$chr)) )),
                  depth_frame %>% filter(file == this_file) %>% filter(pos <= 10000) %>%
                    group_by(chr, pos) %>%
                    summarise(depth = sum(depth, na.rm = TRUE)) %>% select(chr, pos, depth)) %>%
          mutate(depth = ifelse(is.na(depth), 0, depth)) %>%
          ggplot() +
          geom_line(mapping = aes(x = pos, y = depth)) +
          geom_rect(data = tel_29_bp_frame,
                    mapping = aes(xmin = start,
                                  xmax = end,
                                  ymin = 0,
                                  ymax = Inf),
                    alpha = I(1/2),
                    fill = "blue") +
       facet_wrap(~ chr) +
          labs(main = this_file) +
          scale_x_continuous(name = this_file,
                     breaks = c(1, 2000, 4000, 6000, 8000, 10000),
                     limits = c(1, 10000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)),
      width = 20, height = 15, units = "cm", device = png())


}



depth_frame <- left_join(data_frame(chr = rep(factor(rep(unique(depth_frame$chr), each = 500000),
                                  levels = c("1ptel", "1qtel", "2ptel", "2qtel", "3ptel",
                                             "3qtel", "4ptel", "4qtel", "5ptel", "5qtel",
                                             "6ptel", "6qtel", "7qtel", "7ptel", "8ptel", "8qtel",
                                             "9ptel", "9qtel", "10ptel", "10qtel", "11ptel",
                                             "11qtel", "12ptel", "12qtel", "13qtel", "14qtel",
                                             "15qtel", "16ptel", "16qtel", "17ptel", "17qtel",
                                             "18ptel", "18qtel",  "19ptel", "19qtel",
                                             "20ptel", "20qtel", "21qtel",  "22qtel",
                                             "XpYptel", "Xqtel", "Yqtel"),
                                  ordered = TRUE), 2), 
                     pos = rep(rep(1:500000, length(unique(depth_frame$chr)) ), 2 ),
                     cell_type = rep(c("HeLa", "U2OS"),
                                     each = (500000*length(unique(depth_frame$chr))))),
          depth_frame %>%
            mutate(cell_type = ifelse(grepl("HeLA", file), "HeLa", "U2OS")) %>%
            group_by(cell_type, chr, pos) %>%
            summarise(depth = sum(depth, na.rm = TRUE))) %>%
  mutate(depth = ifelse(is.na(depth), 0, depth))

png("coverage_plots/coverage_pooled_U2OS.png", width = 1280, height = 1024)

depth_frame %>%
  filter(cell_type == "U2OS") %>%
  ggplot() +
  geom_rect(data = tel_29_bp_frame,
            mapping = aes(xmin = start,
                          xmax = end,
                          ymin = 0,
                          ymax = Inf),
            alpha = I(1/2),
            fill = "blue") +
  geom_line(mapping = aes(x = pos, y = depth)) +
  facet_wrap( ~ chr) +
  scale_x_continuous(name = "Distance from telomere repeats (bp)",
                     breaks = c(1, 100000, 200000, 300000, 400000, 500000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

png("coverage_plots/coverage_pooled_HeLa.png", width = 1280, height = 1024)

depth_frame %>%
  filter(cell_type == "HeLa") %>%
  ggplot() +
  geom_rect(data = tel_29_bp_frame,
            mapping = aes(xmin = start,
                          xmax = end,
                          ymin = 0,
                          ymax = Inf),
            alpha = I(1/2),
            fill = "blue") +
  geom_line(mapping = aes(x = pos, y = depth)) +
  facet_wrap( ~ chr) +
  scale_x_continuous(name = "Distance from telomere repeats (bp)",
                     breaks = c(1, 100000, 200000, 300000, 400000, 500000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



dev.off()

png("coverage_plots/coverage_pooled_U2OS_zoomed.png", width = 1280, height = 1024)

depth_frame %>%
  filter(cell_type == "U2OS") %>%
  filter(pos <= 10000) %>%
  droplevels() %>%
  ggplot() +
  geom_rect(data = tel_29_bp_frame,
            mapping = aes(xmin = start,
                          xmax = end,
                          ymin = 0,
                          ymax = Inf),
            alpha = I(1/2),
            fill = "blue") +
  geom_line(mapping = aes(x = pos, y = depth)) +
  facet_wrap( ~ chr) +
  scale_x_continuous(name = "Distance from telomere repeats (bp)",
                     breaks = c(1, 2000, 4000, 6000, 8000, 10000),
                     limits = c(1, 10000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

png("coverage_plots/coverage_pooled_HeLa_zoomed.png", width = 1280, height = 1024)

depth_frame %>%
  filter(cell_type == "HeLa") %>%
  filter(pos <= 10000) %>%
  droplevels() %>%
  ggplot() +
  geom_rect(data = tel_29_bp_frame,
            mapping = aes(xmin = start,
                          xmax = end,
                          ymin = 0,
                          ymax = Inf),
            alpha = I(1/2),
            fill = "blue") +
  geom_line(mapping = aes(x = pos, y = depth)) +
  facet_wrap( ~ chr) +
  scale_x_continuous(name = "Distance from telomere repeats (bp)",
                     breaks = c(1, 2000, 4000, 6000, 8000, 10000),
                     limits = c(1, 10000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


