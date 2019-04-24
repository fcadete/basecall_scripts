
library(tidyverse)

rhietman_mapont_files <- list.files(path = "alignment_outputs/", pattern = "*rhietman_mapont_primary.depth")

png("rhietman_coverage_primary.png", width = 1280, height = 1024)

depth_frame <- data_frame()

for (file in rhietman_mapont_files) {
  
  depth_frame <- rbind(depth_frame,
                       cbind(data_frame(file = file),
                             read_tsv(paste0("alignment_outputs/", file),
                                      col_names = c("chr", "pos", "depth"),
                                      col_types = "cii")))
  
}

depth_frame <- as_data_frame(depth_frame)

depth_frame <- depth_frame %>% separate(chr, into = c("chr"), extra = "drop")

depth_frame$chr <- factor(depth_frame$chr,
                          levels = c("1ptel",   "1qtel",  "2ptel",   "2qtel",   "3ptel",
                                     "3qtel",   "4ptel",  "4qtel",   "5ptel",   "5qtel",
                                     "6ptel",   "6qtel",  "7qtel", "7ptel", "8ptel",   "8qtel",
                                     "9ptel",  "9qtel",  "10ptel", "10qtel",  "11ptel",  "11qtel",
                                     "12ptel",  "12qtel", "13qtel", "14qtel", "15qtel",  "16ptel",
                                     "16qtel", "17ptel", "17qtel", "18ptel", "18qtel",  "19ptel", "19qtel",
                                      "20ptel", "20qtel", "21qtel",  "22qtel","XpYptel", "Xqtel", "Yqtel"),
                          ordered = TRUE)
                                       

for (this_file in unique(depth_frame$file)) {
  print(depth_frame %>% filter(file == this_file) %>%
          ggplot(aes(x = pos, y = depth)) +
          geom_line() +
          facet_wrap(~ chr) +
          labs(main = this_file) +
          scale_x_continuous(name = this_file))
}


left_join(data_frame(chr = factor(rep(unique(depth_frame$chr), each = 500000),
                                  levels = c("1ptel",   "1qtel",  "2ptel",   "2qtel",   "3ptel",
                                             "3qtel",   "4ptel",  "4qtel",   "5ptel",   "5qtel",
                                             "6ptel",   "6qtel",  "7qtel", "7ptel", "8ptel",   "8qtel",
                                             "9ptel",  "9qtel",  "10ptel", "10qtel",  "11ptel",  "11qtel",
                                             "12ptel",  "12qtel", "13qtel", "14qtel", "15qtel",  "16ptel",
                                             "16qtel", "17ptel", "17qtel", "18ptel", "18qtel",  "19ptel", "19qtel",
                                             "20ptel", "20qtel", "21qtel",  "22qtel","XpYptel", "Xqtel", "Yqtel"),
                                  ordered = TRUE), 
                     pos = rep(1:500000, length(unique(depth_frame$chr)))),
          depth_frame %>%
            filter(grepl("Joana_TERRA_VNP", file)) %>%
            group_by(chr, pos) %>%
            summarise(depth = sum(depth, na.rm = TRUE))) %>%
  mutate(depth = ifelse(is.na(depth), 0, depth)) %>%
  ggplot(aes(x = pos, y = depth)) + geom_line() + facet_wrap( ~ chr) +
  scale_x_continuous(name = "Distance from telomere repeats (bp)", breaks = c(1, 100000, 200000, 300000, 400000, 500000)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
