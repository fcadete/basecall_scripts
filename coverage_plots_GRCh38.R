
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)

GRCh38_mapont_files <- list.files(path = "alignment_outputs/",
                                    pattern = "*GRCh38_mapont_primary.depth")

GRCh38_mapont_files <- GRCh38_mapont_files[grepl("Joana_TERRA_VNP", GRCh38_mapont_files)]

png("GRCh38_coverage_primary.png", width = 1280, height = 1024)

depth_frame <- data_frame()

for (file in GRCh38_mapont_files) {
  
  depth_frame <- rbind(depth_frame,
                       cbind(data_frame(file = file),
                             read_tsv(paste0("alignment_outputs/", file),
                                      col_names = c("chr", "pos", "depth"),
                                      col_types = "cii")))
  
}

depth_frame <- as_data_frame(depth_frame)

depth_frame <- depth_frame %>% filter(chr %in% c(1:22, "X", "Y"))

depth_frame$chr <- factor(depth_frame$chr,
                          levels = c(1:22, "X", "Y"),
                          ordered = TRUE)


depth_frame <- depth_frame %>% mutate(cell_type = ifelse(grepl("HeL", file),
                                                         "HeLa", "U2OS"))
                                       

depth_frame %>%
  group_by(cell_type, chr, pos) %>%
  summarise(depth = sum(depth, na.rm = TRUE)) %>%
  filter(cell_type == "U2OS") %>%
  droplevels() %>%
  ggplot() +
  geom_point(mapping = aes(x = pos, y = depth)) +
  facet_wrap( ~ chr, scales = "free_x") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


dev.off()


