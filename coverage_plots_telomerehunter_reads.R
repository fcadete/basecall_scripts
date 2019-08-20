
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)


subtelomere_levels <- c("1ptel", "1qtel", "2ptel", "2qtel", "3ptel",
                        "3qtel", "4ptel", "4qtel", "5ptel", "5qtel",
                        "6ptel", "6qtel", "7qtel", "7ptel", "8ptel", "8qtel",
                        "9ptel", "9qtel", "10ptel", "10qtel", "11ptel", "11qtel",
                        "12ptel", "12qtel", "13qtel", "14qtel", "15qtel", "16ptel",
                        "16qtel", "17ptel", "17qtel", "18ptel", "18qtel", "19ptel",
                        "19qtel", "20ptel", "20qtel", "21qtel", "22qtel",
                        "XpYptel", "Xqtel", "Yqtel")

telomerehunter_depth_files <- list.files(path = "telomere_hunter_results/depths_from_telomeric_reads/",
                                    pattern = "*filtered_intratelomeric_primary.depth")

depth_frame <- data_frame()

for (file in telomerehunter_depth_files) {
  
  this_file_depth <- read_tsv(paste0("telomere_hunter_results/depths_from_telomeric_reads/", file),
                                      col_names = c("chr", "pos", "depth"),
                                      col_types = "cii")

  if (length(this_file_depth) > 0) {
  depth_frame <- rbind(depth_frame,
                       cbind(data_frame(file = file),
                             this_file_depth))
  }
}

depth_frame <- as_data_frame(depth_frame)

depth_frame <- depth_frame %>%
                 separate(chr, into = c("chr"), extra = "drop") %>%
                 separate(file, sep = "_filtered_intratelomeric_primary.depth",
                          into = c("file"), extra = "drop")

depth_frame$chr <- factor(depth_frame$chr,
                          levels = subtelomere_levels,
                          ordered = TRUE)
                                       
tel_29_bp_frame <- read_delim("tel_29bp_blast_results.txt",
                              delim = "\t",
                              col_names=c("qseqid", "chr", "pident",
                                          "length", "mismatch", "gapopen",
                                           "qstart", "qend",
                                           "start", "end",
                                           "evalue", "bitscore"))


tel_29_bp_frame <- tel_29_bp_frame %>% separate(chr, into = c("chr"), extra = "drop")

tel_29_bp_frame$chr <- factor(tel_29_bp_frame$chr,
                          levels = subtelomere_levels,
                          ordered = TRUE)

depth_frame <- depth_frame %>%
     mutate(cell_type = ifelse(grepl("HeL", file),
                               "HeLa",
                               ifelse(grepl("GM847", file),
                                      "GM847",
                                      ifelse(grepl("SAOS2", file),
                                             "SAOS2",
                                             ifelse(grepl("HEK", file),
                                                    "HEK293",
                                                    "U2OS")))))

cell_type_file_count <- depth_frame %>% count(cell_type, file)

depth_frame_by_file <- left_join(data_frame(chr = rep(factor(rep(unique(depth_frame$chr), each = 500000),
                                                             levels = subtelomere_levels,
                                                             ordered = TRUE),
                                                      nrow(cell_type_file_count)), 
                                            pos = rep(rep(1:500000, length(unique(depth_frame$chr)) ),
                                                      nrow(cell_type_file_count)),
                                            cell_type = rep(cell_type_file_count$cell_type,
                                                            each = (500000*length(unique(depth_frame$chr)))),
                                            file = rep(cell_type_file_count$file,
                                                       each = (500000*length(unique(depth_frame$chr))))),
          depth_frame %>%
            group_by(cell_type, file, chr, pos) %>%
            summarise(depth = sum(depth, na.rm = TRUE))) %>%
  mutate(depth = ifelse(is.na(depth), 0, depth))

system("wc -l *porechop_finalcall_table > total_reads_per_file.txt")
total_reads <- read_table("total_reads_per_file.txt",
                          col_names = c("total_reads", "file")) %>%
    filter(file != "total") %>% 
    mutate(file = str_remove(file, "guppy_"),
           file = str_remove(file, ".porechop_finalcall_table"))

system("grep Joana_TERRA *.porechop_finalcall_table > reads_with_terra_primer.txt")
primed_reads <- read_table("reads_with_terra_primer.txt",
                           col_names = c("file")) %>%
   separate(file, ":", into = c("file"), extra = "drop") %>%
   mutate(file = str_remove(file, "guppy_"),
          file = str_remove(file, ".porechop_finalcall_table")) %>%
   count(file, name = "primed_reads")

subtelomere_end_reads <- read_delim("soft_clipped_end_counts.txt", col_names = c("file", "subtelomere_end_reads"), delim = " ") %>%
        filter(grepl("Joana_TERRA", file)) %>%
        mutate(file = str_remove(file, "_Joana_TERRA_VNP_no_PCR"),
               file = str_remove(file, "_Joana_TERRA_VNP"),
               file = str_remove(file, "_on_rhietman_mapont_primary_subtelomere_start.bam")) %>%
        group_by(file) %>%
        summarise(subtelomere_end_reads = sum(subtelomere_end_reads))
        
full_join(full_join(total_reads, primed_reads), subtelomere_end_reads) %>%
     select(file, everything()) %>%
     write_tsv("table_for_joana.txt")

pdf("telomere_hunter_results/coverage_plots/coverage_normalized.pdf", width = 16, height = 12)

depth_frame_by_file %>%
   drop_na() %>%
   filter(pos <= 2000,
          file %in% c("VNP_TERRA_20190327",
                      "VNP-TERRA_purified_190403",
                      "VNP_purified_TERRA_U2OS_20190514")) %>%
   left_join(total_reads) %>%
   ggplot() +
      geom_rect(data = tel_29_bp_frame %>% filter(end <= 2000),
                mapping = aes(xmin = start,
                              xmax = end,
                              ymin = 0,
                              ymax = Inf,
                              fill = bitscore)) +
      scale_fill_gradient(low = "pink", high = "red") +
      geom_line(mapping = aes(x = pos,
                              y = depth / total_reads,
                              colour = file),
                alpha = I(3/5)) +
      facet_wrap(~ chr) +
      scale_colour_viridis_d() +
      theme_bw()


depth_frame_by_file %>%
   drop_na() %>%
   filter(pos <= 2000) %>%
   left_join(total_reads) %>%
   ggplot() +
      geom_rect(data = tel_29_bp_frame %>% filter(end <= 2000),
                mapping = aes(xmin = start,
                              xmax = end,
                              ymin = 0,
                              ymax = Inf,
                              fill = bitscore)) +
      scale_fill_gradient(low = "pink", high = "red") +
      geom_line(mapping = aes(x = pos,
                              y = depth / total_reads,
                              colour = file),
                alpha = I(3/5)) +
      facet_wrap(~ chr) +
      scale_colour_viridis_d(option = "plasma") +
      theme_bw()


depth_frame_by_file %>%
   drop_na() %>%
   filter(pos <= 2000) %>%
   left_join(total_reads) %>%
   ggplot() +
      geom_rect(data = tel_29_bp_frame %>% filter(bitscore == 58),
                mapping = aes(xmin = start,
                              xmax = end,
                              ymin = 0,
                              ymax = Inf),
            alpha = I(1/2), colour = "red") +
      geom_line(mapping = aes(x = pos,
                              y = depth / total_reads,
                              colour = file),
                alpha = I(3/5)) +
      facet_wrap(~ chr) +
      scale_colour_viridis_d(option = "plasma") +
      theme_bw()



depth_frame_by_file %>%
   drop_na() %>%
   filter(pos <= 2000) %>%
   left_join(total_reads) %>%
   ggplot() +
      geom_rect(data = tel_29_bp_frame %>% filter(end <= 2000),
                mapping = aes(xmin = start,
                              xmax = end,
                              ymin = 0,
                              ymax = Inf,
                              fill = bitscore)) +
      scale_fill_gradient(low = "pink", high = "red") +
      geom_line(mapping = aes(x = pos,
                              y = depth / total_reads,
                              group = file,
                              colour = cell_type),
                alpha = I(3/5)) +
      facet_wrap(~ chr) +
      scale_colour_viridis_d(option = "plasma") +
      theme_bw()


depth_frame_by_file %>%
   drop_na() %>%
   filter(pos <= 2000) %>%
   left_join(total_reads) %>%
   ggplot() +
      geom_rect(data = tel_29_bp_frame %>% filter(bitscore == 58),
                mapping = aes(xmin = start,
                              xmax = end,
                              ymin = 0,
                              ymax = Inf),
            alpha = I(1/2), colour = "red") +
      geom_line(mapping = aes(x = pos,
                              y = depth / total_reads,
                              group = file,
                              colour = cell_type),
                alpha = I(3/5)) +
      facet_wrap(~ chr) +
      scale_colour_viridis_d(option = "plasma") +
      theme_bw()


depth_frame <- filter(depth_frame, grepl("purified", file))


depth_frame <- left_join(data_frame(chr = rep(factor(rep(unique(depth_frame$chr), each = 500000),
                                                     levels = subtelomere_levels,
                                                     ordered = TRUE),
                                              length(unique(depth_frame$cell_type))), 
                                    pos = rep(rep(1:500000, length(unique(depth_frame$chr)) ),
                                              length(unique(depth_frame$cell_type))),
                                    cell_type = rep(unique(depth_frame$cell_type),
                                                    each = (500000*length(unique(depth_frame$chr))))),
          depth_frame %>%
            group_by(cell_type, chr, pos) %>%
            summarise(depth = sum(depth, na.rm = TRUE))) %>%
  mutate(depth = ifelse(is.na(depth), 0, depth))

total_reads_by_cell_type <- total_reads %>%
     filter(grepl("purified", file)) %>%
     mutate(cell_type = ifelse(grepl("HeL", file),
                               "HeLa",
                               ifelse(grepl("GM847", file),
                                      "GM847",
                                      ifelse(grepl("SAOS2", file),
                                             "SAOS2",
                                             ifelse(grepl("HEK", file),
                                                    "HEK293",
                                                    "U2OS"))))) %>%
     group_by(cell_type) %>%
     summarise(total_reads = sum(total_reads))


depth_frame %>%
   drop_na() %>%
   filter(pos <= 2000) %>%
   left_join(total_reads_by_cell_type) %>%
   ggplot() +
      geom_rect(data = tel_29_bp_frame %>% filter(end <= 2000),
                mapping = aes(xmin = start,
                              xmax = end,
                              ymin = 0,
                              ymax = Inf,
                              fill = bitscore)) +
      scale_fill_gradient(low = "pink", high = "red") +
      geom_line(mapping = aes(x = pos,
                              y = depth / total_reads,
                              colour = cell_type),
                alpha = I(3/5)) +
      facet_wrap(~ chr) +
      scale_colour_viridis_d(option = "plasma") +
      theme_bw()

dev.off()


