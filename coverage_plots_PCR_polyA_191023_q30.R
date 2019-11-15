
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

rhietman_mapont_files <- list.files(path = "PCR_amplified/alignment_outputs/",
                                    pattern = "*rhietman_mapont_q30.depth")

rhietman_mapont_files <- rhietman_mapont_files[grepl("191023", rhietman_mapont_files)]

depth_frame <- data_frame()

for (file in rhietman_mapont_files) {
  

  this_file_depth <- read_tsv(paste0("PCR_amplified/alignment_outputs/", file),
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
                 mutate(barcode = str_extract(file, "BC\\d\\d"),
                        barcode = ifelse(is.na(barcode), "none", barcode))

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

rajika_seq_matches <- read_delim("rajika_sequence_blast.txt",
                              delim = "\t",
                              col_names=c("qseqid", "chr", "pident",
                                          "length", "mismatch", "gapopen",
                                           "qstart", "qend",
                                           "start", "end",
                                           "evalue", "bitscore"))

rajika_seq_matches <- rajika_seq_matches %>% separate(chr, into = c("chr"), extra = "drop")

rajika_seq_matches$chr <- factor(rajika_seq_matches$chr,
                          levels = subtelomere_levels,
                          ordered = TRUE)

rajika_seq_full_matches <- filter(rajika_seq_matches,
                                  length == 19,
                                  pident == 100.000)

#barcode_to_condition <- data.frame(barcode = c("BC01", "BC05", "BC09",
#                                               "BC02", "BC06", "BC10",
#                                               "BC03", "BC07", "BC11",
#                                               "BC04", "BC08", "BC12"),
#                                    sample_ID = c(rep("A", 3),
#                                               rep("B", 3),
#                                               rep("C", 3),
#                                               rep("D", 3)))

barcode_count <- depth_frame %>% count(barcode)

depth_frame_by_barcode <- left_join(data_frame(chr = rep(factor(rep(unique(depth_frame$chr), each = 500000),
                                                             levels = subtelomere_levels,
                                                             ordered = TRUE),
                                                      nrow(barcode_count)), 
                                            pos = rep(rep(1:500000, length(unique(depth_frame$chr)) ),
                                                      nrow(barcode_count)),
                                            barcode = rep(barcode_count$barcode,
                                                            each = (500000*length(unique(depth_frame$chr)))),
                                           ),
          depth_frame %>%
            group_by(barcode, chr, pos) %>%
            summarise(depth = sum(depth, na.rm = TRUE))) %>%
  mutate(depth = ifelse(is.na(depth), 0, depth))

total_reads <- read_tsv("guppy_VNP_polyA_sequencing_barcode_191023.porechop_finalcall_table") %>%
                 dplyr::count(final_barcode_call) %>%
                 select(barcode = final_barcode_call,
                        total_reads = n)

#subtelomere_end_reads <- read_delim("PCR_amplified/soft_clipped_end_counts.txt",
#                                  col_names = c("file", "reads"), delim = " ") %>%
#                             filter(str_detect(file, "191023")) %>%
#                             mutate(barcode = str_extract(file, "BC\\d\\d"),
#                                    barcode = ifelse(is.na(barcode), "none", barcode)) %>%
#                             group_by(barcode) %>%
#                             summarise(subtelomere_end_reads = sum(reads))
#
#full_join(total_reads, subtelomere_end_reads) %>%
#   select(barcode,
#          total_reads,
#          subtelomere_end_reads) %>%
#   write_tsv(path = "PCR_amplified/PCR_amplified_summary_table_191023.txt")


PCR_amplified_coverage_frame <- depth_frame_by_barcode %>%
                           drop_na() %>%
                           filter(pos <= 2000)
save(PCR_amplified_coverage_frame,
     file = "coverage_frames/TALES_coverage_frame_191023.RData")


pdf("PCR_amplified/coverage_plots/coverage_normalized_191023_q30.pdf", width = 16, height = 12)

depth_frame_by_barcode %>%
   drop_na() %>%
   filter(pos <= 2000) %>%
   ggplot() +
      geom_rect(data = rajika_seq_full_matches %>% filter(end <= 2000),
                mapping = aes(xmin = start,
                              xmax = end,
                              ymin = 0,
                              ymax = Inf,
                              fill = bitscore)) +
      scale_fill_gradient(low = "pink", high = "red") +
      geom_line(mapping = aes(x = pos,
                              y = depth,
                              colour = barcode),
                alpha = I(3/5)) +
      facet_wrap(~ chr) +
      scale_colour_viridis_d(option = "plasma") +
      theme_bw()


depth_frame_by_barcode %>%
   drop_na() %>%
   filter(pos <= 2000) %>%
   left_join(total_reads) %>%
   ggplot() +
      geom_rect(data = rajika_seq_full_matches %>% filter(end <= 2000),
                mapping = aes(xmin = start,
                              xmax = end,
                              ymin = 0,
                              ymax = Inf,
                              fill = bitscore)) +
      scale_fill_gradient(low = "pink", high = "red") +
      geom_line(mapping = aes(x = pos,
                              y = depth / total_reads,
                              colour = barcode),
                alpha = I(3/5)) +
      facet_wrap(~ chr) +
      scale_colour_viridis_d(option = "plasma") +
      theme_bw()

#depth_frame_by_barcode %>%
#   drop_na() %>%
#   filter(pos <= 2000) %>%
#   left_join(barcode_to_condition) %>%
#   ggplot() +
#      geom_rect(data = rajika_seq_full_matches %>% filter(end <= 2000),
#                mapping = aes(xmin = start,
#                              xmax = end,
#                              ymin = 0,
#                              ymax = Inf,
#                              fill = bitscore)) +
#      scale_fill_gradient(low = "pink", high = "red") +
#      geom_line(mapping = aes(x = pos,
#                              y = depth,
#                              group = barcode,
#                              colour = sample_ID),
#                alpha = I(3/5)) +
#      facet_wrap(~ chr) +
#      scale_colour_viridis_d(option = "plasma") +
#      theme_bw()
#
#
#depth_frame_by_barcode %>%
#   drop_na() %>%
#   filter(pos <= 2000) %>%
#   left_join(total_reads) %>%
#   left_join(barcode_to_condition) %>%
#   ggplot() +
#      geom_rect(data = rajika_seq_full_matches %>% filter(end <= 2000),
#                mapping = aes(xmin = start,
#                              xmax = end,
#                              ymin = 0,
#                              ymax = Inf,
#                              fill = bitscore)) +
#      scale_fill_gradient(low = "pink", high = "red") +
#      geom_line(mapping = aes(x = pos,
#                              y = depth / total_reads,
#                              group = barcode,
#                              colour = sample_ID),
#                alpha = I(3/5)) +
#      facet_wrap(~ chr) +
#      scale_colour_viridis_d(option = "plasma") +
#      theme_bw()
#
#
#depth_frame_by_barcode %>%
#   drop_na() %>%
#   filter(pos <= 2000) %>%
#   left_join(barcode_to_condition) %>%
#   ggplot() +
#      geom_rect(data = rajika_seq_full_matches %>% filter(end <= 2000),
#                mapping = aes(xmin = start,
#                              xmax = end,
#                              ymin = 0,
#                              ymax = Inf,
#                              fill = bitscore)) +
#      scale_fill_gradient(low = "pink", high = "red") +
#      geom_line(mapping = aes(x = pos,
#                              y = depth,
#                              group = barcode,
#                              colour = sample_ID),
#                alpha = I(3/5)) +
#      facet_wrap(~ chr, scales = "free_y") +
#      scale_colour_viridis_d(option = "plasma") +
#      theme_bw()
#
#
#depth_frame_by_barcode %>%
#   drop_na() %>%
#   filter(pos <= 2000) %>%
#   left_join(total_reads) %>%
#   left_join(barcode_to_condition) %>%
#   ggplot() +
#      geom_rect(data = rajika_seq_full_matches %>% filter(end <= 2000),
#                mapping = aes(xmin = start,
#                              xmax = end,
#                              ymin = 0,
#                              ymax = Inf,
#                              fill = bitscore)) +
#      scale_fill_gradient(low = "pink", high = "red") +
#      geom_line(mapping = aes(x = pos,
#                              y = depth / total_reads,
#                              group = barcode,
#                              colour = sample_ID),
#                alpha = I(3/5)) +
#      facet_wrap(~ chr, scales = "free_y") +
#      scale_colour_viridis_d(option = "plasma") +
#      theme_bw()
#
#PCR_amplified_coverage_frame %>%
#   group_by(sample_ID, chr, pos) %>%
#   summarise(mean = mean(depth), sd = sd(depth)) %>%
#   ggplot(aes(x = pos)) +
#      geom_ribbon(mapping = aes(ymin = mean - sd, ymax = mean + sd, fill = sample_ID),
#                  alpha = I(1/4)) +
#      geom_line(mapping = aes(y = mean, colour = sample_ID)) +
#      facet_wrap(~ chr) +
#      scale_fill_viridis_d(option = "plasma") +
#      scale_colour_viridis_d(option = "plasma") +
#      theme_bw()
#
#PCR_amplified_coverage_frame %>%
#   group_by(sample_ID, chr, pos) %>%
#   summarise(mean = mean(depth), sd = sd(depth)) %>%
#   ggplot(aes(x = pos)) +
#      geom_ribbon(mapping = aes(ymin = mean - sd, ymax = mean + sd, fill = sample_ID),
#                  alpha = I(1/4)) +
#      geom_line(mapping = aes(y = mean, colour = sample_ID)) +
#      facet_wrap(~ chr, scales = "free_y") +
#      scale_fill_viridis_d(option = "plasma") +
#      scale_colour_viridis_d(option = "plasma") +
#      theme_bw()
#
dev.off()


