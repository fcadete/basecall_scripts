
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

rhietman_mapont_files <- list.files(path = "alignment_outputs/",
                                    pattern = "*rhietman_mapont_primary.depth")

depth_frame <- data_frame()

for (file in rhietman_mapont_files) {
  
  if (grepl("Joana_TERRA_VNP", file)) {

  this_file_depth <- read_tsv(paste0("alignment_outputs/", file),
                                      col_names = c("chr", "pos", "depth"),
                                      col_types = "cii")

  if (length(this_file_depth) > 0) {
  depth_frame <- rbind(depth_frame,
                       cbind(data_frame(file = file),
                             this_file_depth))
  }
  }
}

depth_frame <- as_data_frame(depth_frame)

depth_frame <- depth_frame %>%
                 separate(chr, into = c("chr"), extra = "drop") %>%
                 separate(file, sep = "_Joana_TERRA_VNP", into = c("file"), extra = "drop")

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



for (one_file in unique(depth_frame$file)) {
   if (nrow(depth_frame %>% filter(file == one_file, pos <= 2000))>0) {
         png(paste0("coverage_plots/by_sample/", one_file, ".png"), width = 1280, height = 1024)
         p <- depth_frame %>%
                 filter(file == one_file, pos <= 2000) %>%
                 ggplot(aes(x = pos, y = depth)) +
                    geom_line() +
                    facet_wrap(~ chr) +
                    scale_x_continuous(name = "Distance from telomere repeats (bp)",
                                       breaks = c(1, 500, 1000, 1500, 2000), 
                                       limits = c(1, 2000)) +
                    theme_bw() +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                    labs(title=one_file)
         print(p)
         dev.off()
    }
}

#
#for (this_file in unique(depth_frame$file)) {
#  ggsave(paste0("coverage_plots/coverage_", this_file, ".png"),
#         depth_frame %>% filter(file == this_file) %>%
#          group_by(chr, pos) %>%
#          summarise(depth = sum(depth, na.rm = TRUE)) %>%
#          ggplot(aes(x = pos, y = depth)) +
#          geom_line() +
#          facet_wrap(~ chr) +
#          labs(main = this_file) +
#          scale_x_continuous(name = this_file) +
#          theme_bw(),
#      width = 20, height = 15, units = "cm", device = png())
#
#ggsave(paste0("coverage_plots/coverage_", this_file, "_zoomed.png"),
#       left_join(data_frame(chr = factor(rep(unique(depth_frame$chr), each = 10000),
#                                          levels = subtelomere_levels,
#                                           ordered = TRUE),
#                               pos = rep(1:10000, length(unique(depth_frame$chr)) )),
#                  depth_frame %>% filter(file == this_file) %>% filter(pos <= 10000) %>%
#                    group_by(chr, pos) %>%
#                    summarise(depth = sum(depth, na.rm = TRUE)) %>% select(chr, pos, depth)) %>%
#          mutate(depth = ifelse(is.na(depth), 0, depth)) %>%
#          ggplot() +
#          geom_line(mapping = aes(x = pos, y = depth)) +
#          geom_rect(data = tel_29_bp_frame,
#                    mapping = aes(xmin = start,
#                                  xmax = end,
#                                  ymin = 0,
#                                  ymax = Inf,
#                                  fill = bitscore),
#                    alpha = I(1/2)) +
#          scale_fill_viridis_c() +
#       facet_wrap(~ chr) +
#          labs(main = this_file) +
#          scale_x_continuous(name = this_file,
#                     breaks = c(1, 2000, 4000, 6000, 8000, 10000),
#                     limits = c(1, 10000)) +
#  theme_bw() +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)),
#      width = 20, height = 15, units = "cm", device = png())
#
#
#}
#

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


other_coverage_frame <- depth_frame_by_file %>%
                           drop_na() %>%
                           filter(pos <= 2000) %>%
                           left_join(total_reads)
save(other_coverage_frame,
     file = "coverage_frames/other_coverage_frame.RData")

pdf("coverage_plots/coverage_normalized.pdf", width = 16, height = 12)

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


for (this_cell_type in unique(depth_frame$cell_type)) {

   png(paste0("coverage_plots/coverage_pooled_", this_cell_type, ".png"),
       width = 1280, height = 1024)
   
   depth_frame %>%
     filter(cell_type == this_cell_type) %>%
     ggplot() +
     geom_rect(data = tel_29_bp_frame,
               mapping = aes(xmin = start,
                             xmax = end,
                             ymin = 0,
                             ymax = Inf,
                             fill = bitscore),
               alpha = I(1/2)) +
      scale_fill_viridis_c() +
     geom_line(mapping = aes(x = pos, y = depth)) +
     facet_wrap( ~ chr) +
     scale_x_continuous(name = "Distance from telomere repeats (bp)",
                        breaks = c(1, 100000, 200000, 300000, 400000, 500000)) +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     labs(title = this_cell_type)
   
   dev.off()
   
   png(paste0("coverage_plots/coverage_pooled_", this_cell_type, "_zoomed.png"),
       width = 1280, height = 1024)
   
   depth_frame %>%
     filter(cell_type == this_cell_type) %>%
     filter(pos <= 10000) %>%
     droplevels() %>%
     ggplot() +
     geom_rect(data = tel_29_bp_frame,
               mapping = aes(xmin = start,
                             xmax = end,
                             ymin = 0,
                             ymax = Inf,
                             fill = bitscore),
               alpha = I(1/2)) +
     scale_fill_viridis_c() +
     geom_line(mapping = aes(x = pos, y = depth)) +
     facet_wrap( ~ chr) +
     scale_x_continuous(name = "Distance from telomere repeats (bp)",
                        breaks = c(1, 2000, 4000, 6000, 8000, 10000),
                        limits = c(1, 10000)) +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     labs(title = this_cell_type)
   
   dev.off()
   
   png(paste0("coverage_plots/coverage_pooled_", this_cell_type, "_extra_zoomed.png"),
       width = 1280, height = 1024)
   
   depth_frame %>%
     filter(cell_type == this_cell_type) %>%
     filter(pos <= 2000) %>%
     droplevels() %>%
     ggplot() +
     geom_rect(data = tel_29_bp_frame,
               mapping = aes(xmin = start,
                             xmax = end,
                             ymin = 0,
                             ymax = Inf,
                             fill = bitscore),
               alpha = I(1/2)) +
     scale_fill_viridis_c() +
     geom_line(mapping = aes(x = pos, y = depth)) +
     facet_wrap( ~ chr) +
     scale_x_continuous(name = "Distance from telomere repeats (bp)",
                        breaks = c(1, 500, 1000, 1500, 2000),
                        limits = c(1, 2000)) +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     labs(title = this_cell_type)
   
   dev.off()
    
   png(paste0("coverage_plots/coverage_pooled_", this_cell_type, "_extra_zoomed_rajikas_targets.png"),
       width = 1280, height = 1024)
   
   depth_frame %>%
     filter(cell_type == this_cell_type) %>%
     filter(pos <= 2000) %>%
     droplevels() %>%
     ggplot() +
     geom_rect(data = rajika_seq_matches,
               mapping = aes(xmin = start,
                             xmax = end,
                             ymin = 0,
                             ymax = Inf,
                             fill = bitscore),
               alpha = I(1/2)) +
     scale_fill_viridis_c() +
     geom_line(mapping = aes(x = pos, y = depth)) +
     facet_wrap( ~ chr) +
     scale_x_continuous(name = "Distance from telomere repeats (bp)",
                        breaks = c(1, 500, 1000, 1500, 2000),
                        limits = c(1, 2000)) +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     labs(title = this_cell_type)
   
   dev.off()
    
   png(paste0("coverage_plots/coverage_pooled_", this_cell_type, "_extra_zoomed_rajikas_full_match_targets.png"),
       width = 1280, height = 1024)
   
   depth_frame %>%
     filter(cell_type == this_cell_type) %>%
     filter(pos <= 2000) %>%
     droplevels() %>%
     ggplot() +
     geom_rect(data = rajika_seq_full_matches,
               mapping = aes(xmin = start,
                             xmax = end,
                             ymin = 0,
                             ymax = Inf,
                             fill = bitscore),
               alpha = I(1/2)) +
     scale_fill_viridis_c() +
     geom_line(mapping = aes(x = pos, y = depth)) +
     facet_wrap( ~ chr) +
     scale_x_continuous(name = "Distance from telomere repeats (bp)",
                        breaks = c(1, 500, 1000, 1500, 2000),
                        limits = c(1, 2000)) +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     labs(title = this_cell_type)
   
   dev.off()
   
  
  
}



