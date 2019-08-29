
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

barcode_to_conditions <- read_tsv("190807_Seventh_run/Barcodes_clean.txt",
                                  col_names = c("cell_type", "date", "barcode")) %>%
                           mutate(barcode = str_replace(barcode, "NB", "BC"))


depth_frame <- data_frame()


for (barcode in barcode_to_conditions$barcode) {

   rhietman_mapont_files <- list.files(path = paste0("HeLa-HEK_barcoded/terra_primer_separated/", barcode),
                                       pattern = "*rhietman_mapont_primary.depth")
   
   
   for (file in rhietman_mapont_files) {
     
   
     this_file_depth <- read_tsv(paste0("HeLa-HEK_barcoded/terra_primer_separated/", barcode, "/", file),
                                         col_names = c("chr", "pos", "depth"),
                                         col_types = "cii")
   
     if (length(this_file_depth) > 0) {
     depth_frame <- rbind(depth_frame,
                          cbind(data_frame(file = file, barcode = barcode),
                                this_file_depth))
     }
   }
   
}

depth_frame <- as_data_frame(depth_frame)

depth_frame <- depth_frame %>% filter(pos <= 5000)

depth_frame <- depth_frame %>%
                 separate(chr, into = c("chr"), extra = "drop") %>%
                 mutate(file = str_remove(file, "_on_rhietman_mapont_primary.depth"))

depth_frame$chr <- factor(depth_frame$chr,
                          levels = subtelomere_levels,
                          ordered = TRUE)

depth_frame <- left_join(depth_frame, barcode_to_conditions)
 
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

barcode_count <- depth_frame %>% count(barcode)
primer_count <- depth_frame %>% count(file)

depth_frame_by_barcode <- left_join(data_frame(chr = rep(factor(rep(unique(depth_frame$chr), each = 5000),
                                                             levels = subtelomere_levels,
                                                             ordered = TRUE),
                                                      nrow(barcode_count)*nrow(primer_count)), 
                                            pos = rep(rep(1:5000, length(unique(depth_frame$chr)) ),
                                                      nrow(barcode_count)*nrow(primer_count)),
                                            primer = rep(rep(primer_count$file,
                                                         each = (5000*length(unique(depth_frame$chr)))), nrow(barcode_count)),
                                            barcode = rep(rep(barcode_count$barcode,
                                                           each = (5000*length(unique(depth_frame$chr)))), nrow(primer_count))
                                           ),
          depth_frame %>%
            filter(pos <= 5000) %>%
            group_by(barcode, file, chr, pos) %>%
            summarise(depth = sum(depth, na.rm = TRUE)) %>%
            select(barcode, primer = file, chr, pos, depth)) %>%
  mutate(depth = ifelse(is.na(depth), 0, depth))

#
#total_reads <- read_tsv("guppy_VNP_purified_TERRA_HeLa-HEK293_20190807.porechop_finalcall_table") %>%
#                 dplyr::count(final_barcode_call) %>%
#                 select(barcode = final_barcode_call,
#                        total_reads = n) %>%
#                 left_join(barcode_to_conditions)
#
#subtelomere_end_reads <- read_delim("HeLa-HEK_barcoded/soft_clipped_end_counts.txt",
#                                  col_names = c("file", "reads"), delim = " ") %>%
#                             mutate(barcode = str_extract(file, "BC\\d\\d"),
#                                    barcode = ifelse(is.na(barcode), "none", barcode)) %>%
#                             group_by(barcode) %>%
#                             summarise(subtelomere_end_reads = sum(reads))
#
#full_join(total_reads, subtelomere_end_reads) %>%
#   select(barcode, condition,
#          sample, total_reads,
#          subtelomere_end_reads) %>%
#   write_tsv(path = "HeLa-HEK_barcoded/HeLa-HEK_barcoded_summary_table.txt")
#

pdf("HeLa-HEK_barcoded/coverage_plots/coverage_by_primer.pdf", width = 16, height = 12)

depth_frame_by_barcode %>%
   drop_na() %>%
   filter(pos <= 2000) %>%
   left_join(barcode_to_conditions) %>%
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
                              group = paste(primer, barcode),
                              colour = cell_type)) +
      facet_wrap(~ chr + primer) +
      scale_colour_viridis_d(option = "plasma") +
      theme_bw()

dev.off()

