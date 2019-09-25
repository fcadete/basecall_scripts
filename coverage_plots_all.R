
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
                 separate(file, sep = "_Joana_TERRA_VNP", into = c("file"), extra = "drop") %>%
                 filter(grepl("purified", file),
                        (grepl("TALE", file)) == FALSE) 

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


depth_frame <- depth_frame %>%
     mutate(cell_type = ifelse(grepl("HeL", file),
                               "HeLa",
                               ifelse(grepl("GM847", file),
                                      "GM847",
                                      ifelse(grepl("SAOS2", file),
                                             "SAOS2",
                                             ifelse(grepl("HEK", file),
                                                    "HEK293T",
                                                    "U2OS")))))

cell_type_file_count <- depth_frame %>% count(cell_type, file)

depth_frame_by_file <- left_join(data_frame(chr = rep(factor(rep(unique(depth_frame$chr), each = 5000),
                                                             levels = subtelomere_levels,
                                                             ordered = TRUE),
                                                      nrow(cell_type_file_count)), 
                                            pos = rep(rep(1:5000, length(unique(depth_frame$chr)) ),
                                                      nrow(cell_type_file_count)),
                                            cell_type = rep(cell_type_file_count$cell_type,
                                                            each = (5000*length(unique(depth_frame$chr)))),
                                            file = rep(cell_type_file_count$file,
                                                       each = (5000*length(unique(depth_frame$chr))))),
          depth_frame %>%
            filter(pos <= 5000) %>%
            group_by(cell_type, file, chr, pos) %>%
            summarise(depth = sum(depth, na.rm = TRUE))) %>%
  mutate(depth = ifelse(is.na(depth), 0, depth))

system("wc -l *porechop_finalcall_table > total_reads_per_file.txt")
total_reads <- read_table("total_reads_per_file.txt",
                          col_names = c("total_reads", "file")) %>%
    filter(file != "total") %>% 
    mutate(file = str_remove(file, "guppy_"),
           file = str_remove(file, ".porechop_finalcall_table"))

total_reads_by_cell_type <- total_reads %>%
     filter(grepl("purified", file),
            (grepl("TALE", file)) == FALSE) %>%
     mutate(cell_type = ifelse(grepl("HeL", file),
                               "HeLa",
                               ifelse(grepl("GM847", file),
                                      "GM847",
                                      ifelse(grepl("SAOS2", file),
                                             "SAOS2",
                                             ifelse(grepl("HEK", file),
                                                    "HEK293T",
                                                    "U2OS"))))) %>%
     group_by(cell_type) %>%
     summarise(total_reads = sum(total_reads))




barcode_to_conditions <- read_tsv("190807_Seventh_run/Barcodes_clean.txt",
                                  col_names = c("cell_type", "date", "barcode")) %>%
                           mutate(barcode = str_replace(barcode, "NB", "BC"))

depth_frame_barcoded_HH <- data_frame()

for (barcode in barcode_to_conditions$barcode) {

   rhietman_mapont_files <- list.files(path = paste0("HeLa-HEK_barcoded/terra_primer_separated/", barcode),
                                       pattern = "*rhietman_mapont_primary.depth")
   
   
   for (file in rhietman_mapont_files) {
     
   
     this_file_depth <- read_tsv(paste0("HeLa-HEK_barcoded/terra_primer_separated/", barcode, "/", file),
                                         col_names = c("chr", "pos", "depth"),
                                         col_types = "cii")
   
     if (length(this_file_depth) > 0) {
     depth_frame_barcoded_HH <- rbind(depth_frame_barcoded_HH,
                          cbind(data_frame(file = file, barcode = barcode),
                                this_file_depth))
     }
   }
   
}

depth_frame_barcoded_HH <- as_data_frame(depth_frame_barcoded_HH)

depth_frame_barcoded_HH <- depth_frame_barcoded_HH %>% filter(pos <= 5000)

depth_frame_barcoded_HH <- depth_frame_barcoded_HH %>%
                 separate(chr, into = c("chr"), extra = "drop") %>%
                 mutate(file = str_remove(file, "_on_rhietman_mapont_primary.depth"))

depth_frame_barcoded_HH$chr <- factor(depth_frame_barcoded_HH$chr,
                          levels = subtelomere_levels,
                          ordered = TRUE)

depth_frame_barcoded_HH <- left_join(depth_frame_barcoded_HH, barcode_to_conditions)
 
depth_frame_barcoded_HH <- filter(depth_frame_barcoded_HH,
                                  str_detect(file, "Joana_TERRA_VNP"))


barcode_count <- depth_frame_barcoded_HH %>% count(barcode)
primer_count <- depth_frame_barcoded_HH %>% count(file)

depth_frame_barcoded_HH_by_barcode <- left_join(data_frame(chr = rep(factor(rep(unique(depth_frame_barcoded_HH$chr), each = 5000),
                                                             levels = subtelomere_levels,
                                                             ordered = TRUE),
                                                      nrow(barcode_count)*nrow(primer_count)), 
                                            pos = rep(rep(1:5000, length(unique(depth_frame_barcoded_HH$chr)) ),
                                                      nrow(barcode_count)*nrow(primer_count)),
                                            primer = rep(rep(primer_count$file,
                                                         each = (5000*length(unique(depth_frame_barcoded_HH$chr)))), nrow(barcode_count)),
                                            barcode = rep(rep(barcode_count$barcode,
                                                           each = (5000*length(unique(depth_frame_barcoded_HH$chr)))), nrow(primer_count))
                                           ),
          depth_frame_barcoded_HH %>%
            filter(pos <= 5000) %>%
            group_by(barcode, file, chr, pos) %>%
            summarise(depth = sum(depth, na.rm = TRUE)) %>%
            select(barcode, primer = file, chr, pos, depth)) %>%
  mutate(depth = ifelse(is.na(depth), 0, depth))

depth_frame_barcoded_HH_by_barcode <- left_join(depth_frame_barcoded_HH_by_barcode,
                                                barcode_to_conditions)

total_reads_barcoded_HH <- read_tsv("guppy_VNP_purified_TERRA_HeLa-HEK293_20190807.porechop_finalcall_table") %>%
                 dplyr::count(final_barcode_call) %>%
                 select(barcode = final_barcode_call,
                        total_reads_barcoded_HH = n) %>%
                 left_join(barcode_to_conditions)


depth_frame_all <- rbind(select(depth_frame_by_file,
                                chr, pos, cell_type, depth),
                         select(depth_frame_barcoded_HH_by_barcode,
                                chr, pos, cell_type, depth))

total_reads_all <- rbind(total_reads_by_cell_type,
                         drop_na(select(total_reads_barcoded_HH,
                                        cell_type,
                                        total_reads = total_reads_barcoded_HH))) %>%
                      group_by(cell_type) %>%
                      summarise(total_reads = sum(total_reads))


pdf("coverage_primed_reads_all_samples_pooled.pdf",
    width = 16,
    height = 12)

depth_frame_all %>%
   group_by(cell_type, chr, pos) %>%
   summarise(depth = sum(depth)) %>%
   ggplot(aes(x = pos, y = depth, group = cell_type, colour = cell_type)) + 
      geom_line() +
      facet_wrap(~ chr)

depth_frame_all %>%
   group_by(cell_type, chr, pos) %>%
   summarise(depth = sum(depth)) %>%
   left_join(total_reads_all) %>%
   ggplot(aes(x = pos, y = depth / total_reads, group = cell_type, colour = cell_type)) + 
      geom_line() +
      facet_wrap(~ chr)

dev.off()


depth_frame_files <- rbind(select(depth_frame_by_file,
                                  chr, pos, cell_type, file, depth),
                           mutate(depth_frame_barcoded_HH_by_barcode,
                                  file = paste(cell_type, date, sep = "_")) %>%
                              select(chr, pos, cell_type, file, depth))



total_reads_by_file <- total_reads %>%
     filter(grepl("purified", file),
            (grepl("TALE", file)) == FALSE) %>%
     mutate(cell_type = ifelse(grepl("HeL", file),
                               "HeLa",
                               ifelse(grepl("GM847", file),
                                      "GM847",
                                      ifelse(grepl("SAOS2", file),
                                             "SAOS2",
                                             ifelse(grepl("HEK", file),
                                                    "HEK293T",
                                                    "U2OS"))))) %>%
   select(file, total_reads)

total_reads_by_sample <- rbind(total_reads_by_file,
                               mutate(total_reads_barcoded_HH,
                                      file = paste(cell_type, date, sep = "_")) %>%
                                 drop_na() %>%
                                 select(file,
                                        total_reads = total_reads_barcoded_HH))


for (this_cell_type in unique(depth_frame_files$cell_type)) {

   n_files <- length(unique(filter(depth_frame_files, cell_type == this_cell_type)$file))

   p <- depth_frame_files %>%
         filter(cell_type == this_cell_type) %>%
         group_by(cell_type, file, chr, pos) %>%
         summarise(depth = sum(depth)) %>%
         ggplot(aes(x = pos, y = depth, group = file)) + 
            geom_line(alpha = I(1/n_files)) +
            facet_wrap(~ chr) +
            labs(title = this_cell_type)
   print(p)

   p <- depth_frame_files %>%
           filter(cell_type == this_cell_type) %>%
           group_by(cell_type, file, chr, pos) %>%
           summarise(depth = sum(depth)) %>%
           left_join(total_reads_by_sample) %>%
           ggplot(aes(x = pos, y = depth / total_reads, group = file)) +
              geom_line(alpha = I(1/n_files)) +
              facet_wrap(~ chr) +
              labs(title = this_cell_type)
   print(p)

}



