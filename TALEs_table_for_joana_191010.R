
library(tidyverse)


barcode_to_condition <- data.frame(barcode = c("BC01", "BC05", "BC09",
                                               "BC02", "BC06", "BC10",
                                               "BC03", "BC07", "BC11",
                                               "BC04", "BC08", "BC12"),
                                    sample_ID = c(rep("A", 3),
                                                  rep("B", 3),
                                                  rep("C", 3),
                                                  rep("D", 3)),
                                    replicate = rep(c("1", "2", "3"),
                                                    4),
                                    condition = c(rep("SID4_dox+", 3),
                                                  rep("SID4_dox-", 3),
                                                  rep("NLS3_dox+", 3),
                                                  rep("NLS3_dox-", 3)))


TALEs_barcode_counts <- read_tsv("guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010.porechop_finalcall_table") %>%
                           dplyr::count(final_barcode_call) %>%
                           select(barcode = final_barcode_call, total_reads = n) %>%
                           filter(barcode != "none") %>%
                           left_join(barcode_to_condition)


primer_call_files <- list.files("TALEs/primer_separation/", "porechop_finalcall_table", full.names = TRUE)

primer_counts <- primer_call_files %>%
                   map(~ mutate(read_tsv(.),
                                barcode = str_extract(., "BC\\d\\d"),
                                primer = final_barcode_call)) %>%
                   reduce(rbind) %>%
                   count(barcode, primer) %>%
                   filter(grepl("Joana_TERRA", primer)) %>%
                   group_by(barcode) %>%
                   summarise(n = sum(n)) %>%
                   select(barcode, primered_reads = n)

subtelomere_end_reads <- read_delim("TALEs/soft_clipped_end_primered_counts.txt",
                                    col_names = c("file",
                                    "subtelomere_end_reads"), delim = " ") %>%
        filter(grepl("Joana_TERRA", file)) %>%
        mutate(file = str_remove(file, "_Joana_TERRA_VNP_no_PCR"),
               file = str_remove(file, "_Joana_TERRA_VNP"),
               file = str_remove(file, "_on_rhietman_mapont_primary_subtelomere_start.bam")) %>%
        group_by(file) %>%
        summarise(subtelomere_end_reads = sum(subtelomere_end_reads)) %>%
        mutate(barcode = str_extract(file, "BC\\d\\d"))
        

left_join(TALEs_barcode_counts, primer_counts) %>%
   left_join(subtelomere_end_reads) %>%
   select(barcode, sample_ID, replicate, condition, total_reads, primered_reads, subtelomere_end_reads) %>%
   write_tsv(path = "TALEs/191010_primered_table.tsv")

