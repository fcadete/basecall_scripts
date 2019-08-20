
library(ggplot2)
library(tidyr)
library(readr)
library(dplyr)

full_call_table <- read_tsv("guppy_VNP_purified_TERRA_U2OS_20190514.porechop_finalcall_table",
                            col_names = c("readid", "runid",
                                          "sampleid", "read",
                                          "ch", "start_time", 
                                          "final_barcode_call"),
                            skip = 1)

full_call_table %>%
   arrange(final_barcode_call, start_time) %>%
   group_by(final_barcode_call) %>%
   mutate(cum_number_reads = row_number(final_barcode_call)) %>%
   ggplot(aes(x = start_time, y = cum_number_reads, colour = final_barcode_call)) +
     geom_line() +
     scale_colour_viridis_d() +
     theme_classic()

full_call_table %>%
   arrange(final_barcode_call, start_time) %>%
   group_by(final_barcode_call) %>%
   mutate(cum_number_reads = row_number(final_barcode_call)) %>%
   ggplot(aes(x = start_time, y = cum_number_reads, colour = final_barcode_call)) +
     geom_line() +
     scale_colour_viridis_d() +
     theme_classic() +
     scale_y_log10()

full_call_table %>%
   arrange(final_barcode_call, start_time) %>%
   group_by(final_barcode_call) %>%
   mutate(cum_number_reads = row_number(final_barcode_call),
          total_reads = max(cum_number_reads)) %>%
   ggplot(aes(x = start_time, y = cum_number_reads / total_reads, colour = final_barcode_call)) +
     geom_line() +
     scale_colour_viridis_d() +
     theme_classic()




