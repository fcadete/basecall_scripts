
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

tvr_files <- list.files(path = "telomere_hunter_results",
                        pattern = "*normalized_TVR_counts.tsv",
                        recursive = TRUE)
names(tvr_files) <- str_split_fixed(tvr_files, "/", 2)[,1]

tvr_frames <- lapply(as.list(paste0("telomere_hunter_results/", tvr_files)),
                     read_tsv)

tvr_frames <- lapply(1:length(tvr_frames),
                     function(i) bind_cols(tibble(sample = rep(names(tvr_files)[[i]], nrow(tvr_frames[[i]]))),
                                           tvr_frames[[i]]))

tvr_frames <- bind_rows(tvr_frames)

tcgj_type_tvrs <- c("TTAGGG", "TCAGGG", "TGAGGG", "TTGGGG")

ggplot(tvr_frames,
       aes(x = sample, y = Count_norm_by_all_reads_T, fill = Pattern)) +
   geom_col(position = "stack") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(tvr_frames,
       aes(x = sample, y = Count_norm_by_intratel_reads_T, fill = Pattern)) +
   geom_col(position = "stack") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(tvr_frames,
       aes(x = sample, y = Count_per_100_bp_intratel_read_T, fill = Pattern)) +
   geom_col(position = "stack") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(tvr_frames,
       aes(x = sample, y = Count_norm_by_all_reads_T, fill = Pattern %in% tcgj_type_tvrs)) +
   geom_col(position = "stack") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(tvr_frames,
       aes(x = sample, y = Count_norm_by_intratel_reads_T, fill = Pattern %in% tcgj_type_tvrs)) +
   geom_col(position = "stack") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(tvr_frames,
       aes(x = sample, y = Count_per_100_bp_intratel_read_T, fill = Pattern %in% tcgj_type_tvrs)) +
   geom_col(position = "stack") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

