

library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Biostrings)
library(RcppRoll)

subtelomere_levels <- c("1ptel", "1qtel", "2ptel", "2qtel", "3ptel",
                                     "3qtel", "4ptel", "4qtel", "5ptel", "5qtel",
                                     "6ptel", "6qtel", "7qtel", "7ptel", "8ptel", "8qtel",
                                     "9ptel", "9qtel", "10ptel", "10qtel", "11ptel", "11qtel",
                                     "12ptel", "12qtel", "13qtel", "14qtel", "15qtel", "16ptel",
                                     "16qtel", "17ptel", "17qtel", "18ptel", "18qtel", "19ptel",
                                     "19qtel", "20ptel", "20qtel", "21qtel", "22qtel",
                                     "XpYptel", "Xqtel", "Yqtel")

subtelomeres <- readDNAStringSet("references/ConcatenatedFASTAAassemblies_hTel.txt", "fasta")

CpG_matches <- vmatchPattern("CG", subtelomeres)

CpG_islands <- data_frame()

for (subtelomere in names(subtelomeres)) {

      G_and_C_frequency_in_subtelomere <- letterFrequencyInSlidingView(subtelomeres[[subtelomere]], 200, c("G", "C"))

      CpG_expected_in_subtelomere <- (G_and_C_frequency_in_subtelomere[,"C"] * G_and_C_frequency_in_subtelomere[,"G"]) / 200

      GC_content_in_subtelomere <- (G_and_C_frequency_in_subtelomere[,"C"] + G_and_C_frequency_in_subtelomere[,"G"]) / 200


      two_hundred_bp_windows <- GRanges(seqnames = subtelomere,
                                        ranges = IRanges(start = seq(200, length(subtelomeres[[subtelomere]]), by = 1) - 199,
                                        end = seq(200, length(subtelomeres[[subtelomere]]), by = 1)),
                                        strand = "*")

      CpG_observed_in_subtelomere <- countOverlaps(two_hundred_bp_windows,
                                                   GRanges(seqnames = subtelomere, ranges =  CpG_matches[[subtelomere]]))

      CpG_islands_in_subtelomere <- ((CpG_observed_in_subtelomere / CpG_expected_in_subtelomere) > 0.6) &
                                      (GC_content_in_subtelomere > 0.5)

      CpG_islands_in_subtelomere[is.na(CpG_islands_in_subtelomere)] <- FALSE 

      two_hundred_bp_windows$CpG_islands <- CpG_islands_in_subtelomere

      CpG_islands <- bind_rows(CpG_islands,
                               as_tibble(reduce(two_hundred_bp_windows[two_hundred_bp_windows$CpG_islands])))
}


CpG_islands <- CpG_islands %>%
                  separate(seqnames, into = c("subtelomere"), extra = "drop") %>%
                  mutate(subtelomere = factor(subtelomere,
                                              levels = c(subtelomere_levels),
                                              ordered = TRUE)) %>%
                  arrange(subtelomere)

write_tsv(CpG_islands,
          path = "CpG_islands_in_subtelomeres.tsv")


