
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)

subtelomere_levels <- c("1p", "1q", "2p", "2q", "3p",
                        "3q", "4p", "4q", "5p", "5q",
                        "6p", "6q", "7q", "7p", "8p", "8q",
                        "9p", "9q", "10p", "10q", "11p", "11q",
                        "12p", "12q", "13q", "14q", "15q", "16p",
                        "16q", "17p", "17q", "18p", "18q", "19p",
                        "19q", "20p", "20q", "21q", "22q",
                        "XpYp", "Xp", "Xq", "Yp", "Yq")

                                       
tel_29_bp_frame <- read_delim("tel_29bp_blast_results.txt",
                              delim = "\t",
                              col_names=c("qseqid", "subtelomere", "pident",
                                          "length", "mismatch", "gapopen",
                                           "qstart", "qend",
                                           "start", "end",
                                           "evalue", "bitscore"))


tel_29_bp_frame <- tel_29_bp_frame %>%
                      separate(subtelomere,
                               into = c("subtelomere"),
                               extra = "drop",
                               sep = "tel_")

tel_29_bp_frame$subtelomere <- factor(tel_29_bp_frame$subtelomere,
                                      levels = subtelomere_levels,
                                      ordered = TRUE)

tel_29bp_0.05_evalue <- filter(tel_29_bp_frame,
                               evalue <= 0.05,
                               end <= 2000)


tel_29bp_three_or_less_mismatches <- filter(tel_29_bp_frame,
                                 length >= 29,
                                 pident >= ((29 - 3) / 29) * 100,
                                 end <= 2000)

rajika_seq_matches <- read_delim("rajika_sequence_blast.txt",
                              delim = "\t",
                              col_names=c("qseqid", "subtelomere", "pident",
                                          "length", "mismatch", "gapopen",
                                           "qstart", "qend",
                                           "start", "end",
                                           "evalue", "bitscore"))

rajika_seq_matches <- rajika_seq_matches %>%
                         separate(subtelomere,
                                  into = c("subtelomere"),
                                  extra = "drop",
                                  sep = "tel_")

rajika_seq_matches$subtelomere <- factor(rajika_seq_matches$subtelomere,
                                         levels = subtelomere_levels,
                                         ordered = TRUE)

rajika_seq_full_matches <- filter(rajika_seq_matches,
                                  length == 19,
                                  pident == 100.000,
                                  end <= 2000)



nergadze_info <- read_tsv("nergadze_tels_with_repeats.txt")

nergadze_info <- mutate(nergadze_info,
                        subtelomere = factor(subtelomere,
                                             levels = subtelomere_levels,
                                             ordered = TRUE))

tel_29bp_0.05_matches_per_subtelomere <- count(tel_29bp_0.05_evalue,
                                               subtelomere,
                                               name = "n_29bp_0.05_evalue",
                                               .drop = FALSE)


tel_29bp_full_matches_per_subtelomere <- count(tel_29bp_three_or_less_mismatches,
                                               subtelomere,
                                               name = "n_29bp_three_or_less_mismatches",
                                               .drop = FALSE)

rajika_seq_full_matches_per_subtelomere <- count(rajika_seq_full_matches,
                                                 subtelomere,
                                                 name = "n_rajika_full_matches",
                                                 .drop = FALSE)


all_info <- full_join(nergadze_info, tel_29bp_full_matches_per_subtelomere)

all_info <- full_join(all_info, tel_29bp_0.05_matches_per_subtelomere)

all_info <- full_join(all_info, rajika_seq_full_matches_per_subtelomere)


pdf("comparison_of_nergadze_and_blast.pdf")
                   
p <- all_info %>%
   filter(subtelomere != "XpYp") %>%
   gather(-subtelomere,
          key = "feature",
          value = "count") %>%
   mutate(chromosome = str_extract(subtelomere, "\\d+"),
          chromosome = ifelse(is.na(chromosome),
                              str_extract(subtelomere, "[[:upper:]]"),
                              chromosome),
          chromosome = factor(chromosome,
                              levels = c("Y", "X", as.character(22:1)),
                              ordered = TRUE),
          end = str_extract(subtelomere, "[[:lower:]]+"),
          feature = factor(feature,
                           levels = c("FISH", "nergadze_blast", "n_29bp_0.05_evalue",
                                      "n_29bp_three_or_less_mismatches", "n_rajika_full_matches"),
                           ordered = TRUE)) %>%
   ggplot(aes(x = end, y = chromosome, fill = count > 0)) +
      geom_tile(colour = "black") +
      geom_text(data = . %>% filter(feature %in% c("n_29bp_0.05_evalue",
                                                   "n_29bp_three_or_less_mismatches",
                                                   "n_rajika_full_matches")),
                mapping = aes(label = count)) +
      facet_wrap(~ feature,
                 nrow = 1,
                 labeller = as_labeller(c("FISH" = "Nergadze, et al\nFISH",
                                          "nergadze_blast" = "Nergadze, et al\nBLAST",
                                          "n_29bp_0.05_evalue" = "29bp repeat\ne-value <= 0.05",
                                          "n_29bp_three_or_less_mismatches" = "29bp repeat\n3 or less\nmismatches",
                                          "n_rajika_full_matches" = "TALE target\nexact matches"))) +
      scale_fill_brewer(name = "",
                        labels = c("Absent", "Present"),
                        type = "qual",
                        palette = "Set2") +
      theme_bw()

print(p)

dev.off()

