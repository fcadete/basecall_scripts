

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

CpG_matches <- vmatchPattern("GC", subtelomeres)

CpG_sliding_windows <- data_frame()

GC_content_sliding_windows <- data_frame()

for (subtelomere in names(subtelomeres)) {

      CpG_sliding_windows <- rbind(CpG_sliding_windows,
                                   data_frame(subtelomere = subtelomere,
                                              CpG_mean_over_window = roll_mean(as.vector(coverage(CpG_matches[[subtelomere]], width = 500000) / 2),
                                                                               n = 100, by = 10, fill = NA)))

      GC_content_sliding_windows <- rbind(GC_content_sliding_windows,
                                          data_frame(subtelomere = subtelomere,
                                                     GC_content_over_window = roll_mean(letterFrequencyInSlidingView(subtelomeres[[subtelomere]], 1, letters = c("G", "C")) %>%
                                                                                             rowSums(),
                                                                                        n = 100, by = 10, fill = NA)))

}

CpG_sliding_windows <- CpG_sliding_windows %>%
    group_by(subtelomere) %>%
    mutate(pos = 1:500000) %>%
    drop_na()

CpG_sliding_windows <- CpG_sliding_windows %>%
                 separate(subtelomere, into = c("subtelomere"), extra = "drop")

CpG_sliding_windows$subtelomere <- factor(CpG_sliding_windows$subtelomere,
                                          levels = subtelomere_levels,
                                          ordered = TRUE)


subtelomere_aligned_counts <- read_tsv("subtelomere_aligned_counts_all.tsv")

subtelomere_aligned_counts$seqnames <- factor(subtelomere_aligned_counts$seqnames,
                                              levels = subtelomere_levels,
                                              ordered = TRUE)

TSS_maxima <- read_tsv("telomere_distal_alignment_maxima.tsv")

TSS_maxima <- separate(TSS_maxima, long_seqnames, into = c("subtelomere"), extra = "drop")
TSS_maxima$subtelomere <- factor(TSS_maxima$subtelomere,
                                 levels = subtelomere_levels,
                                 ordered = TRUE)

png("coverage_plots/GpC_islands_zoomed.png", width = 1280, height = 1024)

ggplot(CpG_sliding_windows %>% filter(pos <= 50000),
       aes(x = pos, y = CpG_mean_over_window)) +
   geom_line() +
   facet_wrap(~ subtelomere) +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

png("coverage_plots/GpC_islands_extra_zoomed.png", width = 1280, height = 1024)

ggplot(CpG_sliding_windows %>% filter(pos <= 2000)) +
   geom_rect(data = left_join(TSS_maxima,
                              select(subtelomere_aligned_counts, subtelomere = seqnames, mapped_reads)),
             mapping = aes(xmin = start_coord, xmax = end_coord,
                           ymin = -Inf, ymax = Inf,
                           fill = mapped_reads),
             alpha = I(1/2)) +
   scale_fill_viridis_c() +
   geom_line(mapping = aes(x = pos, y = CpG_mean_over_window)) +
   facet_wrap(~ subtelomere) +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()



GC_content_sliding_windows <- GC_content_sliding_windows %>%
    group_by(subtelomere) %>%
    mutate(pos = 1:width(subtelomeres[subtelomere])) %>%
    drop_na()

GC_content_sliding_windows <- GC_content_sliding_windows %>%
                 separate(subtelomere, into = c("subtelomere"), extra = "drop")

GC_content_sliding_windows$subtelomere <- factor(GC_content_sliding_windows$subtelomere,
                          levels = subtelomere_levels,
                          ordered = TRUE)

png("coverage_plots/gc_content_zoomed.png", width = 1280, height = 1024)

ggplot(GC_content_sliding_windows %>% filter(pos <= 50000),
       aes(x = pos, y = GC_content_over_window)) +
   geom_line() +
   facet_wrap(~ subtelomere) +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


png("coverage_plots/gc_content_extra_zoomed.png", width = 1280, height = 1024)

ggplot(GC_content_sliding_windows %>% filter(pos <= 2000)) +
   geom_rect(data = left_join(TSS_maxima,
                              select(subtelomere_aligned_counts, subtelomere = seqnames, mapped_reads)),
             mapping = aes(xmin = start_coord, xmax = end_coord,
                           ymin = -Inf, ymax = Inf,
                           fill = mapped_reads),
             alpha = I(1/2)) +
   scale_fill_viridis_c() +
   geom_line(mapping = aes(x = pos, y = GC_content_over_window)) +
   facet_wrap(~ subtelomere) +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()




