
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(GenomicAlignments)
library(DESeq2)


subtelomere_levels <- c("1ptel", "1qtel", "2ptel", "2qtel", "3ptel",
                        "3qtel", "4ptel", "4qtel", "5ptel", "5qtel",
                        "6ptel", "6qtel", "7qtel", "7ptel", "8ptel", "8qtel",
                        "9ptel", "9qtel", "10ptel", "10qtel", "11ptel", "11qtel",
                        "12ptel", "12qtel", "13qtel", "14qtel", "15qtel", "16ptel",
                        "16qtel", "17ptel", "17qtel", "18ptel", "18qtel", "19ptel",
                        "19qtel", "20ptel", "20qtel", "21qtel", "22qtel",
                        "XpYptel", "Xqtel", "Yqtel")


barcode_to_condition <- read_tsv("190807_Seventh_run/Barcodes_clean.txt",
                                 col_names = c("cell_type", "date", "barcode")) %>%
                          mutate(barcode = str_replace(barcode, "NB", "BC"))

bam_files <- list.files(path = "HeLa-HEK_barcoded/soft_clipped_ends", pattern = "*bam")

subtelomere_bed <- read_tsv("subtelomere_starts_1_5000.bed",
                            col_names = c("seqnames", "start", "end"))

subtelomere_bed <- subtelomere_bed %>% mutate(end = 2000)

subtelomere_ranges <- GRanges(subtelomere_bed)

bam_records <- lapply(as.list(paste0("HeLa-HEK_barcoded/soft_clipped_ends/", bam_files)),
                      readGAlignments)

count_matrix <- do.call(cbind,
                        lapply(bam_records, function(ranges) countOverlaps(subtelomere_ranges,
                                                                           ranges)))

colnames(count_matrix) <- str_remove(bam_files, "_on_rhietman_mapont_primary_subtelomere_start.bam")
rownames(count_matrix) <- as.character(seqnames(subtelomere_ranges))


count_frame <- count_matrix %>%
                  as.data.frame() %>%
                  tibble::rownames_to_column(var = "subtelomere") %>%
                  gather(key = "file",
                         value = "counts",
                         -subtelomere)

file_to_samples <- tibble(file = colnames(count_matrix)) %>%
                      mutate(barcode = str_extract(file, "BC\\d\\d"),
                             barcode = ifelse(is.na(barcode), "none", barcode)) %>%
                      left_join(barcode_to_condition) %>%
                      mutate(sample = paste(cell_type, date, sep = "_"))

sample_count_frame <- left_join(count_frame, file_to_samples) %>%
                         group_by(subtelomere, sample) %>%
                         summarise(counts = sum(counts)) %>%
                         spread(key = sample, value = counts)


sample_count_matrix <- sample_count_frame %>%
                          ungroup() %>%
                          select(2:9) %>%
                          as.matrix()
rownames(sample_count_matrix) <- sample_count_frame$subtelomere

coldata <- file_to_samples %>%
              select(sample, cell_type, date) %>%
              unique() %>%
              drop_na() %>%
              as.data.frame()
rownames(coldata) <- coldata$sample
coldata <- coldata[colnames(sample_count_matrix), ]

dds <- DESeqDataSetFromMatrix(countData = sample_count_matrix,
                              colData = coldata,
                              design = ~ cell_type)

dds <- DESeq(dds)

results(dds, name = "cell_type_HeLa_vs_HEK293T") %>%
   as_tibble(rownames = "subtelomere") %>%
   filter(padj <= 0.05)


ddsLRT_minusDate <- DESeqDataSetFromMatrix(countData = sample_count_matrix,
                                           colData = coldata,
                                           design = ~ cell_type + date)

ddsLRT_minusDate <- DESeq(ddsLRT_minusDate,
                          test = "LRT",
                          reduced = ~ date)

results(ddsLRT_minusDate) %>%
   as_tibble(rownames = "subtelomere") %>%
   filter(padj <= 0.05)



ddsLRT_minusCell <- DESeqDataSetFromMatrix(countData = sample_count_matrix,
                                           colData = coldata,
                                           design = ~ cell_type + date)

ddsLRT_minusCell <- DESeq(ddsLRT_minusCell,
                          test = "LRT",
                          reduced = ~ cell_type)

results(ddsLRT_minusCell) %>%
   as_tibble(rownames = "subtelomere") %>%
   filter(padj <= 0.05)



sample_count_frame <- sample_count_frame %>%
                         separate(subtelomere, into = c("subtelomere"), extra = "drop")

sample_count_frame$subtelomere <- factor(sample_count_frame$subtelomere,
                          levels = subtelomere_levels,
                          ordered = TRUE)

pdf("HeLa-HEK_barcoded/level_comparison_plots.pdf", width = 16, height = 12)

sample_count_frame %>%
   gather(key = "sample", value = "counts", 2:9) %>%
   left_join(file_to_samples) %>%
   mutate(subtelomere = factor(subtelomere,
                          levels = subtelomere_levels,
                          ordered = TRUE) ) %>%
   ggplot(aes(x = sample, y = counts, colour = cell_type)) +
      geom_point() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

sample_count_frame %>%
   gather(key = "sample", value = "counts", 2:9) %>%
   left_join(file_to_samples) %>%
   mutate(subtelomere = factor(subtelomere,
                          levels = subtelomere_levels,
                          ordered = TRUE) ) %>%
   ggplot(aes(x = sample, y = counts, colour = cell_type)) +
      geom_jitter(height = 0) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

sample_count_frame %>%
   gather(key = "sample", value = "counts", 2:9) %>%
   left_join(file_to_samples) %>%
   mutate(subtelomere = factor(subtelomere,
                          levels = subtelomere_levels,
                          ordered = TRUE) ) %>%
   ggplot(aes(x = sample, y = counts, colour = cell_type)) +
      geom_point() +
      facet_wrap(~ subtelomere, scales = "free_y") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

sample_count_frame %>%
   gather(key = "sample", value = "counts", 2:9) %>%
   left_join(file_to_samples) %>%
   mutate(subtelomere = factor(subtelomere,
                          levels = subtelomere_levels,
                          ordered = TRUE) ) %>%
   ggplot(aes(x = sample, y = counts, colour = cell_type)) +
      geom_point() +
      facet_wrap(~ subtelomere) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))


sample_count_frame %>%
   gather(key = "sample", value = "counts", 2:9) %>%
   left_join(file_to_samples) %>%
   mutate(subtelomere = factor(subtelomere,
                          levels = subtelomere_levels,
                          ordered = TRUE) ) %>%
   ggplot(aes(x = sample, y = counts, colour = date)) +
      geom_point() +
      facet_wrap(~ subtelomere) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))



dev.off()

