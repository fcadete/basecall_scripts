
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(GenomicAlignments)
library(DESeq2)
library(ggrepel)

pdf("TALEs/TALES_DESeq2_191010.pdf")

theme_set(theme_bw())

subtelomere_levels <- c("1ptel", "1qtel", "2ptel", "2qtel", "3ptel",
                        "3qtel", "4ptel", "4qtel", "5ptel", "5qtel",
                        "6ptel", "6qtel", "7qtel", "7ptel", "8ptel", "8qtel",
                        "9ptel", "9qtel", "10ptel", "10qtel", "11ptel", "11qtel",
                        "12ptel", "12qtel", "13qtel", "14qtel", "15qtel", "16ptel",
                        "16qtel", "17ptel", "17qtel", "18ptel", "18qtel", "19ptel",
                        "19qtel", "20ptel", "20qtel", "21qtel", "22qtel",
                        "XpYptel", "Xqtel", "Yqtel")


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

bam_files <- list.files(path = "TALEs/soft_clipped_ends", pattern = "*bam")

bam_files <- bam_files[grepl("191010", bam_files)]

subtelomere_bed <- read_tsv("subtelomere_starts_1_5000.bed",
                            col_names = c("seqnames", "start", "end"))

subtelomere_bed <- subtelomere_bed %>% mutate(end = 2000)

subtelomere_ranges <- GRanges(subtelomere_bed)

bam_records <- lapply(as.list(paste0("TALEs/soft_clipped_ends/", bam_files)),
                      readGAlignments)

count_matrix <- do.call(cbind,
                        lapply(bam_records, function(ranges) countOverlaps(subtelomere_ranges, ranges)))

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
                      left_join(barcode_to_condition)

sample_count_frame <- left_join(count_frame, file_to_samples)


coldata <- file_to_samples %>% select(file, barcode, replicate, sample_ID) %>% unique() %>% drop_na() %>% as.data.frame()
rownames(coldata) <- coldata$file

count_matrix <- count_matrix[, rownames(coldata)]


dds_trn <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ sample_ID)

TALEs_barcode_counts <- read_tsv("guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010.porechop_finalcall_table") %>%
                           dplyr::count(final_barcode_call) %>%
                           select(barcode = final_barcode_call, total_reads = n) %>%
                           filter(barcode != "none") %>%
                           left_join(file_to_samples)

total_reads <- TALEs_barcode_counts$total_reads
names(total_reads) <- TALEs_barcode_counts$file


sizeFactors(dds_trn) <- total_reads[colnames(count_matrix)] / mean(total_reads[colnames(count_matrix)])

rlogs_trn_blind <- rlog(dds_trn, blind = TRUE)

p <- plotPCA(rlogs_trn_blind, intgroup = "sample_ID") +
   scale_colour_viridis_d() +
   labs(title = "Blind PCA")
print(p)

p <- plotPCA(rlogs_trn_blind, intgroup = "replicate") +
   scale_colour_viridis_d() +
   labs(title = "Blind PCA")
print(p)

rlogs_trn_design <- rlog(dds_trn, blind = FALSE)

p <- plotPCA(rlogs_trn_design, intgroup = "sample_ID") +
   scale_colour_viridis_d() +
   labs(title = "Non-blind PCA")
print(p)

p <- plotPCA(rlogs_trn_design, intgroup = "replicate") +
   scale_colour_viridis_d() +
   labs(title = "Non-blind PCA")
print(p)

rlogs_trn_design_frame <- as.data.frame(assay(rlogs_trn_design))
rlogs_trn_design_frame$subtelomere <- rownames(rlogs_trn_design_frame)
rlogs_trn_design_frame <- rlogs_trn_design_frame %>% gather(key = "file", value = "rlogs_trn_design", -subtelomere)
rlogs_trn_design_frame <- left_join(rlogs_trn_design_frame, file_to_samples)

rlogs_trn_design_frame <- rlogs_trn_design_frame %>%
                            separate(subtelomere, into = "subtelomere", extra = "drop") %>%
                            mutate(subtelomere = factor(subtelomere,
                                                        levels = subtelomere_levels,
                                                        ordered = TRUE))

p <- rlogs_trn_design_frame %>%
   ggplot(aes(y = rlogs_trn_design,
              x = sample_ID,
              colour = sample_ID)) +
   geom_point() +
   facet_wrap(~ subtelomere) +
   scale_colour_viridis_d() +
   labs(title = "Non-blind rlogs")
print(p)

p <- rlogs_trn_design_frame %>%
   ggplot(aes(fill = rlogs_trn_design,
              x = barcode,
              y = subtelomere)) +
   geom_tile() +
   scale_fill_viridis_c() +
   facet_wrap(~ sample_ID, nrow = 1, scales = "free_x") +
   labs(title = "Non-blind rlogs")
print(p)


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

p <- rlogs_trn_design_frame %>%
   ggplot(aes(x = sample_ID, y = rlogs_trn_design, colour = sample_ID,
              shape = subtelomere %in% filter(rajika_seq_full_matches, end <= 2000)$chr)) +
      geom_point() +
      scale_shape_discrete("Has TALE target") +
      facet_wrap(~ subtelomere) +
      scale_colour_viridis_d() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)

dds_trn <- DESeq(dds_trn)

results_A_B <- results(dds_trn,
                       contrast = c("sample_ID", "A", "B")) %>%
                  as_tibble(rownames = c("subtelomere")) %>%
                  separate(subtelomere, into = "subtelomere", extra = "drop") %>%
                  mutate(subtelomere = factor(subtelomere,
                                              levels = subtelomere_levels,
                                              ordered = TRUE),
                         has_rajika = subtelomere %in% filter(rajika_seq_full_matches, end <= 2000)$chr)


p <- ggplot(results_A_B,
       aes(x = log2FoldChange,
           y = -log10(padj),
           colour = has_rajika)) +
   geom_point() +
   scale_colour_viridis_d(name = "TALE target",
                          labels = c("Absent", "Present")) +
   labs(title = "SID4 dox+ vs. Mock")
print(p)


p <- ggplot(results_A_B,
       aes(x = baseMean,
           y = log2FoldChange,
           colour = has_rajika)) +
   geom_point() +
   scale_colour_viridis_d(name = "TALE target",
                          labels = c("Absent", "Present")) +
   labs(title = "SID4 dox+ vs. Mock")
print(p)


results_C_D <- results(dds_trn,
                       contrast = c("sample_ID", "C", "D")) %>%
                  as_tibble(rownames = c("subtelomere")) %>%
                  separate(subtelomere, into = "subtelomere", extra = "drop") %>%
                  mutate(subtelomere = factor(subtelomere,
                                              levels = subtelomere_levels,
                                              ordered = TRUE),
                         has_rajika = subtelomere %in% filter(rajika_seq_full_matches, end <= 2000)$chr)


p <- ggplot(results_C_D,
       aes(x = log2FoldChange,
           y = -log10(padj),
           colour = has_rajika)) +
   geom_point() +
   scale_colour_viridis_d(name = "TALE target",
                          labels = c("Absent", "Present")) +
   labs(title = "NLS3 dox+ vs. Mock")
print(p)

p <- ggplot(results_C_D,
       aes(x = baseMean,
           y = log2FoldChange,
           colour = has_rajika)) +
   geom_point() +
   scale_colour_viridis_d(name = "TALE target",
                          labels = c("Absent", "Present")) +
   labs(title = "NLS3 dox+ vs. Mock")
print(p)


results_both <- bind_rows(mutate(results_A_B, comparison = "SID4 dox+ vs. Mock"),
                          mutate(results_C_D, comparison = "NLS3 dox+ vs. Mock"))



p <- ggplot(results_both,
       aes(x = log2FoldChange,
           y = -log10(padj),
           colour = has_rajika)) +
   geom_point() +
   scale_colour_viridis_d(name = "TALE target",
                          labels = c("Absent", "Present")) +
   facet_wrap(~ comparison) +
   theme(aspect.ratio = 1)
print(p)

p <- ggplot(results_both,
       aes(x = baseMean,
           y = log2FoldChange,
           colour = has_rajika)) +
   geom_point() +
   scale_colour_viridis_d(name = "TALE target",
                          labels = c("Absent", "Present")) +
   facet_wrap(~ comparison) +
   theme(aspect.ratio = 1)
print(p)


p <- left_join(results_A_B,
          rajika_seq_full_matches %>% dplyr::count(chr, .drop = FALSE),
          by = c("subtelomere" = "chr")) %>%
   ggplot(aes(x = n, y = log2FoldChange)) +
   geom_point()
print(p)

p <- left_join(results_A_B,
          rajika_seq_full_matches %>% dplyr::count(chr, .drop = FALSE),
          by = c("subtelomere" = "chr")) %>%
   ggplot(aes(x = n, y = log2FoldChange, colour = log2(baseMean))) +
   geom_point() +
   scale_colour_viridis_c()
print(p)

p <- left_join(results_A_B,
          rajika_seq_full_matches %>% dplyr::count(chr, .drop = FALSE),
          by = c("subtelomere" = "chr")) %>%
   ggplot(aes(x = n, y = log2FoldChange, colour = -log10(padj))) +
   geom_point() +
   scale_colour_viridis_c()
print(p)


results_B_D <- results(dds_trn,
                       contrast = c("sample_ID", "B", "D")) %>%
                  as_tibble(rownames = c("subtelomere")) %>%
                  separate(subtelomere, into = "subtelomere", extra = "drop") %>%
                  mutate(subtelomere = factor(subtelomere,
                                              levels = subtelomere_levels,
                                              ordered = TRUE),
                         has_rajika = subtelomere %in% filter(rajika_seq_full_matches, end <= 2000)$chr)


p <- ggplot(results_C_D,
       aes(x = log2FoldChange,
           y = -log10(padj),
           colour = has_rajika)) +
   geom_point() +
   scale_colour_viridis_d(name = "TALE target",
                          labels = c("Absent", "Present")) +
   labs(title = "SID4 Mock vs. NLS3 Mock")
print(p)

p <- ggplot(results_C_D,
       aes(x = baseMean,
           y = log2FoldChange,
           colour = has_rajika)) +
   geom_point() +
   scale_colour_viridis_d(name = "TALE target",
                          labels = c("Absent", "Present")) +
   labs(title = "SID4 Mock vs. NLS3 Mock")
print(p)

results_all <- bind_rows(mutate(results_A_B, comparison = "SID4 dox+ vs. Mock"),
                         mutate(results_C_D, comparison = "NLS3 dox+ vs. Mock"),
                         mutate(results_B_D, comparison = "SID4 Mock vs. NLS3 Mock"))

p <- ggplot(results_all,
       aes(x = log2FoldChange,
           y = -log10(padj),
           colour = has_rajika)) +
   geom_point() +
   scale_colour_viridis_d(name = "TALE target",
                          labels = c("Absent", "Present")) +
   facet_wrap(~ comparison) +
   theme(aspect.ratio = 1)
print(p)

p <- ggplot(results_all,
       aes(x = baseMean,
           y = log2FoldChange,
           colour = has_rajika)) +
   geom_point() +
   geom_label_repel(data = filter(results_all, baseMean > 100),
                    mapping = aes(label = subtelomere)) +
   scale_colour_viridis_d(name = "TALE target",
                          labels = c("Absent", "Present")) +
   facet_wrap(~ comparison) +
   theme(aspect.ratio = 1)
print(p)

dev.off()
