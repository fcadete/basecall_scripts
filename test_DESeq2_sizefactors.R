
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(DESeq2)


# Set common subtelomere levels and ranges to be used --------------------------

subtelomere_levels <- c("1ptel", "1qtel", "2ptel", "2qtel", "3ptel",
                        "3qtel", "4ptel", "4qtel", "5ptel", "5qtel",
                        "6ptel", "6qtel", "7qtel", "7ptel", "8ptel", "8qtel",
                        "9ptel", "9qtel", "10ptel", "10qtel", "11ptel", "11qtel",
                        "12ptel", "12qtel", "13qtel", "14qtel", "15qtel", "16ptel",
                        "16qtel", "17ptel", "17qtel", "18ptel", "18qtel", "19ptel",
                        "19qtel", "20ptel", "20qtel", "21qtel", "22qtel",
                        "XpYptel", "Xqtel", "Yqtel")

subtelomere_bed <- read_tsv("subtelomere_starts_1_5000.bed", col_names = c("seqnames", "start", "end"))

subtelomere_bed <- subtelomere_bed %>% mutate(end = 2000)

subtelomere_ranges <- GRanges(subtelomere_bed)


# Load coverage and total reads frames------------------------------------------

load("coverage_frames/HEK_HeLa_coverage_frame.RData")

HEK_HeLa_total_reads <- read_tsv("guppy_VNP_purified_TERRA_HeLa-HEK293_20190807.porechop_finalcall_table") %>%
                 dplyr::count(final_barcode_call) %>%
                 select(barcode = final_barcode_call,
                        total_reads = n)


# Load and do DESEq on barcoded HeLa and HEK data ------------------------------

HH_barcode_to_condition <- read_tsv("190807_Seventh_run/Barcodes_clean.txt",
                                    col_names = c("cell_type", "date", "barcode")) %>%
                              mutate(barcode = str_replace(barcode, "NB", "BC"))

HH_bam_files <- list.files(path = "HeLa-HEK_barcoded/soft_clipped_ends", pattern = "*bam")

HH_bam_records <- lapply(as.list(paste0("HeLa-HEK_barcoded/soft_clipped_ends/", HH_bam_files)),
                         readGAlignments)

HH_count_matrix <- do.call(cbind,
                           lapply(HH_bam_records, function(ranges) countOverlaps(subtelomere_ranges,
                                                                                 ranges)))
colnames(HH_count_matrix) <- str_remove(HH_bam_files, "_on_rhietman_mapont_primary_subtelomere_start.bam")
rownames(HH_count_matrix) <- as.character(seqnames(subtelomere_ranges))

HH_count_frame <- HH_count_matrix %>%
                     as.data.frame() %>%
                     tibble::rownames_to_column(var = "subtelomere") %>%
                     gather(key = "file",
                            value = "counts",
                            -subtelomere)

HH_file_to_samples <- tibble(file = colnames(HH_count_matrix)) %>%
                         mutate(barcode = str_extract(file, "BC\\d\\d"),
                                barcode = ifelse(is.na(barcode), "none", barcode)) %>%
                         left_join(HH_barcode_to_condition) %>%
                         mutate(sample = paste(cell_type, date, sep = "_"))

HH_sample_count_frame <- left_join(HH_count_frame, HH_file_to_samples) %>%
                            group_by(subtelomere, sample) %>%
                            summarise(counts = sum(counts)) %>%
                            spread(key = sample, value = counts)

HH_sample_count_matrix <- HH_sample_count_frame %>%
                             ungroup() %>%
                             select(2:9) %>%
                             as.matrix()
rownames(HH_sample_count_matrix) <- HH_sample_count_frame$subtelomere

HH_coldata <- HH_file_to_samples %>%
                 select(sample, cell_type, date) %>%
                 unique() %>%
                 drop_na() %>%
                 as.data.frame()
rownames(HH_coldata) <- HH_coldata$sample
HH_coldata <- HH_coldata[colnames(HH_sample_count_matrix), ]

dds_HH <- DESeqDataSetFromMatrix(countData = HH_sample_count_matrix,
                                 colData = HH_coldata,
                                 design = ~ cell_type)

dds_HH_deseq_size_factors <- DESeq(dds_HH)

diff_subtelomeres_HH_deseq_size_factors <- results(dds_HH_deseq_size_factors) %>%
                                              as_tibble(rownames = "subtelomere") %>%
                                              filter(padj <= 0.05) %>%
                                              separate(subtelomere, into = "subtelomere")


HH_total_reads <- filter(left_join(HH_file_to_samples, HEK_HeLa_total_reads), barcode != "none")

HH_sizes_vector <- HH_total_reads$total_reads
names(HH_sizes_vector) <- HH_total_reads$sample
HH_sizes_vector <- HH_sizes_vector[names(sizeFactors(dds_HH_deseq_size_factors))]

dds_HH_total_reads_size_factors <- DESeqDataSetFromMatrix(countData = HH_sample_count_matrix,
                                                          colData = HH_coldata,
                                                          design = ~ cell_type)

sizeFactors(dds_HH_total_reads_size_factors) <- HH_sizes_vector / median(HH_sizes_vector)

dds_HH_total_reads_size_factors <- DESeq(dds_HH_total_reads_size_factors)

diff_subtelomeres_HH_total_reads_size_factors <- results(dds_HH_total_reads_size_factors) %>%
                                                    as_tibble(rownames = "subtelomere") %>%
                                                    filter(padj <= 0.05) %>%
                                                    separate(subtelomere, into = "subtelomere")



dds_HH_iterated_size_factors <- DESeqDataSetFromMatrix(countData = HH_sample_count_matrix,
                                                       colData = HH_coldata,
                                                       design = ~ cell_type)

dds_HH_iterated_size_factors <- estimateSizeFactors(dds_HH_iterated_size_factors,
                                                    type = "iterate")

dds_HH_iterated_size_factors <- DESeq(dds_HH_iterated_size_factors)

diff_subtelomeres_HH_iterated_size_factors <- results(dds_HH_iterated_size_factors) %>%
                                                    as_tibble(rownames = "subtelomere") %>%
                                                    filter(padj <= 0.05) %>%
                                                    separate(subtelomere, into = "subtelomere")



dds_HH_poscounts_size_factors <- DESeqDataSetFromMatrix(countData = HH_sample_count_matrix,
                                                        colData = HH_coldata,
                                                        design = ~ cell_type)

dds_HH_poscounts_size_factors <- estimateSizeFactors(dds_HH_poscounts_size_factors,
                                                     type = "poscounts")

dds_HH_poscounts_size_factors <- DESeq(dds_HH_poscounts_size_factors)

diff_subtelomeres_HH_poscounts_size_factors <- results(dds_HH_poscounts_size_factors) %>%
                                                  as_tibble(rownames = "subtelomere") %>%
                                                  filter(padj <= 0.05) %>%
                                                  separate(subtelomere, into = "subtelomere")



