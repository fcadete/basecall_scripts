
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(DESeq2)

pdf("DESeq2_dataset_comparisons.pdf", width = 16, height = 12)

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

# Load covreage frames for later plotting

load("coverage_frames/TALES_coverage_frame.RData")
load("coverage_frames/HEK_HeLa_coverage_frame.RData")
load("coverage_frames/other_coverage_frame.RData")

TALEs_total_reads <- bind_rows(read_tsv("guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded.porechop_finalcall_table"),
                         read_tsv("guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_2.porechop_finalcall_table"),
                         read_tsv("guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_2b.porechop_finalcall_table")) %>%
                 dplyr::count(final_barcode_call) %>%
                 select(barcode = final_barcode_call,
                        total_reads = n)

HEK_HeLa_total_reads <- read_tsv("guppy_VNP_purified_TERRA_HeLa-HEK293_20190807.porechop_finalcall_table") %>%
                 dplyr::count(final_barcode_call) %>%
                 select(barcode = final_barcode_call,
                        total_reads = n)

other_total_reads <- read_table("total_reads_per_file.txt",
                          col_names = c("total_reads", "file")) %>%
    filter(file != "total") %>% 
    mutate(file = str_remove(file, "guppy_"),
           file = str_remove(file, ".porechop_finalcall_table"))


# Load and do DESeq on TALEs data ----------------------------------------------

TALEs_barcode_to_condition <- read_tsv("TALEs/barcode_to_condition.tsv",
                                       col_names = c("barcode", "condition", "sample"))

TALEs_bam_files <- list.files(path = "TALEs/soft_clipped_ends", pattern = "*bam")

TALEs_bam_records <- lapply(as.list(paste0("TALEs/soft_clipped_ends/", TALEs_bam_files)),
                            readGAlignments)

TALEs_count_matrix <- do.call(cbind,
                              lapply(TALEs_bam_records,
                                     function(ranges) countOverlaps(subtelomere_ranges, ranges)))

colnames(TALEs_count_matrix) <- str_remove(TALEs_bam_files, "_on_rhietman_mapont_primary_subtelomere_start.bam")
rownames(TALEs_count_matrix) <- as.character(seqnames(subtelomere_ranges))


TALEs_count_frame <- TALEs_count_matrix %>%
                        as.data.frame() %>%
                        tibble::rownames_to_column(var = "subtelomere") %>%
                        gather(key = "file",
                               value = "counts",
                               -subtelomere)

TALEs_file_to_samples <- tibble(file = colnames(TALEs_count_matrix)) %>%
                            mutate(barcode = str_extract(file, "BC\\d\\d"),
                                   barcode = ifelse(is.na(barcode), "none", barcode)) %>%
                            left_join(TALEs_barcode_to_condition)

TALEs_sample_count_frame <- left_join(TALEs_count_frame, TALEs_file_to_samples) %>%
                               group_by(subtelomere, sample) %>%
                               summarise(counts = sum(counts)) %>%
                               spread(key = sample, value = counts)

TALEs_sample_count_matrix <- TALEs_sample_count_frame %>%
                                ungroup() %>%
                                select(2:5) %>%
                                as.matrix()

rownames(TALEs_sample_count_matrix) <- TALEs_sample_count_frame$subtelomere
colnames(TALEs_sample_count_matrix) <- str_replace(colnames(TALEs_sample_count_matrix),
                                                   "\\+", "_")

TALEs_coldata <- TALEs_file_to_samples %>%
                    select(sample, condition) %>%
                    unique() %>%
                    drop_na() %>%
                    as.data.frame() %>%
                    mutate(sample = str_replace(as.character(sample), "\\+", "_"),
                           condition = str_replace(as.character(condition), "\\+", "_"))
rownames(TALEs_coldata) <- TALEs_coldata$sample
TALEs_coldata <- TALEs_coldata[colnames(TALEs_sample_count_matrix), ]

dds_TALEs <- DESeqDataSetFromMatrix(countData = TALEs_sample_count_matrix,
                                    colData = TALEs_coldata,
                                    design = ~ condition)

dds_TALEs <- DESeq(dds_TALEs)

diff_subtelomeres_TALEs <- results(dds_TALEs) %>%
                           as_tibble(rownames = "subtelomere") %>%
                           filter(padj <= 0.05) %>%
                           separate(subtelomere, into = "subtelomere")

# Plot TALEs differentially-expressed subtelomeres------------------------------
left_join(filter(TALEs_coverage_frame, barcode != "none"),
          TALEs_file_to_samples) %>%
   group_by(chr, pos, condition, sample) %>%
   summarise(depth = sum(depth)) %>%
   ggplot(aes(x = pos,
              y = depth,
              group = sample,
              colour = condition)) +
      geom_line(alpha = I(1/2)) +
      facet_wrap(~ chr + (chr %in% diff_subtelomeres_TALEs$subtelomere))

TALEs_total_reads_per_sample <- left_join(TALEs_total_reads,
                                          select(TALEs_file_to_samples, barcode, condition, sample) %>%
                                             unique()) %>%
                                             group_by(condition, sample) %>%
                                             summarise(total_reads = sum(total_reads)) %>%
                                             drop_na()

left_join(filter(TALEs_coverage_frame, barcode != "none"),
          TALEs_file_to_samples) %>%
   group_by(chr, pos, condition, sample) %>%
   summarise(depth = sum(depth)) %>%
   left_join(TALEs_total_reads_per_sample) %>%
   ggplot(aes(x = pos,
              y = depth / total_reads,
              group = sample,
              colour = condition)) +
      geom_line(alpha = I(1/2)) +
      facet_wrap(~ chr + (chr %in% diff_subtelomeres_TALEs$subtelomere))

   

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

dds_HH <- DESeq(dds_HH)

diff_subtelomeres_HH <- results(dds_HH) %>%
                           as_tibble(rownames = "subtelomere") %>%
                           filter(padj <= 0.05) %>%
                           separate(subtelomere, into = "subtelomere")

# Plot HH differentially-expressed subtelomeres
left_join(filter(HEK_HeLa_coverage_frame, barcode != "none"),
          HH_file_to_samples) %>%
   ggplot(aes(x = pos,
              y = depth,
              group = barcode,
              colour = cell_type)) +
      geom_line(alpha = I(1/2)) +
      facet_wrap(~ chr + (chr %in% diff_subtelomeres_HH$subtelomere))

left_join(filter(HEK_HeLa_coverage_frame, barcode != "none"),
          HH_file_to_samples) %>%
   left_join(HEK_HeLa_total_reads) %>%
   ggplot(aes(x = pos,
              y = depth / total_reads,
              group = barcode,
              colour = cell_type)) +
      geom_line(alpha = I(1/2)) +
      facet_wrap(~ chr + (chr %in% diff_subtelomeres_HH$subtelomere))

# Load and do DESeq on original HeLa and U2OS samples --------------------------

HU_bam_files <- list.files(path = "soft_clipped_ends",
                           pattern = "*_on_rhietman_mapont_primary_subtelomere_start.bam")

HU_bam_files <- HU_bam_files[(str_detect(HU_bam_files, "GM|HEK|SAOS") == FALSE) &
                                str_detect(HU_bam_files, "Joana") &
                                str_detect(HU_bam_files, "purified") &
                                (str_detect(HU_bam_files, "pA") == FALSE)]

HU_bam_records <- lapply(as.list(paste0("soft_clipped_ends/", HU_bam_files)),
                         readGAlignments)

HU_count_matrix <- do.call(cbind,
                           lapply(HU_bam_records, function(ranges) countOverlaps(subtelomere_ranges,
                                                                                 ranges)))
colnames(HU_count_matrix) <- str_remove(HU_bam_files, "_on_rhietman_mapont_primary_subtelomere_start.bam")
rownames(HU_count_matrix) <- as.character(seqnames(subtelomere_ranges))

HU_file_to_samples <- tibble(file = colnames(HU_count_matrix),
                             sample = str_split_fixed(colnames(HU_count_matrix),
                                                      "_Joana",
                                                      n = 2)[, 1],
                             cell_type = ifelse(str_detect(colnames(HU_count_matrix),
                                                           "HeL"),
                                                "HeLa",
                                                "U2OS"))

HU_count_frame <- HU_count_matrix %>%
                     as.data.frame() %>%
                     tibble::rownames_to_column(var = "subtelomere") %>%
                     gather(key = "file",
                            value = "counts",
                            -subtelomere)

HU_sample_count_frame <- left_join(HU_count_frame, HU_file_to_samples) %>%
                            group_by(subtelomere, sample) %>%
                            summarise(counts = sum(counts)) %>%
                            spread(key = sample, value = counts)

HU_sample_count_matrix <- HU_sample_count_frame %>%
                             ungroup() %>%
                             select(2:6) %>%
                             as.matrix()
rownames(HU_sample_count_matrix) <- HU_sample_count_frame$subtelomere

HU_coldata <- HU_file_to_samples %>%
                 select(sample, cell_type) %>%
                 unique() %>%
                 drop_na() %>%
                 as.data.frame()
rownames(HU_coldata) <- HU_coldata$sample
HU_coldata <- HU_coldata[colnames(HU_sample_count_matrix), ]

dds_HU <- DESeqDataSetFromMatrix(countData = HU_sample_count_matrix,
                                 colData = HU_coldata,
                                 design = ~ cell_type)

dds_HU <- DESeq(dds_HU)

diff_subtelomeres_HU <- results(dds_HU) %>%
                           as_tibble(rownames = "subtelomere") %>%
                           filter(padj <= 0.05) %>%
                           separate(subtelomere, into = "subtelomere")

# Plot HU differentially-expressed subtelomeres
filter(other_coverage_frame,
       (str_detect(file, "GM|HEK|SAOS") == FALSE),
       str_detect(file, "purified"),
       (str_detect(file, "pA") == FALSE)) %>%
   ggplot(aes(x = pos,
              y = depth,
              group = file,
              colour = cell_type)) +
      geom_line(alpha = I(1/2)) +
      facet_wrap(~ chr + (chr %in% diff_subtelomeres_HU$subtelomere))

filter(other_coverage_frame,
       (str_detect(file, "GM|HEK|SAOS") == FALSE),
       str_detect(file, "purified"),
       (str_detect(file, "pA") == FALSE)) %>%
   left_join(other_total_reads) %>%
   ggplot(aes(x = pos,
              y = depth / total_reads,
              group = file,
              colour = cell_type)) +
      geom_line(alpha = I(1/2)) +
      facet_wrap(~ chr + (chr %in% diff_subtelomeres_HU$subtelomere))

# Compare the barcoded HEK to the U2OS samples

HEK_U2OS_sample_count_matrix <- cbind(HH_sample_count_matrix[, rownames(HH_coldata)[HH_coldata$cell_type == "HEK293T"]],
                                      HU_sample_count_matrix[, rownames(HU_coldata)[HU_coldata$cell_type == "U2OS"]])

HEK_U2OS_coldata <- rbind(HH_coldata[HH_coldata$cell_type == "HEK293T", c("sample", "cell_type")],
                          HU_coldata[HU_coldata$cell_type == "U2OS", c("sample", "cell_type")])

dds_HEK_U2OS <- DESeqDataSetFromMatrix(countData = HEK_U2OS_sample_count_matrix,
                                 colData = HEK_U2OS_coldata,
                                 design = ~ cell_type)

dds_HEK_U2OS <- DESeq(dds_HEK_U2OS)

diff_subtelomeres_HEK_U2OS <- results(dds_HEK_U2OS) %>%
                                 as_tibble(rownames = "subtelomere") %>%
                                 filter(padj <= 0.05) %>%
                                 separate(subtelomere, into = "subtelomere")

HEK_U2OS_frame_for_plotting <- rbind(
                    filter(other_coverage_frame,
                       (str_detect(file, "GM|HEK|SAOS|HeL") == FALSE),
                       str_detect(file, "purified"),
                       (str_detect(file, "pA") == FALSE)) %>%
                     select(chr, pos, cell_type, depth, total_reads, sample = file),
                    HEK_HeLa_coverage_frame %>%
                       filter(cell_type == "HEK293T") %>%
                       left_join(HEK_HeLa_total_reads) %>%
                       select(chr, pos, cell_type, depth, total_reads, sample = barcode))

ggplot(data = HEK_U2OS_frame_for_plotting,
       mapping = aes(x = pos, y = depth, colour = cell_type)) +
   geom_line(data = filter(other_coverage_frame,
                       (str_detect(file, "GM|HEK|SAOS|HeL") == FALSE),
                       str_detect(file, "purified"),
                       (str_detect(file, "pA") == FALSE)),
             mapping = aes(group = file)) +
   geom_line(data = HEK_HeLa_coverage_frame %>%
                       filter(cell_type == "HEK293T"),
             mapping = aes(group = barcode)) +
   facet_wrap(~ chr + (chr %in% diff_subtelomeres_HEK_U2OS$subtelomere))

ggplot(mapping = aes(x = pos, y = depth / total_reads, colour = cell_type)) +
   geom_line(data = filter(other_coverage_frame,
                       (str_detect(file, "GM|HEK|SAOS|HeL") == FALSE),
                       str_detect(file, "purified"),
                       (str_detect(file, "pA") == FALSE)),
             mapping = aes(group = file)) +
   geom_line(data = HEK_HeLa_coverage_frame %>%
                       filter(cell_type == "HEK293T") %>%
                       left_join(HEK_HeLa_total_reads),
             mapping = aes(group = barcode)) +
   facet_wrap(~ chr )



dev.off()



subtelomere_small_levels <- c("1p", "1q", "2p", "2q", "3p",
                              "3q", "4p", "4q", "5p", "5q",
                              "6p", "6q", "7q", "7p", "8p", "8q",
                              "9p", "9q", "10p", "10q", "11p", "11q",
                              "12p", "12q", "13q", "14q", "15q", "16p",
                              "16q", "17p", "17q", "18p", "18q", "19p",
                              "19q", "20p", "20q", "21q", "22q",
                              "XpYp", "Xq", "Yq")


HEK_HeLa_total_reads <- left_join(HH_file_to_samples,
                                  HEK_HeLa_total_reads)

HEK_HeLa_total_reads <- as.data.frame(HEK_HeLa_total_reads,
                                      stringsAsFactors = FALSE)
rownames(HEK_HeLa_total_reads) <- HEK_HeLa_total_reads$sample

pdf("HH_count_comparison.pdf", width = 16, height = 12)

HH_sample_count_frame %>%
   gather(-subtelomere,
          key = sample,
          value = counts) %>%
   left_join(HH_file_to_samples) %>%
   separate(subtelomere, into = c("subtelomere"), sep = "tel_") %>%
   mutate(subtelomere = factor(subtelomere,
                               levels = subtelomere_small_levels,
                               ordered = TRUE),
          diff_expressed = subtelomere %in% str_remove(diff_subtelomeres_HH$subtelomere, "tel")) %>%
   drop_na() %>%
   ggplot(aes(x = cell_type,
              y = counts + 0.1,
              colour = cell_type)) + 
      geom_jitter(height = 0) +
      facet_wrap(~ subtelomere + diff_expressed) +
      scale_y_log10()

HH_sample_count_frame %>%
   gather(-subtelomere,
          key = sample,
          value = counts) %>%
   left_join(HEK_HeLa_total_reads) %>%
   separate(subtelomere, into = c("subtelomere"), sep = "tel_") %>%
   mutate(subtelomere = factor(subtelomere,
                               levels = subtelomere_small_levels,
                               ordered = TRUE),
          diff_expressed = subtelomere %in% str_remove(diff_subtelomeres_HH$subtelomere, "tel")) %>%
   drop_na() %>%
   ggplot(aes(x = cell_type,
              y = (counts + 0.1) / total_reads,
              colour = cell_type)) + 
      geom_jitter(height = 0) +
      facet_wrap(~ subtelomere + diff_expressed) +
      scale_y_log10()

dev.off()

diff_subtelomeres_loo <- list()


for (sample_to_remove in colnames(HH_sample_count_matrix)) {

   print(sample_to_remove)
   
   HH_sample_count_matrix_loo <- HH_sample_count_matrix[, (colnames(HH_sample_count_matrix) != sample_to_remove)]

   HH_coldata_loo <- HH_coldata[rownames(HH_coldata) != sample_to_remove, ]

   print(HH_coldata_loo)

   dds_HH_loo <- DESeqDataSetFromMatrix(countData = HH_sample_count_matrix_loo,
                                        colData = HH_coldata_loo,
                                        design = ~ cell_type)
   
   dds_HH_loo <- DESeq(dds_HH_loo)
   
   diff_subtelomeres_loo <- append(diff_subtelomeres_loo,
                                   list((results(dds_HH_loo) %>%
                                       as_tibble(rownames = "subtelomere") %>%
                                       filter(padj <= 0.05) %>%
                                       separate(subtelomere, into = "subtelomere"))$subtelomere))
   
}

names(diff_subtelomeres_loo) <- colnames(HH_sample_count_matrix)

pdf("TALEs_count_comparison.pdf", width = 16, height = 12)

TALEs_sample_count_frame %>%
   gather(-subtelomere,
          key = sample,
          value = counts) %>%
   mutate(sample = str_replace(sample, "\\+", "_")) %>%
   left_join(TALEs_coldata) %>%
   separate(subtelomere,
            into = c("subtelomere"),
            sep = "tel_") %>%
   mutate(subtelomere = factor(subtelomere,
                               levels = subtelomere_small_levels,
                               ordered = TRUE)) %>%
   drop_na() %>%
   ggplot(aes(x = condition,
              y = counts + 0.1,
              colour = condition)) + 
      geom_jitter(height = 0) +
      facet_wrap(~ subtelomere) +
      scale_y_log10()


TALEs_total_reads_per_sample <- ungroup(TALEs_total_reads_per_sample) %>%
                                   mutate(condition = str_replace(condition, "\\+", "_"),
                                          sample = str_replace(sample, "\\+", "_"))

TALEs_sample_count_frame %>%
   gather(-subtelomere,
          key = sample,
          value = counts) %>%
   mutate(sample = str_replace(sample, "\\+", "_")) %>%
   left_join(TALEs_total_reads_per_sample) %>%
   separate(subtelomere, into = c("subtelomere"), sep = "tel_") %>%
   mutate(subtelomere = factor(subtelomere,
                               levels = subtelomere_small_levels,
                               ordered = TRUE)) %>%
   drop_na() %>%
   ggplot(aes(x = condition,
              y = (counts + 0.1) / total_reads,
              colour = condition)) + 
      geom_jitter(height = 0) +
      facet_wrap(~ subtelomere) +
      scale_y_log10()

dev.off()

