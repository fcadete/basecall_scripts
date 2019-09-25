
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

# Read in total read frames

HEK_HeLa_total_reads <- read_tsv("guppy_VNP_purified_TERRA_HeLa-HEK293_20190807.porechop_finalcall_table") %>%
                 dplyr::count(final_barcode_call) %>%
                 select(barcode = final_barcode_call,
                        total_reads = n)

other_total_reads <- read_table("total_reads_per_file.txt",
                          col_names = c("total_reads", "file")) %>%
    filter(file != "total") %>% 
    mutate(file = str_remove(file, "guppy_"),
           file = str_remove(file, ".porechop_finalcall_table"))


# Load and do DESEq on barcoded HeLa and HEK data ------------------------------

HH_barcode_to_condition <- read_tsv("190807_Seventh_run/Barcodes_clean.txt",
                                    col_names = c("cell_type", "date", "barcode")) %>%
                              mutate(barcode = str_replace(barcode, "NB", "BC"),
                                     file = paste(cell_type, date, sep = "_"))

HH_barcode_to_condition <- left_join(HH_barcode_to_condition,
                                     HEK_HeLa_total_reads)

HH_bam_files <- paste0("HeLa-HEK_barcoded/terra_primer_separated/",
                       HH_barcode_to_condition$barcode,
                       "/Joana_TERRA_VNP_no_PCR_on_rhietman_mapont_primary_subtelomere_start.bam")

HH_bam_records <- lapply(as.list(HH_bam_files),
                         readGAlignments)

HH_count_matrix <- do.call(cbind,
                           lapply(HH_bam_records, function(ranges) countOverlaps(subtelomere_ranges,
                                                                                 ranges)))

colnames(HH_count_matrix) <- paste0(HH_barcode_to_condition$cell_type,
                                    "_",
                                    HH_barcode_to_condition$date)

rownames(HH_count_matrix) <- as.character(seqnames(subtelomere_ranges))

HH_count_frame <- HH_count_matrix %>%
                     as.data.frame() %>%
                     tibble::rownames_to_column(var = "subtelomere") %>%
                     gather(key = "file",
                            value = "counts",
                            -subtelomere)

HH_count_frame <- left_join(HH_count_frame, HH_barcode_to_condition) 

HH_count_frame <- select(HH_count_frame,
                         cell_type, subtelomere, sample = file, counts)

HH_count_frame <- left_join(HH_count_frame,
                            select(HH_barcode_to_condition,
                                   sample = file, total_reads))


# Load and do DESeq on original HeLa and U2OS samples --------------------------

HU_bam_files <- list.files(path = "soft_clipped_ends",
                           pattern = "*_on_rhietman_mapont_primary_subtelomere_start.bam")

HU_bam_files <- HU_bam_files[str_detect(HU_bam_files, "Joana_TERRA") &
                             str_detect(HU_bam_files, "purified") &
                             (str_detect(HU_bam_files, "pA") == FALSE)]


HU_bam_records <- lapply(as.list(paste0("soft_clipped_ends/", HU_bam_files)),
                         readGAlignments)

HU_count_matrix <- do.call(cbind,
                           lapply(HU_bam_records, function(ranges) countOverlaps(subtelomere_ranges,
                                                                                 ranges)))
colnames(HU_count_matrix) <- str_remove(HU_bam_files, "_on_rhietman_mapont_primary_subtelomere_start.bam")
rownames(HU_count_matrix) <- as.character(seqnames(subtelomere_ranges))

HU_count_frame <- HU_count_matrix %>%
                     as.data.frame() %>%
                     tibble::rownames_to_column(var = "subtelomere") %>%
                     gather(key = "file",
                            value = "counts",
                            -subtelomere)


HU_count_frame <- HU_count_frame %>%
                     mutate(sample = str_split_fixed(file,
                                                      "_Joana",
                                                      n = 2)[, 1],
                            cell_type = ifelse(grepl("HeL", file),
                                               "HeLa",
                                               ifelse(grepl("GM847", file),
                                                      "GM847",
                                                      ifelse(grepl("SAOS2", file),
                                                             "SAOS2",
                                                             ifelse(grepl("HEK", file),
                                                                    "HEK293T",
                                                                    "U2OS"))))) 


HU_count_frame <- HU_count_frame %>%
                     group_by(cell_type, subtelomere, sample) %>%
                     summarise(counts = sum(counts)) %>%
                     ungroup()


HU_count_frame <- left_join(HU_count_frame,
                            select(other_total_reads,
                                   sample = file, total_reads))


full_count_frame <- bind_rows(HH_count_frame, HU_count_frame) %>%
                       separate(subtelomere, into = c("chr"))


tel_29_bp_frame <- read_delim("tel_29bp_blast_results.txt",
                              delim = "\t",
                              col_names=c("qseqid", "chr", "pident",
                                          "length", "mismatch", "gapopen",
                                           "qstart", "qend",
                                           "start", "end",
                                           "evalue", "bitscore"))


tel_29_bp_frame <- tel_29_bp_frame %>% separate(chr, into = c("chr"), extra = "drop")

tel_29_bp_frame$chr <- factor(tel_29_bp_frame$chr,
                          levels = subtelomere_levels,
                          ordered = TRUE)


tel_29_bp_summary <- tel_29_bp_frame %>%
                        filter(start <= 2000) %>%
                        group_by(chr) %>%
                        summarise(tel29bp_repeats      = n(),
                                  tel29bp_full_matches = sum(bitscore == 58),
                                  tel29bp_bitscore_sum = sum(bitscore),
                                  tel29bp_bitscore_avg = mean(bitscore))

tel_29_bp_summary <- left_join(tibble(chr = subtelomere_levels),
                               tel_29_bp_summary) %>%
                        tidyr::replace_na(list(tel29bp_repeats      = 0,
                                               tel29bp_full_matches = 0,
                                               tel29bp_bitscore_sum = 0,
                                               tel29bp_bitscore_avg = 0))


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


rajika_summary <- rajika_seq_full_matches %>%
                     filter(start <= 2000) %>%
                     group_by(chr) %>%
                     summarise(rajika_repeats      = n(),
                               rajika_bitscore_sum = sum(bitscore),
                               rajika_bitscore_avg = mean(bitscore))

rajika_summary <- left_join(tibble(chr = subtelomere_levels),
                            rajika_summary) %>%
                     tidyr::replace_na(list(rajika_repeats      = 0,
                                            rajika_bitscore_sum = 0,
                                            rajika_bitscore_avg = 0))



count_repeat_frame <- left_join(full_count_frame,
                                tel_29_bp_summary)


count_repeat_frame <- left_join(count_repeat_frame,
                                rajika_summary)


count_matrix <- count_repeat_frame %>%
                   select(chr, sample, counts) %>%
                   spread(key = sample, value = counts)
rownames(count_matrix) <- count_matrix$chr
count_matrix <- count_matrix[, -1]

count_matrix <- count_matrix[, colSums(count_matrix) != 0]

coldata <- unique(select(count_repeat_frame, sample, cell_type))
rownames(coldata) <- coldata$sample

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata[colnames(count_matrix),],
                              design = ~ cell_type)

rlogs <- rlog(dds, blind = FALSE)

rlog_frame <- as_tibble(assay(rlogs))
rlog_frame$chr <-rownames(count_matrix)

rlog_frame <- gather(rlog_frame, -chr,
                     key = "sample", value = "rlog")


count_repeat_frame <- left_join(count_repeat_frame, rlog_frame)

total_reads_frame <- select(count_repeat_frame,
                            sample, total_reads) %>%
                        unique()

total_reads <- total_reads_frame$total_reads
names(total_reads) <- total_reads_frame$sample

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata[colnames(count_matrix),],
                              design = ~ cell_type)

sizeFactors(dds) <- total_reads[colnames(count_matrix)] / mean(total_reads[colnames(count_matrix)])

rlogs_trn <- rlog(dds, blind = FALSE)

rlog_frame <- as_tibble(assay(rlogs_trn))
rlog_frame$chr <-rownames(count_matrix)

rlog_frame <- gather(rlog_frame, -chr,
                     key = "sample", value = "rlog_trn")


count_repeat_frame <- left_join(count_repeat_frame, rlog_frame)

write_tsv(count_repeat_frame,
          path = "count_repeat_frame.tsv")



