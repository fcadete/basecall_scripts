
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


# Load barcoded HeLa and HEK data ------------------------------

HH_barcode_to_condition <- read_tsv("190807_Seventh_run/Barcodes_clean.txt",
                                    col_names = c("cell_type", "date", "barcode")) %>%
                              mutate(barcode = str_replace(barcode, "NB", "BC"),
                                     file = paste(cell_type, date, sep = "_"))

HH_barcode_to_condition <- left_join(HH_barcode_to_condition,
                                     HEK_HeLa_total_reads)

HH_bam_files <- paste0("reads_mapping_to_very_end/HeLa-HEK_barcoded/",
                       HH_barcode_to_condition$barcode,
                       "/Joana_TERRA_VNP_no_PCR_on_rhietman_mapont_primary_very_end.bam")

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


# Load original HeLa and U2OS samples --------------------------

HU_bam_files <- list.files(path = "reads_mapping_to_very_end",
                           pattern = "*_primary_very_end.bam")

HU_bam_files <- HU_bam_files[str_detect(HU_bam_files, "Joana_TERRA") &
                             str_detect(HU_bam_files, "purified") &
                             (str_detect(HU_bam_files, "pA") == FALSE)]


HU_bam_records <- lapply(as.list(paste0("reads_mapping_to_very_end/", HU_bam_files)),
                         readGAlignments)

HU_count_matrix <- do.call(cbind,
                           lapply(HU_bam_records, function(ranges) countOverlaps(subtelomere_ranges,
                                                                                 ranges)))

colnames(HU_count_matrix) <- str_remove(HU_bam_files, "on_rhietman_mapont_primary_very_end.bam")
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


count_matrix <- full_count_frame %>%
                   select(chr, sample, counts) %>%
                   spread(key = sample, value = counts)
rownames(count_matrix) <- count_matrix$chr
count_matrix <- count_matrix[, -1]

count_matrix <- count_matrix[, colSums(count_matrix) != 0]

coldata <- unique(select(full_count_frame, sample, cell_type))
rownames(coldata) <- coldata$sample


total_reads_frame <- select(full_count_frame,
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


count_repeat_frame <- left_join(full_count_frame, rlog_frame)

write_tsv(count_repeat_frame,
          path = "count_repeat_frame_reads_mapping_very_end_primary.tsv")



