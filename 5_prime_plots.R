
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(Rsamtools)
library(GenomicAlignments)
library(seqLogo)

alignment_mapont_files <- list.files(path = "alignment_outputs/",
                                    pattern = "*rhietman_mapont_primary.bam")

alignment_mapont_files <- alignment_mapont_files[grepl("bam.bai", alignment_mapont_files) == FALSE]

mapping_frame <- data_frame()

for (file in alignment_mapont_files) {
  
  if (grepl("Joana_TERRA_VNP", file)) {

     alignments <- readGAlignments(paste0("alignment_outputs/", file)) %>% as.data.frame()

     if (nrow(alignments) > 0) {
        mapping_frame <- rbind(mapping_frame,
                               cbind(data_frame(file = file),
                                     alignments))
     }
  }
}

mapping_frame <- as_data_frame(mapping_frame)

mapping_frame <- mapping_frame %>%
                 mutate(long_seqnames = seqnames) %>%
                 separate(seqnames, into = c("seqnames"), extra = "drop") %>%
                 separate(file, sep = "_Joana_TERRA_VNP", into = c("file"), extra = "drop")

long_subtelomere_names <- mapping_frame %>%
                            select(seqnames, long_seqnames) %>%
                            distinct()

mapping_frame$seqnames <- factor(mapping_frame$seqnames,
                          levels = c("1ptel", "1qtel", "2ptel", "2qtel", "3ptel",
                                     "3qtel", "4ptel", "4qtel", "5ptel", "5qtel",
                                     "6ptel", "6qtel", "7qtel", "7ptel", "8ptel", "8qtel",
                                     "9ptel", "9qtel", "10ptel", "10qtel", "11ptel", "11qtel",
                                     "12ptel", "12qtel", "13qtel", "14qtel", "15qtel", "16ptel",
                                     "16qtel", "17ptel", "17qtel", "18ptel", "18qtel", "19ptel",
                                     "19qtel", "20ptel", "20qtel", "21qtel", "22qtel",
                                     "XpYptel", "Xqtel", "Yqtel"),
                          ordered = TRUE)

mapping_ranges <- GRanges(mapping_frame)

subtelomere_end_ranges <- with(long_subtelomere_names,
                               GRanges(seqnames = seqnames,
                                       ranges = IRanges(start = 1,
                                                        end = 2000),
                                       long_seqnames = long_seqnames))

subtelomere_end_ranges$mapped_reads <- countOverlaps(subtelomere_end_ranges, mapping_ranges)

as_tibble(subtelomere_end_ranges) %>%
   write_tsv(path = "subtelomere_aligned_counts_all.tsv")


tel_29_bp_frame <- read_delim("tel_29bp_blast_results.txt",
                              delim = "\t",
                              col_names=c("qseqid", "seqnames", "pident",
                                          "length", "mismatch", "gapopen",
                                           "qstart", "qend",
                                           "start", "end",
                                           "evalue", "bitscore"))

tel_29_bp_frame <- tel_29_bp_frame %>% separate(seqnames, into = c("seqnames"), extra = "drop")

tel_29_bp_frame$seqnames <- factor(tel_29_bp_frame$seqnames,
                          levels = c("1ptel", "1qtel", "2ptel", "2qtel", "3ptel",
                                     "3qtel", "4ptel", "4qtel", "5ptel", "5qtel",
                                     "6ptel", "6qtel", "7qtel", "7ptel", "8ptel", "8qtel",
                                     "9ptel", "9qtel", "10ptel", "10qtel", "11ptel", "11qtel",
                                     "12ptel", "12qtel", "13qtel", "14qtel", "15qtel", "16ptel",
                                     "16qtel", "17ptel", "17qtel", "18ptel", "18qtel", "19ptel",
                                     "19qtel", "20ptel", "20qtel", "21qtel", "22qtel",
                                     "XpYptel", "Xqtel", "Yqtel"),
                          ordered = TRUE)


for (this_file in unique(mapping_frame$file)) {
  ggsave(paste0("coverage_plots/alignment_ends_", this_file, ".png"),
         rbind(mapping_frame %>% filter(file == this_file) %>%
                 group_by(file, seqnames) %>%
                 dplyr::rename(pos = start) %>%
                 dplyr::count(pos) %>%
                 mutate(end = "telomere_proximal"),
               mapping_frame %>% filter(file == this_file) %>%
                 group_by(file, seqnames) %>%
                 dplyr::rename(pos = end) %>%
                 dplyr::count(pos) %>%
                 mutate(end = "telomere_distal")) %>%
                 filter(pos <= 50000) %>%
           ggplot(aes(x = pos, y = n, colour = end)) +
              geom_point() +
              theme_bw() +
              facet_wrap(~ seqnames),
      width = 20, height = 15, units = "cm", device = png())
   dev.off()
}



mapping_frame <- rbind(mapping_frame %>%
                        mutate(cell_type = ifelse(grepl("HeL", file), "HeLa", "U2OS")) %>%
                        group_by(cell_type, seqnames) %>%
                        dplyr::rename(pos = start) %>%
                        dplyr::count(pos) %>%
                        mutate(end = "telomere_proximal"),
                      mapping_frame %>%
                        mutate(cell_type = ifelse(grepl("HeL", file), "HeLa", "U2OS")) %>%
                        group_by(cell_type, seqnames) %>%
                        dplyr::rename(pos = end) %>%
                        dplyr::count(pos) %>%
                        mutate(end = "telomere_distal"))

mapping_frame$seqnames <- factor(mapping_frame$seqnames,
                          levels = c("1ptel", "1qtel", "2ptel", "2qtel", "3ptel",
                                     "3qtel", "4ptel", "4qtel", "5ptel", "5qtel",
                                     "6ptel", "6qtel", "7qtel", "7ptel", "8ptel", "8qtel",
                                     "9ptel", "9qtel", "10ptel", "10qtel", "11ptel", "11qtel",
                                     "12ptel", "12qtel", "13qtel", "14qtel", "15qtel", "16ptel",
                                     "16qtel", "17ptel", "17qtel", "18ptel", "18qtel", "19ptel",
                                     "19qtel", "20ptel", "20qtel", "21qtel", "22qtel",
                                     "XpYptel", "Xqtel", "Yqtel"),
                          ordered = TRUE)



png("coverage_plots/alignment_ends_pooled_U2OS.png", width = 1280, height = 1024)

mapping_frame %>%
  filter(cell_type == "U2OS") %>%
  ggplot() +
  geom_rect(data = tel_29_bp_frame,
            mapping = aes(xmin = start,
                          xmax = end,
                          ymin = 0,
                          ymax = Inf,
                          fill = bitscore),
            alpha = I(1/2)) +
  scale_fill_viridis_c() +
  geom_point(mapping = aes(x = pos, y = n, colour = end)) +
  facet_wrap( ~ seqnames) +
  scale_x_continuous(name = "Distance from telomere repeats (bp)",
                     breaks = c(1, 100000, 200000, 300000, 400000, 500000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

png("coverage_plots/alignment_ends_pooled_HeLa.png", width = 1280, height = 1024)

mapping_frame %>%
  filter(cell_type == "HeLa") %>%
  ggplot() +
  geom_rect(data = tel_29_bp_frame,
            mapping = aes(xmin = start,
                          xmax = end,
                          ymin = 0,
                          ymax = Inf,
                          fill = bitscore),
            alpha = I(1/2)) +
  scale_fill_viridis_c() +
  geom_point(mapping = aes(x = pos, y = n, colour = end)) +
  facet_wrap( ~ seqnames) +
  scale_x_continuous(name = "Distance from telomere repeats (bp)",
                     breaks = c(1, 100000, 200000, 300000, 400000, 500000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



dev.off()

png("coverage_plots/alignment_ends_pooled_U2OS_zoomed.png", width = 1280, height = 1024)

mapping_frame %>%
  filter(cell_type == "U2OS") %>%
  filter(pos <= 10000) %>%
  droplevels() %>%
  ggplot() +
  geom_rect(data = tel_29_bp_frame,
            mapping = aes(xmin = start,
                          xmax = end,
                          ymin = 0,
                          ymax = Inf,
                          fill = bitscore),
            alpha = I(1/2)) +
  scale_fill_viridis_c() +
  geom_point(mapping = aes(x = pos, y = n, colour = end)) +
  facet_wrap( ~ seqnames, scales = "free_y") +
  scale_x_continuous(name = "Distance from telomere repeats (bp)",
                     breaks = c(1, 2000, 4000, 6000, 8000, 10000),
                     limits = c(1, 10000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

png("coverage_plots/alignment_ends_pooled_HeLa_zoomed.png", width = 1280, height = 1024)

mapping_frame %>%
  filter(cell_type == "HeLa") %>%
  filter(pos <= 10000) %>%
  droplevels() %>%
  ggplot() +
  geom_rect(data = tel_29_bp_frame,
            mapping = aes(xmin = start,
                          xmax = end,
                          ymin = 0,
                          ymax = Inf,
                          fill = bitscore),
            alpha = I(1/2)) +
  scale_fill_viridis_c() +
  geom_point(mapping = aes(x = pos, y = n, colour = end)) +
  facet_wrap( ~ seqnames, scales = "free_y") +
  scale_x_continuous(name = "Distance from telomere repeats (bp)",
                     breaks = c(1, 2000, 4000, 6000, 8000, 10000),
                     limits = c(1, 10000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


png("coverage_plots/alignment_ends_pooled_U2OS_extra_zoomed.png", width = 1280, height = 1024)

mapping_frame %>%
  filter(cell_type == "U2OS") %>%
  filter(pos <= 2000) %>%
  droplevels() %>%
  ggplot() +
  geom_rect(data = tel_29_bp_frame,
            mapping = aes(xmin = start,
                          xmax = end,
                          ymin = 0,
                          ymax = Inf,
                          fill = bitscore),
            alpha = I(1/2)) +
  scale_fill_viridis_c() +
  geom_point(mapping = aes(x = pos, y = n, colour = end)) +
  facet_wrap( ~ seqnames, scales = "free_y") +
  scale_x_continuous(name = "Distance from telomere repeats (bp)",
                     breaks = c(1, 500, 1000, 1500, 2000),
                     limits = c(1, 2000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

png("coverage_plots/alignment_ends_pooled_HeLa_extra_zoomed.png", width = 1280, height = 1024)

mapping_frame %>%
  filter(cell_type == "HeLa") %>%
  filter(pos <= 2000) %>%
  droplevels() %>%
  ggplot() +
  geom_rect(data = tel_29_bp_frame,
            mapping = aes(xmin = start,
                          xmax = end,
                          ymin = 0,
                          ymax = Inf,
                          fill = bitscore),
            alpha = I(1/2)) +
  scale_fill_viridis_c() +
  geom_point(mapping = aes(x = pos, y = n, colour = end)) +
  facet_wrap( ~ seqnames, scales = "free_y") +
  scale_x_continuous(name = "Distance from telomere repeats (bp)",
                     breaks = c(1, 500, 1000, 1500, 2000),
                     limits = c(1, 2000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


mapping_frame %>%
   filter(cell_type == "U2OS") %>%
   filter(pos <= 2000) %>%
   group_by(end, seqnames) %>%
   summarise(max_coord = pos[which.max(n)],
             start_coord = max_coord - 100,
             end_coord = max_coord + 100) %>%
   ungroup() %>%
   filter(end == "telomere_distal") %>%
   left_join(long_subtelomere_names) %>%
   select(long_seqnames, start_coord, end_coord) %>%
   write_tsv(path = "telomere_distal_alignment_maxima.tsv")


five_prime_max_coords <- mapping_frame %>%
   filter(cell_type == "U2OS") %>%
   filter(pos <= 2000) %>%
   group_by(end, seqnames) %>%
   summarise(max_coord = pos[which.max(n)],
             start_coord = max_coord - 100,
             end_coord = max_coord + 100) %>%
   ungroup() %>%
   filter(end == "telomere_distal") %>%
   left_join(long_subtelomere_names) %>%
   select(long_seqnames, start_coord, end_coord)

subtelomere_sequences <- readDNAStringSet("references/ConcatenatedFASTAAassemblies_hTel.txt",
                                          format = "FASTA")

TSS_sequences <- subseq(subtelomere_sequences[five_prime_max_coords$long_seqnames],
                        start = five_prime_max_coords$start_coord,
                        end = five_prime_max_coords$end_coord)


pdf("base_content_around_TSS.pdf")

consensusMatrix(TSS_sequences, as.prob = TRUE)[1:4,] %>%
   as.data.frame() %>%
   mutate(nucleotide = rownames(.)) %>%
   gather(key = "pos", value = "prob", -nucleotide) %>%
   mutate(pos = as.integer(str_remove(pos, "V"))) %>%
   ggplot(aes(x = pos, y = prob, fill = nucleotide)) + 
      geom_col(position = "stack")

consensusMatrix(TSS_sequences, as.prob = TRUE)[1:4,] %>%
   makePWM() %>%
   seqLogo()

consensusMatrix(TSS_sequences, as.prob = TRUE)[1:4,] %>%
   makePWM() %>%
   seqLogo(ic.scale = FALSE)

consensusMatrix(TSS_sequences, as.prob = TRUE)[1:4, 90:110] %>%
   makePWM() %>%
   seqLogo()

consensusMatrix(TSS_sequences, as.prob = TRUE)[1:4, 90:110] %>%
   makePWM() %>%
   seqLogo(ic.scale = FALSE)


gather.matrix <- function(mat) {
  if (is.null(dimnames(mat))) {
    grid <- expand.grid(seq.int(nrow(mat)), seq.int(ncol(mat)))
  } else {
    grid <- expand.grid(dimnames(mat))
  }
  cbind(grid, value = as.vector(mat))
}

base_frequencies_TSS <- consensusMatrix(TSS_sequences, as.prob = TRUE)[1:4, ]
colnames(base_frequencies_TSS) <- -100:100
base_frequencies_TSS <- gather.matrix(base_frequencies_TSS)
colnames(base_frequencies_TSS) <- c("Nucleotide", "Position", "Proportion")
base_frequencies_TSS$Position <- as.integer(as.character(base_frequencies_TSS$Position))

ggplot(base_frequencies_TSS,
       aes(x = Position,
           y = Proportion,
           group = Nucleotide,
           colour = Nucleotide)) +
   geom_line(alpha = I(1/3)) +
   geom_smooth() +
   scale_colour_viridis_d()


promoter_sequences <- subseq(subtelomere_sequences[five_prime_max_coords$long_seqnames],
                             start = five_prime_max_coords$start_coord,
                             end = five_prime_max_coords$end_coord + 2900)


base_frequencies_promoter <- consensusMatrix(promoter_sequences)[1:4, ]
colnames(base_frequencies_promoter) <- -100:3000
base_frequencies_promoter <- gather.matrix(base_frequencies_promoter)
colnames(base_frequencies_promoter) <- c("Nucleotide", "Position", "Counts")
base_frequencies_promoter$Position <- as.integer(as.character(base_frequencies_promoter$Position))

base_frequencies_promoter %>%
    spread(key = "Nucleotide", value = "Counts") %>%
    group_by(Position) %>%
    mutate(GC_content = sum(C + G)) %>%
    ggplot(aes(x = Position,
               y = GC_content)) +
      geom_line(alpha = I(1/3)) +
      geom_smooth()

ggplot(base_frequencies_promoter,
       aes(x = Position,
	   y = Counts,
	   group = Nucleotide,
	   colour = Nucleotide)) +
   geom_line(alpha = I(1/3)) +
   geom_smooth() +
   scale_colour_viridis_d()

dev.off()

library(TFBSTools)
library(JASPAR2014)

opts <- list()
opts[["species"]] <- 9606

jaspar_pfm <- getMatrixSet(JASPAR2014, opts)

jaspar_pwm <- toPWM(jaspar_pfm)

siteseqlist <- searchSeq(jaspar_pwm, promoter_sequences, min.score = "95%")

siteseqframe <- as(siteseqlist, "data.frame")

siteseqframe %>%
    ggplot(aes(x = TF)) +
      geom_bar() +
      theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1))

siteseqframe %>%
   ggplot(aes(x = start)) + geom_bar() + facet_grid(seqnames ~ TF)

siteseqframe %>%
   ggplot(aes(x = start, fill = TF)) + geom_bar() + facet_wrap(~ seqnames)

siteseqframe %>%
   ggplot(aes(x = start, fill = TF)) +
      geom_bar() +
      theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1))

siteseqframe %>%
   dplyr::count(TF, seqnames) %>%
   dplyr::count(TF) %>%
   ggplot(aes(x = TF, y = n)) +
      geom_col() +
      labs(title = "Subtelomeres with one PWM match",
           subtitle = "Searched in 30 subtelomeres that have a TERRA TSS") +
      theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1))

TFs_in_all_expressed <- siteseqframe %>%
   dplyr::count(TF, seqnames) %>%
   dplyr::count(TF) %>%
   filter(n == 30) %>%
   select(TF)
TFs_in_all_expressed <- as.character(TFs_in_all_expressed$TF)


siteseqframe %>%
   filter(TF %in% TFs_in_all_expressed) %>%
   ggplot(aes(x = start)) + geom_bar() + facet_grid(seqnames ~ TF)

siteseqframe %>%
   filter(TF %in% TFs_in_all_expressed) %>%
   ggplot(aes(x = start, fill = TF)) + geom_bar() + facet_wrap(~ seqnames)

siteseqframe %>%
   filter(TF %in% TFs_in_all_expressed) %>%
   ggplot(aes(x = start, fill = TF)) +
      geom_bar() +
      theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1))


dev.off()


