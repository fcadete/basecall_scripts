
library(Rsamtools)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

pdf("soft_clipped_sequence_analysis_telomereHunter.pdf", width = 16, height = 12)

subtelomere_end_bams <- list.files("soft_clipped_ends_telomereHunter/", pattern = "intratelomeric_subtelomere_start.bam")

full_frame <- data.frame()

for (file in subtelomere_end_bams) {

    this_bam <- scanBam(paste0("soft_clipped_ends_telomereHunter/", file))

    if (length(this_bam[[1]]$pos) > 0) {
        full_frame <- rbind(full_frame,
                            cbind(file, as.data.frame(this_bam)))
    }
}


full_frame <- full_frame %>%
                 mutate(cell_type = ifelse(grepl("HeL", file),
                               "HeLa",
                               ifelse(grepl("GM847", file),
                                      "GM847",
                                      ifelse(grepl("SAOS2", file),
                                             "SAOS2",
                                             ifelse(grepl("HEK", file),
                                                    "HEK293",
                                                    "U2OS")))),
                        adapter = str_extract(file, "Joana_[:alnum:]+"),
                        adapter = ifelse(is.na(adapter), "none", adapter),
                        sample = ifelse(str_detect(file, "VNP-TERRA_purified") |
                                         str_detect(file, "VNP_purified_TERRA"),
                                        "TERRA_purified",
                                        ifelse(str_detect(file, c("VNP_TERRA")),
                                               "TERRA",
                                               ifelse(str_detect(file, "VNP_polya"),
                                                      "polyA",
                                                      ifelse(str_detect(file, "VNP-pA-TERRA_purified") |
                                                              str_detect(file, "VNP_purified_pA_TERRA"),
                                                             "pA_TERRA_purified",
                                                             ifelse(str_detect(file, "VNP_pA-TERRA"),
                                                                       "pA_TERRA",
                                                                        NA))))))

full_frame <- full_frame %>%
   group_by(rname) %>%
   mutate(min_pos = min(pos)) %>%
   ungroup() %>%
   filter(abs(pos - min_pos) <= 100) %>%
   mutate(pos_to_min_pos = abs(pos - min_pos),
          soft_clipped_start = str_match(cigar, "\\d+S") %>% str_remove("S") %>% as.integer(),
          soft_clipped_seq = str_sub(seq, 1, soft_clipped_start - pos_to_min_pos),
          soft_clipped_quals = str_sub(qual, 1, soft_clipped_start - pos_to_min_pos),
          soft_clipped_seq_padded = str_pad(soft_clipped_seq, width = 2000, side = "left", pad = "N"),
          soft_clipped_seq_right_padded = str_pad(soft_clipped_seq, width = 2000, side = "right", pad = "N"))


full_frame <- full_frame %>%
   separate(rname, into = c("rname")) %>%
   mutate(rname = fct_reorder(rname, as.integer(str_extract(rname, "\\d+"))))

#p <- full_frame %>%
#   select(rname, soft_clipped_seq_padded) %>%
#   drop_na() %>%
#   group_by(rname) %>%
#   filter(n() >= 10) %>%
#   ungroup() %>%
#   separate(soft_clipped_seq_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
#   gather(key = "pos", value = "nucleotide", -rname, convert = TRUE) %>%
#   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
#   filter(pos >= 1750, nucleotide != "N") %>%
#   ggplot(aes(x = pos, fill = nucleotide, colour = nucleotide)) +
#   geom_bar() +
#   scale_fill_viridis_d() +
#   facet_wrap(~ rname, scales = "free_y" ) +
#   theme_bw() +
#   labs(title = "Number of soft-clipped nucleotides by subtelomere",
#        subtitle = "Aligned on the right, at the end of the subtelomere")
#print(p)
#
p <- full_frame %>%
   select(qname, rname, soft_clipped_seq_padded) %>%
   drop_na() %>%
   arrange(stringi::stri_reverse(soft_clipped_seq_padded)) %>%
   group_by(rname) %>%
   filter(n() >= 10) %>%
   mutate(index = row_number()) %>%
   ungroup() %>%
   separate(soft_clipped_seq_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -rname, -qname, -index, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos >= 1750) %>%
   ggplot(aes(x = pos, y = index, fill = nucleotide)) +
   geom_tile() +
   facet_wrap(~ rname, scales = "free_y") +
   theme_bw() +
   theme(axis.title.y = element_blank(),
         axis.text.y = element_blank()) +
   labs(title = "Heatmap of soft-clipped nucleotides by subtelomere",
        subtitle = "Aligned on the right, at the end of the subtelomere\nEach line corresponds to one read")
print(p)

p <- full_frame %>%
   select(qname, rname, soft_clipped_seq_padded) %>%
   drop_na() %>%
   arrange(stringi::stri_reverse(soft_clipped_seq_padded)) %>%
   group_by(rname) %>%
   filter(n() >= 10) %>%
   mutate(index = row_number()) %>%
   ungroup() %>%
   separate(soft_clipped_seq_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -rname, -qname, -index, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos >= 1950) %>%
   ggplot(aes(x = pos, y = index, fill = nucleotide)) +
   geom_tile() +
   facet_wrap(~ rname, scales = "free_y") +
   theme_bw() +
   theme(axis.title.y = element_blank(),
         axis.text.y = element_blank()) +
   labs(title = "Heatmap of soft-clipped nucleotides by subtelomere",
        subtitle = "Aligned on the right, at the end of the subtelomere\nEach line corresponds to one read")

   
print(p)

p <- full_frame %>%
   select(qname, rname, soft_clipped_seq_padded) %>%
   drop_na() %>%
   arrange(stringi::stri_reverse(soft_clipped_seq_padded)) %>%
   group_by(rname) %>%
   filter(n() >= 10) %>%
   mutate(index = row_number()) %>%
   ungroup() %>%
   separate(soft_clipped_seq_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -rname, -qname, -index, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos >= 1975) %>%
   ggplot(aes(x = pos, y = index, fill = nucleotide)) +
   geom_tile() +
   facet_wrap(~ rname, scales = "free_y") +
   theme_bw() +
   theme(axis.title.y = element_blank(),
         axis.text.y = element_blank()) +
   labs(title = "Heatmap of soft-clipped nucleotides by subtelomere",
        subtitle = "Aligned on the right, at the end of the subtelomere\nEach line corresponds to one read")

print(p)

p <- full_frame %>%
   select(qname, rname, soft_clipped_seq_padded) %>%
   drop_na() %>%
   arrange(stringi::stri_reverse(soft_clipped_seq_padded)) %>%
   mutate(index = row_number()) %>%
   separate(soft_clipped_seq_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -rname, -qname, -index, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos >= 1975) %>%
   ggplot(aes(x = pos, y = index, fill = nucleotide)) +
   geom_tile() +
   theme_bw() +
   theme(axis.title.y = element_blank(),
         axis.text.y = element_blank()) +
   labs(title = "Heatmap of soft-clipped nucleotides",
        subtitle = "Aligned on the right, at the end of the subtelomere\nEach line corresponds to one read")

print(p)


p <- full_frame %>%
   select(qname, rname, sample, soft_clipped_seq_padded) %>%
   drop_na() %>%
   group_by(sample) %>%
   arrange(stringi::stri_reverse(soft_clipped_seq_padded)) %>%
   mutate(index = row_number()) %>%
   ungroup() %>%
   separate(soft_clipped_seq_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -rname, -qname, -index, -sample, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos >= 1975) %>%
   ggplot(aes(x = pos, y = index, fill = nucleotide)) +
   facet_wrap(~ sample, scales = "free_y") +
   geom_tile() +
   theme_bw() +
   theme(axis.title.y = element_blank(),
         axis.text.y = element_blank()) +
   labs(title = "Heatmap of soft-clipped nucleotides by sample type",
        subtitle = "Aligned on the right, at the end of the subtelomere\nEach line corresponds to one read")


print(p)

p <- full_frame %>%
   select(qname, rname, cell_type, soft_clipped_seq_padded) %>%
   drop_na() %>%
   group_by(cell_type) %>%
   arrange(stringi::stri_reverse(soft_clipped_seq_padded)) %>%
   mutate(index = row_number()) %>%
   ungroup() %>%
   separate(soft_clipped_seq_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -rname, -qname, -index, -cell_type, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos >= 1975) %>%
   ggplot(aes(x = pos, y = index, fill = nucleotide)) +
   facet_wrap(~ cell_type, scales = "free_y") +
   geom_tile() +
   theme_bw() +
   theme(axis.title.y = element_blank(),
         axis.text.y = element_blank()) +
   labs(title = "Heatmap of soft-clipped nucleotides by cell type",
        subtitle = "Aligned on the right, at the end of the subtelomere\nEach line corresponds to one read")


print(p)



p <- full_frame %>%
   select(rname, soft_clipped_seq_padded) %>%
   drop_na() %>%
   group_by(rname) %>%
   filter(n() >= 10) %>%
   separate(soft_clipped_seq_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -rname, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos >= 1750,  nucleotide != "N") %>%
   ggplot(aes(x = pos, fill = nucleotide, colour = nucleotide)) +
   geom_bar(position = "fill") +
   scale_fill_viridis_d() +
   facet_wrap(~ rname, scales = "free_y") +
   theme_bw() +
   labs(title = "Proportion of soft-clipped nucleotides by subtelomere",
        subtitle = "Aligned on the right, at the end of the subtelomere")

print(p)

p <- full_frame %>%
   select(sample, soft_clipped_seq_padded) %>%
   drop_na() %>%
   separate(soft_clipped_seq_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -sample, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos >= 1750, nucleotide != "N") %>%
   ggplot(aes(x = pos, fill = nucleotide, colour = nucleotide)) +
   geom_bar() +
   scale_fill_viridis_d() +
   facet_wrap(~ sample, scales = "free_y") +
   theme_bw() +
   labs(title = "Number of soft-clipped nucleotides by sample type",
        subtitle = "Aligned on the right, at the end of the subtelomere")


print(p) 

p <- full_frame %>%
   select(sample, soft_clipped_seq_padded) %>%
   drop_na() %>%
   separate(soft_clipped_seq_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -sample, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos >= 1750, nucleotide != "N") %>%
   ggplot(aes(x = pos, fill = nucleotide, colour = nucleotide)) +
   geom_bar(position = "fill") +
   scale_fill_viridis_d() +
   facet_wrap(~ sample, scales = "free_y") +
   theme_bw() +
   labs(title = "Proportion of soft-clipped nucleotides by subtelomere",
        subtitle = "Aligned on the right, at the end of the subtelomere")


print(p)

p <- full_frame %>%
   select(sample, soft_clipped_seq_padded) %>%
   drop_na() %>%
   separate(soft_clipped_seq_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -sample, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos >= 1750, nucleotide != "N") %>%
   count(sample, pos, nucleotide) %>%
   ggplot(aes(x = pos, colour = nucleotide, y = n)) +
   geom_line() +
   scale_colour_viridis_d() +
   facet_wrap(~ sample, scales = "free_y") +
   theme_bw() +
   labs(title = "Number of soft-clipped nucleotides by sample",
        subtitle = "Aligned on the right, at the end of the subtelomere")


print(p)

p <- full_frame %>%
   select(sample, soft_clipped_seq_padded) %>%
   drop_na() %>%
   separate(soft_clipped_seq_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -sample, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos >= 1950, nucleotide != "N") %>%
   count(sample, pos, nucleotide) %>%
   ggplot(aes(x = pos, colour = nucleotide, y = n)) +
   geom_line() +
   scale_colour_viridis_d() +
   facet_wrap(~ sample, scales = "free_y") +
   theme_bw() +
   labs(title = "Number of soft-clipped nucleotides by sample",
        subtitle = "Aligned on the right, at the end of the subtelomere")


print(p)

p <- full_frame %>%
   select(sample, soft_clipped_seq_padded) %>%
   drop_na() %>%
   separate(soft_clipped_seq_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -sample, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos >= 1975, nucleotide != "N") %>%
   count(sample, pos, nucleotide) %>%
   ggplot(aes(x = pos, colour = nucleotide, y = n)) +
   geom_line() +
   scale_colour_viridis_d() +
   facet_wrap(~ sample, scales = "free_y") +
   theme_bw() +
   labs(title = "Number of soft-clipped nucleotides by sample",
        subtitle = "Aligned on the right, at the end of the subtelomere")


print(p)


# Now visualise the soft-clipped sequences aligning by their telomere-proximal end
p <- full_frame %>%
   select(sample, soft_clipped_seq_right_padded) %>%
   drop_na() %>%
   separate(soft_clipped_seq_right_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -sample, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(nucleotide != "N") %>%
   count(sample, pos, nucleotide) %>%
   ggplot(aes(x = pos, colour = nucleotide, y = n)) +
   geom_line() +
   scale_colour_viridis_d() +
   facet_wrap(~ sample, scales = "free_y") +
   theme_bw() +
   labs(title = "Number of soft-clipped nucleotides by sample",
        subtitle = "Aligned on the left, at the telomeric-prixmal end of the soft-clipped sequence")


print(p)

# Same as above, but zoomed in
p <- full_frame %>%
   select(sample, soft_clipped_seq_right_padded) %>%
   drop_na() %>%
   separate(soft_clipped_seq_right_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -sample, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos <= 500, nucleotide != "N") %>%
   count(sample, pos, nucleotide) %>%
   ggplot(aes(x = pos, colour = nucleotide, y = n)) +
   geom_line() +
   scale_colour_viridis_d() +
   facet_wrap(~ sample, scales = "free_y") +
   theme_bw() +
   labs(title = "Number of soft-clipped nucleotides by sample",
        subtitle = "Aligned on the left, at the telomeric-prixmal end of the soft-clipped sequence")

print(p)

# And even more zoomed in
p <- full_frame %>%
   select(sample, soft_clipped_seq_right_padded) %>%
   drop_na() %>%
   separate(soft_clipped_seq_right_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -sample, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos <= 100, nucleotide != "N") %>%
   count(sample, pos, nucleotide) %>%
   ggplot(aes(x = pos, colour = nucleotide, y = n)) +
   geom_line() +
   scale_colour_viridis_d() +
   facet_wrap(~ sample, scales = "free_y") +
   theme_bw() +
   labs(title = "Number of soft-clipped nucleotides by sample",
        subtitle = "Aligned on the left, at the telomeric-prixmal end of the soft-clipped sequence")

print(p)

p <- full_frame %>%
   select(qname, rname, soft_clipped_seq_right_padded) %>%
   drop_na() %>%
   arrange(soft_clipped_seq_right_padded) %>%
   group_by(rname) %>%
   filter(n() >= 10) %>%
   mutate(index = row_number()) %>%
   ungroup() %>%
   separate(soft_clipped_seq_right_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -rname, -qname, -index, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   ggplot(aes(x = pos, y = index, fill = nucleotide)) +
   geom_tile() +
   facet_wrap(~ rname, scales = "free_y") +
   theme_bw() +
   theme(axis.title.y = element_blank(),
         axis.text.y = element_blank()) +
   labs(title = "Heatmap of soft-clipped nucleotides by subtelomere",
        subtitle = "Aligned on the left, at the telomeric-prixmal end of the soft-clipped sequence")

print(p)

p <- full_frame %>%
   select(qname, rname, soft_clipped_seq_right_padded) %>%
   drop_na() %>%
   arrange(soft_clipped_seq_right_padded) %>%
   group_by(rname) %>%
   filter(n() >= 10) %>%
   mutate(index = row_number()) %>%
   ungroup() %>%
   separate(soft_clipped_seq_right_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -rname, -qname, -index, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos <= 500) %>%
   ggplot(aes(x = pos, y = index, fill = nucleotide)) +
   geom_tile() +
   facet_wrap(~ rname, scales = "free_y") +
   theme_bw() +
   theme(axis.title.y = element_blank(),
         axis.text.y = element_blank()) +
   labs(title = "Heatmap of soft-clipped nucleotides by subtelomere",
        subtitle = "Aligned on the left, at the telomeric-prixmal end of the soft-clipped sequence")

print(p)

p <- full_frame %>%
   select(qname, rname, soft_clipped_seq_right_padded) %>%
   drop_na() %>%
   arrange(soft_clipped_seq_right_padded) %>%
   group_by(rname) %>%
   filter(n() >= 10) %>%
   mutate(index = row_number()) %>%
   ungroup() %>%
   separate(soft_clipped_seq_right_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -rname, -qname, -index, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos <= 100) %>%
   ggplot(aes(x = pos, y = index, fill = nucleotide)) +
   geom_tile() +
   facet_wrap(~ rname, scales = "free_y") +
   theme_bw() +
   theme(axis.title.y = element_blank(),
         axis.text.y = element_blank()) +
   labs(title = "Heatmap of soft-clipped nucleotides by subtelomere",
        subtitle = "Aligned on the left, at the telomeric-prixmal end of the soft-clipped sequence")

print(p)

p <- full_frame %>%
   select(qname, rname, soft_clipped_seq_right_padded) %>%
   drop_na() %>%
   arrange(soft_clipped_seq_right_padded) %>%
   mutate(index = row_number()) %>%
   separate(soft_clipped_seq_right_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -rname, -qname, -index, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos <= 100) %>%
   ggplot(aes(x = pos, y = index, fill = nucleotide)) +
   geom_tile() +
   theme_bw() +
   theme(axis.title.y = element_blank(),
         axis.text.y = element_blank()) +
   labs(title = "Heatmap of soft-clipped nucleotides",
        subtitle = "Aligned on the left, at the telomeric-prixmal end of the soft-clipped sequence")

print(p)

p <- full_frame %>%
   select(qname, rname, sample, soft_clipped_seq_right_padded) %>%
   drop_na() %>%
   group_by(sample) %>%
   arrange(soft_clipped_seq_right_padded) %>%
   mutate(index = row_number()) %>%
   ungroup() %>%
   separate(soft_clipped_seq_right_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -rname, -qname, -index, -sample, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos <= 100) %>%
   ggplot(aes(x = pos, y = index, fill = nucleotide)) +
   geom_tile() +
   facet_wrap(~ sample, scales = "free_y") +
   theme_bw() +
   theme(axis.title.y = element_blank(),
         axis.text.y = element_blank()) +
   labs(title = "Heatmap of soft-clipped nucleotides by sample",
        subtitle = "Aligned on the left, at the telomeric-prixmal end of the soft-clipped sequence")

print(p)

p <- full_frame %>%
   select(qname, rname, cell_type, soft_clipped_seq_right_padded) %>%
   drop_na() %>%
   group_by(cell_type) %>%
   arrange(soft_clipped_seq_right_padded) %>%
   mutate(index = row_number()) %>%
   ungroup() %>%
   separate(soft_clipped_seq_right_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -rname, -qname, -index, -cell_type, convert = TRUE) %>%
   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
   filter(pos <= 100) %>%
   ggplot(aes(x = pos, y = index, fill = nucleotide)) +
   geom_tile() +
   facet_wrap(~ cell_type, scales = "free_y") +
   theme_bw() +
   theme(axis.title.y = element_blank(),
         axis.text.y = element_blank()) +
   labs(title = "Heatmap of soft-clipped nucleotides by cell type",
        subtitle = "Aligned on the left, at the telomeric-prixmal end of the soft-clipped sequence")

print(p)





dev.off()
#
#clipped_seq_50bp_end_dists <- full_frame$soft_clipped_seq_padded %>%
#   str_sub(1950, 2000) %>%
#   stringdistmatrix(method = "jaccard")
#
#clipped_seq_50bp_end_dists[is.na(clipped_seq_50bp_end_dists)] <- 0
#
#clipped_seq_50bp_end_clust <- hclust(clipped_seq_50bp_end_dists)
#
#
#full_frame %>%
#   select(qname, rname, soft_clipped_seq_padded) %>%
#   mutate(order_clustered = clipped_seq_50bp_end_clust$order) %>%
#   drop_na() %>%
#   arrange(stringi::stri_reverse(soft_clipped_seq_padded)) %>%
#   mutate(index = row_number()) %>%
#   separate(soft_clipped_seq_padded, into = as.character(c(1:2000)), sep = "(?<=.)") %>%
#   gather(key = "pos", value = "nucleotide", -rname, -qname, -index, -order_clustered, convert = TRUE) %>%
#   mutate(nucleotide = factor(nucleotide, levels = c("N", "C", "A", "T", "G"), ordered = TRUE)) %>%
#   filter(pos >= 1975) %>%
#   ggplot(aes(x = pos, y = order_clustered, fill = nucleotide)) +
#   geom_tile() +
#   theme_bw() +
#   theme(axis.title.y = element_blank(),
#         axis.text.y = element_blank()) +
#   labs(title = "Heatmap of soft-clipped nucleotides",
#        subtitle = "Aligned on the right, at the end of the subtelomere\nEach line corresponds to one read\ny-axis sorted by clustering")
#
#


