
library(Rsamtools)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(purrr)


pdf("soft_clipped_alignments_plots.pdf", height = 12, width = 24)

alns_folder <- "soft_clipped_ends_no_barcoding_q30/by_subtelomere/"

soft_clipped_alns <- list.files(alns_folder, pattern = ".aln")

alignments <- map_df(soft_clipped_alns,
                     function(aln_path) {
                          this_aln <- readDNAMultipleAlignment(paste0(alns_folder, aln_path), format = "FASTA")
                          this_aln <- as.character(this_aln)
                          this_aln <- tibble(file = aln_path, aln = this_aln)
                          return(this_aln)
                     }) %>%
   separate(file, sep = "VNP[-_]", into = c("subtelomere", "file")) %>%
   separate(subtelomere, sep = "_", into = c("subtelomere"), extra = "drop") %>%
   separate(file, sep = "_on_rhietman", into = c("file"), extra = "drop")



plot_frame <- alignments %>%
   group_by(file, subtelomere) %>%
   mutate(index = row_number()) %>%
   filter(n() > 10) %>%
   ungroup() %>%
   separate(aln, into = as.character(c(1:3000)), sep = "(?<=.)") %>%
   gather(key = "pos", value = "nucleotide", -file, -subtelomere, -index, convert = TRUE) %>%
   drop_na() %>%
   mutate(nucleotide = factor(nucleotide, levels = c("-", "C", "A", "T", "G"), ordered = TRUE))

plot_frame %>%
   ggplot(aes(x = pos, y = index, fill = nucleotide)) +
      geom_tile() +
      facet_wrap(~ subtelomere + file, scales = "free") +
      theme_bw()

for (this_file in unique(plot_frame$file)) {

   for(this_subtelomere in unique(plot_frame$subtelomere)) {

      p <- plot_frame %>%
          filter(file == this_file,
                 subtelomere == this_subtelomere) %>%
          ggplot(aes(x = pos, y = index, fill = nucleotide)) +
            geom_tile() +
            ggtitle(paste(this_file, this_subtelomere)) +
            theme_bw()

       print(p)

    }
}
         
dev.off()

