
library(tidyverse)

kmer_common_files <- list.files(pattern = "*kmer_common", recursive=TRUE)

kmer_table <- data.frame()

for (file in kmer_common_files) {

   purified_TERRA_kmers <- read_tsv(file, col_names = FALSE) %>%
      gather(-X1, key = "column", value = "kmer") %>%
      separate(kmer, into = c("kmer", "count"), convert = TRUE) %>%
      separate(X1, sep = " ",
               into = c("readid", "runid", "sampleid",
                        "readnumber", "channel", "starttime")) %>%
      select(-column)

      kmer_table <- rbind(kmer_table, purified_TERRA_kmers)

}

kmer_table %>%
   filter(kmer == "CCCTAA") %>%
   group_by(kmer, sampleid) %>%
   summarise(counts = sum(count, na.rm = TRUE),
             nreads = n(),
             counts_per_read = counts /nreads,
             nreadids = n_distinct(readid))

