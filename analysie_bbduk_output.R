
library(tidyverse)
library(Shortread)

bbduk_files <- list.files(pattern = "*3joint_repeats.TTAGGG_counts", recursive=TRUE)

kmer_table <- data.frame()

for (file in bbduk_files) {
    
    reads <- readFastq(str_replace(file, "TTAGGG_counts", "fastq"))
    kmer_table <- rbind(kmer_table,
                        as.tibble(data.frame(read_table(file, col_names = c("counts", "readid")),
                                             sample = file,
                                             read_length = width(sread(reads)))))
}

kmer_table <- kmer_table %>% separate(sample, into = c("run", "sample"), sep = "/")

ggplot(kmer_table,
       aes(x = counts)) +
   geom_bar() +
   facet_wrap(~ sample)


ggplot(kmer_table,
       aes(x = read_length, y = counts, colour = sample)) +
   geom_point(alpha = 1/10) +
   facet_wrap(~ sample)


