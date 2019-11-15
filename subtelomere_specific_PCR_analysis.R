
library(ShortRead)
library(tidyverse)

barcode_to_primer <- read_tsv("PCR_amplified/barcode_to_primers.tsv")

system("wc -l 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_verbatim/*/* > PCR_amplified/jo_primers_verbatim_total_reads_per_file.txt")
jo_primers_verbatim_total_reads <- read_table("PCR_amplified/jo_primers_verbatim_total_reads_per_file.txt",
                                              col_names = c("total_reads", "file")) %>%
                                      mutate(total_reads = total_reads / 4) %>%
                                      filter(file != "total") %>%
                                      mutate(barcode = str_extract(file, "BC\\d+"),
                                             barcode = ifelse(is.na(barcode), "none", barcode),
                                             primer_detected = str_extract(file, "SSP_.+"),
                                             primer_detected = str_remove(primer_detected, ".fastq"),
                                             primer_detected = ifelse(is.na(primer_detected), "none", primer_detected),
                                             primer_orientation = "verbatim") 


system("wc -l 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_revcomp/*/* > PCR_amplified/jo_primers_revcomp_total_reads_per_file.txt")
jo_primers_revcomp_total_reads <- read_table("PCR_amplified/jo_primers_revcomp_total_reads_per_file.txt",
                                              col_names = c("total_reads", "file"))  %>%
                                      mutate(total_reads = total_reads / 4)  %>%
                                      filter(file != "total") %>%
                                      mutate(barcode = str_extract(file, "BC\\d+"),
                                             barcode = ifelse(is.na(barcode), "none", barcode),
                                             primer_detected = str_extract(file, "SSP_.+"),
                                             primer_detected = str_remove(primer_detected, ".fastq"),
                                             primer_detected = ifelse(is.na(primer_detected), "none", primer_detected),
                                             primer_orientation = "revcomp") 



jo_primers_total_reads <- bind_rows(jo_primers_verbatim_total_reads,
                                    jo_primers_revcomp_total_reads) %>%
                             left_join(select(barcode_to_primer,
                                              barcode, primer_used = primer))

ggplot(jo_primers_total_reads %>% filter(barcode != "none",
                                         primer_detected != "none"),
       aes(x = primer_detected,
           y = total_reads,
           group = primer_orientation,
           fill = primer_detected == primer_used)) +
   geom_col() +
   facet_wrap(~ barcode)


ggplot(jo_primers_total_reads %>% filter(barcode != "none",
                                         primer_detected != "none"),
       aes(x = primer_detected,
           y = total_reads,
           fill = primer_detected == primer_used)) +
   geom_col() +
   facet_wrap(~ barcode + primer_orientation, ncol = 6)

load_primered_reads <- function(this_primer, this_barcode, verbatim_or_revcomp) {

      reads <- readFastq(paste0("191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_",
                         verbatim_or_revcomp, "/",
                         this_barcode, "/",
                         this_primer, ".fastq"))

      reads_as_frames <- tibble(barcode = this_barcode,
                                primer = this_primer,
                                orientation = verbatim_or_revcomp,
                                sequence = as.character(sread(reads)),
                                read_number = 1:length(reads))

      return(reads_as_frames)

   }



reads_verbatim_primer <- map2(barcode_to_primer$primer, barcode_to_primer$barcode,
                              load_primered_reads,
                              "verbatim") %>%
                            bind_rows()


reads_revcomp_primer <- map2(barcode_to_primer$primer, barcode_to_primer$barcode,
                              load_primered_reads,
                              "revcomp") %>%
                            bind_rows()

#
#reads_verbatim_primer %>%
#   bind_rows() %>%
#   separate(sequence,
#            sep = "(?<=.)",
#            into = as.character(1:max(nchar(reads_verbatim_primer$sequence))),
#            fill = "left") %>%
#   gather(key = "position", value = "nucleotide",
#          -barcode, -primer, -read_number) %>%
#   ggplot(aes(x = as.integer(position),
#              y = read_number,
#              fill = nucleotide)) +
#      geom_tile(colour = NA) +
#      facet_wrap(~ barcode + primer, scales = "free_y")
#
#reads_revcomp_primer %>%
#   bind_rows() %>%
#   separate(sequence,
#            sep = "(?<=.)",
#            into = as.character(1:max(nchar(reads_verbatim_primer$sequence))),
#            fill = "right") %>%
#   gather(key = "position", value = "nucleotide",
#          -barcode, -primer, -read_number) %>%
#   ggplot(aes(x = as.integer(position),
#              y = read_number,
#              fill = nucleotide)) +
#      geom_tile(colour = NA) +
#      facet_wrap(~ barcode + primer, scales = "free_y")
#

reads_verbatim_primer %>%
   bind_rows() %>%
   separate(sequence,
            sep = "(?<=.)",
            into = as.character(1:max(nchar(reads_verbatim_primer$sequence))),
            fill = "right") %>%
   gather(key = "position", value = "nucleotide",
          -barcode, -primer, -read_number, -orientation) %>%
   group_by(barcode, primer, read_number, orientation, position) %>%
   filter(n() >= 10) %>%
   ungroup() %>%
   ggplot(aes(x = as.integer(position),
              fill = nucleotide)) +
      geom_bar() +
      facet_wrap(~ barcode + primer, scales = "free")


reads_verbatim_primer %>%
   bind_rows() %>%
   separate(sequence,
            sep = "(?<=.)",
            into = as.character(1:max(nchar(reads_verbatim_primer$sequence))),
            fill = "right") %>%
   gather(key = "position", value = "nucleotide",
          -barcode, -primer, -read_number, -orientation) %>%
   drop_na() %>%
   group_by(barcode, primer, read_number, orientation, position) %>%
   filter(n() >= 10) %>%
   ungroup() %>%
   ggplot(aes(x = as.integer(position),
              fill = nucleotide)) +
      geom_bar(position = "fill", na.rm = TRUE) +
      facet_wrap(~ barcode + primer, scales = "free")



reads_revcomp_primer %>%
   bind_rows() %>%
   separate(sequence,
            sep = "(?<=.)",
            into = as.character(1:max(nchar(reads_revcomp_primer$sequence))),
            fill = "left") %>%
   gather(key = "position", value = "nucleotide",
          -barcode, -primer, -read_number, -orientation) %>%
   ggplot(aes(x = as.integer(position),
              fill = nucleotide)) +
      geom_bar() +
      facet_wrap(~ barcode + primer, scales = "free")


reads_revcomp_primer %>%
   bind_rows() %>%
   separate(sequence,
            sep = "(?<=.)",
            into = as.character(1:max(nchar(reads_revcomp_primer$sequence))),
            fill = "left") %>%
   gather(key = "position", value = "nucleotide",
          -barcode, -primer, -read_number, -orientation) %>%
   drop_na() %>%
   ggplot(aes(x = as.integer(position),
              fill = nucleotide)) +
      geom_bar(position = "fill", na.rm = TRUE) +
      facet_wrap(~ barcode + primer, scales = "free")




reads_all_primer <- bind_rows(reads_verbatim_primer %>% bind_rows(),
                              mutate(reads_revcomp_primer %>% bind_rows(),
                                     sequence = as.character(reverseComplement(DNAString(sequence)))))




reads <- readFastq("191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/BC01/SSP_1q.fastq")

reads_as_frame <- tibble(barcode = "BC01",
                          primer = "SSP_1q",
                          sequence = as.character(sread(reads)),
                          read_number = 1:length(reads))

reads_as_frame %>%
   arrange(sequence) %>%
   mutate(read_number = 1:n()) %>%
   separate(sequence,
            sep = "(?<=.)",
            into = as.character(1:max(nchar(reads_as_frame$sequence))),
            fill = "right") %>%
   gather(key = "position", value = "nucleotide",
          -barcode, -primer, -read_number) %>%
   drop_na() %>%
   filter(as.integer(position) <= 60) %>%
   ggplot(aes(x = as.integer(position),
              y = read_number,
              fill = nucleotide)) +
      geom_tile()




reads_50bp <- subseq(sread(reads), 1, 50)


reads_50bp_msa <- as.character(msa(reads_50bp))


reads_msa_as_frame <- tibble(barcode = "BC01",
                             primer = "SSP_1q",
                             sequence = paste0(reads_50bp_msa, subseq(sread(reads), 51)),
                             read_number = 1:length(reads))

reads_msa_as_frame %>%
   arrange(sequence) %>%
   mutate(read_number = 1:n()) %>%
   separate(sequence,
            sep = "(?<=.)",
            into = as.character(1:max(nchar(reads_msa_as_frame$sequence))),
            fill = "right") %>%
   gather(key = "position", value = "nucleotide",
          -barcode, -primer, -read_number) %>%
   drop_na() %>%
   ggplot(aes(x = as.integer(position),
              y = read_number,
              fill = nucleotide)) +
      geom_tile()


reads_as_frame %>%
   arrange(desc(reverse(str_replace_all(sequence, "AA[CTG]AA", "AAAAA")))) %>%
   mutate(read_number = 1:n()) %>%
   separate(sequence,
            sep = "(?<=.)",
            into = as.character(1:max(nchar(reads_as_frame$sequence))),
            fill = "left") %>%
   gather(key = "position", value = "nucleotide",
          -barcode, -primer, -read_number) %>%
   drop_na() %>%
   ggplot(aes(x = as.integer(position),
              y = read_number,
              fill = nucleotide)) +
      geom_tile()




load_reads <- function(this_primer, this_barcode, verbatim_or_revcomp) {

      reads <- readFastq(paste0("191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/",
                         this_barcode, "/",
                         this_primer, ".fastq"))

      reads <- reads[width(reads) >= 50]

      reads_50bp <- subseq(sread(reads), 1, 50)

      reads_50bp_msa <- as.character(msa(reads_50bp))


      reads_msa_as_frame <- tibble(barcode = this_barcode,
                                   primer = this_primer,
                                   sequence = as.character(sread(reads)),
                                   sequence_msa = paste0(reads_50bp_msa, subseq(sread(reads), 51)),
                                   read_number = 1:length(reads))


      return(reads_msa_as_frame)

   }


reads <- map2(barcode_to_primer$primer, barcode_to_primer$barcode,
                              load_reads) %>%
                            bind_rows()




for (this_barcode in unique(reads$barcode)) {

pdf(paste0("PCR_amplified/plots/", this_barcode, ".pdf"))

p <- reads %>%
   filter(barcode == this_barcode) %>%
   group_by(barcode) %>%
   arrange(barcode, desc(reverse(str_replace_all(sequence, "AA[CTG]AA", "AAAAA")))) %>%
   mutate(read_number = 1:n()) %>%
   select(-sequence_msa) %>%
   separate(sequence,
            sep = "(?<=.)",
            into = as.character(1:max(nchar(reads$sequence))),
            fill = "left") %>%
   gather(key = "position", value = "nucleotide",
          -barcode, -primer, -read_number) %>%
   drop_na() %>%
   filter((max(as.integer(position)) - as.integer(position)) <= 150) %>%
   ggplot(aes(x = as.integer(position),
              y = read_number,
              fill = nucleotide)) +
      geom_tile() +
      facet_wrap(~ barcode + primer, scales = "free")

print(p)

p <- reads %>%
   filter(barcode == this_barcode) %>%
   group_by(barcode) %>%
   arrange(barcode, sequence_msa) %>%
   mutate(read_number = 1:n()) %>%
   select(-sequence) %>%
   separate(sequence_msa,
            sep = "(?<=.)",
            into = as.character(1:max(nchar(reads$sequence_msa))),
            fill = "right") %>%
   gather(key = "position", value = "nucleotide",
          -barcode, -primer, -read_number) %>%
   drop_na() %>%
   filter(as.integer(position) <= 150) %>%
   ggplot(aes(x = as.integer(position),
              y = read_number,
              fill = nucleotide)) +
      geom_tile() +
      facet_wrap(~ barcode + primer, scales = "free")

print(p)


dev.off()

}



