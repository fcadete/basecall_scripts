
library(tidyverse)
library(Rsamtools)
library(GenomicAlignments)
library(seqLogo)
library(TFBSTools)
library(JASPAR2014)

subtelomere_sequences <- readDNAStringSet("references/ConcatenatedFASTAAassemblies_hTel.txt",
                                          format = "FASTA")

subtelomere_starts <- substring(subtelomere_sequences, 1, 5000)

opts <- list()
opts[["species"]] <- 9606

jaspar_pfm <- getMatrixSet(JASPAR2014, opts)

jaspar_pwm <- toPWM(jaspar_pfm)

siteseqlist <- searchSeq(jaspar_pwm,
                         DNAStringSet(subtelomere_starts),
                         min.score = "95%")

siteseqframe <- as(siteseqlist, "data.frame")

write_tsv(siteseqframe, path = "JASPAR_TFBS_subtelomere_matches.tsv")

