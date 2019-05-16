
library(tidyverse)

call_files <- list.files(pattern = "./*porechop_finalcall_table")
call_table <- data_frame()
for (file in call_files) {
  call_table <- rbind(call_table,
                         cbind(file = file, read_tsv(file,
                                  col_names = c("readid", "runid", "sampleid",
                                                "read", "ch", "start_time",
                                                "final_barcode_call"),
                                  skip = 1)))
}


adapter_files <- list.files(pattern = "./*porechop_adapters_table")
adapters_table <- data_frame()
for (file in adapter_files) {
  adapters_table <- rbind(adapters_table,
                          read_tsv(file))
}


barcode_files <- list.files(pattern = "./*porechop_barcodes_table")
barcodes_table <- data_frame()
for (file in barcode_files) {
  barcodes_table <- rbind(barcodes_table,
                          read_tsv(file))
}

call_table %>%
  count(sampleid, final_barcode_call) %>%
  spread(key = sampleid, value = n)

call_table %>%
  count(sampleid, final_barcode_call) %>%
  group_by(sampleid) %>%
  mutate(prop = (n / sum(n)) * 100) %>%
  select(-n) %>%
  spread(key = sampleid, value = prop)

call_table %>%
  count(file, final_barcode_call) %>%
  spread(key = file, value = n)

call_table %>%
  count(file, final_barcode_call) %>%
  group_by(file) %>%
  mutate(prop = (n / sum(n)) * 100) %>%
  select(-n) %>%
  spread(key = file, value = prop)



barcodes_table <- barcodes_table %>%
  group_by(sampleid, readid) %>%
  arrange(sampleid, readid, desc(score)) %>%
  mutate(difference_to_next_adapter = score - lead(score))

barcodes_table %>% filter(score >= 60) %>%
  ggplot(aes(x = score, y = difference_to_next_adapter)) +
  geom_hex() +
  facet_wrap(~ sampleid)

top_barcodes_table <- barcodes_table %>%
  group_by(sampleid, readid) %>%
  filter(score == max(score)) %>%
  mutate(n_top_barcodes = n())

top_barcodes_table %>%
  filter(score > 60) %>%
  ggplot(aes(x = score, y = difference_to_next_adapter)) +
  geom_hex() +
  facet_wrap(~ sampleid)


top_barcodes_table %>%
  filter(score > 60) %>%
  count(sampleid, readid) %>%
  dim()

call_table %>%
  filter(final_barcode_call != "none") %>%
  dim()

top_barcodes_table %>%
  filter(score > 60) %>%
  count(sampleid, readid) %>%
  ungroup() %>%
  count(sampleid, n) %>%
  spread(key = n, value = nn)

left_join(top_barcodes_table %>%
            filter(n_top_barcodes == 1) %>%
            select(readid, sampleid, barcode, score),
          call_table %>%
            select(readid, final_barcode_call)) %>%
  filter(final_barcode_call == "none") %>%
  ggplot(aes(x = score, colour = sampleid)) +
  geom_density()


left_join(top_barcodes_table %>%
            select(readid, sampleid, barcode, score),
          call_table %>%
            select(readid, final_barcode_call)) %>%
  group_by(readid) %>%
  slice(1) %>%
  filter(final_barcode_call == "none") %>%
  ungroup() %>%
  count(sampleid, score > 60)

left_join(top_barcodes_table %>%
            filter(n_top_barcodes == 1) %>%
            select(readid, sampleid, barcode, score),
          call_table %>%
            select(readid, final_barcode_call)) %>%
  filter(final_barcode_call == "none") %>%
  ungroup() %>%
  count(sampleid, score > 60)

left_join(top_barcodes_table %>%
            filter(n_top_barcodes > 1) %>%
            select(readid, sampleid, barcode, score),
          call_table %>%
            select(readid, final_barcode_call)) %>%
  group_by(readid) %>%
  slice(1) %>%
  filter(final_barcode_call == "none") %>%
  ungroup() %>%
  count(sampleid, score > 60)

barcodes_table %>%
  ggplot(aes(x = score, colour = barcode)) +
  geom_vline(xintercept = 60) +
  geom_density() +
  facet_wrap(~ sampleid, scales = "free_y")
  
adapters_table %>%
  select(readid, sampleid, which_end, adapter_name, `full score`) %>%
  spread(key = adapter_name, value = `full score`) %>%
  ggplot(aes(x = `Barcode Joana_polyA_VNP (forward)`, y = `PCR adapters 1`)) +
  geom_hex() +
  facet_wrap(~ sampleid)
