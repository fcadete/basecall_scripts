
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)

pdf("subsample_results_2.pdf")

theme_set(theme_bw())

sample_percentage <- seq(5, 95, by = 5)
resample_index <- seq(1, 100)


#call_frame <- tibble()
#
#for (percentage in sample_percentage) {
#  for (resample in resample_index) {
#    if (file.exists(paste0("~/minion/subsample_analysis/porechop/",
#                           percentage, "/", resample,
#                           ".porechop_finalcall_table.gz"))) {
#       call_frame <- rbind(call_frame,
#                           cbind(percentage_sampled = percentage,
#                                 resample_index = resample,
#                                 read_tsv(paste0("~/minion/subsample_analysis/porechop/",
#                                                 percentage, "/", resample, ".porechop_finalcall_table.gz"),
#                                          col_names = c("readid", "runid", "sampleid",
#                                                        "read", "ch", "start_time",
#                                                        "final_barcode_call"),
#                                          skip = 1)))
#    }
#  }
#}
#
#call_frame %>%
#   count(percentage_sampled,
#         resample_index,
#         final_barcode_call) %>%
#   ggplot(aes(x = percentage_sampled,
#              colour = final_barcode_call,
#              y = n)) +
#     geom_jitter(height=0) +
#     scale_colour_viridis_d()
#
#call_frame %>%
#   count(percentage_sampled,
#         resample_index,
#         final_barcode_call) %>%
#   group_by(percentage_sampled, resample_index) %>%
#   mutate(total = sum(n)) %>%
#   ggplot(aes(x = percentage_sampled,
#              colour = final_barcode_call,
#              y = n / total)) +
#     geom_jitter(height=0) +
#     scale_colour_viridis_d()
#




depth_frame <- tibble()

for (percentage in sample_percentage) {
  for (resample in resample_index) {
    if (file.exists(paste0("~/minion/subsample_analysis/depths/",
                           percentage, "/", resample,
                           "/Joana_TERRA_VNP_no_PCR_on_rhietman_mapont_primary.depth.gz"))) {

       this_frame <- read_tsv(paste0("~/minion/subsample_analysis/depths/",
                                                 percentage, "/", resample,
                                                 "/Joana_TERRA_VNP_no_PCR_on_rhietman_mapont_primary.depth.gz"),
                                          col_names = c("chr", "pos", "depth")) %>% filter(pos <= 50000) 

       if (nrow(this_frame) > 0) {
          depth_frame <- rbind(depth_frame,
                               cbind(percentage_sampled = percentage,
                                     resample_index = resample,
                                     this_frame))
       }


       this_frame <- read_tsv(paste0("~/minion/subsample_analysis/depths/",
                                                 percentage, "/", resample,
                                                 "/Joana_TERRA_VNP_on_rhietman_mapont_primary.depth.gz"),
                                          col_names = c("chr", "pos", "depth")) %>% filter(pos <= 10000)

       if (nrow(this_frame) > 0) {
          depth_frame <- rbind(depth_frame,
                               cbind(percentage_sampled = percentage,
                                     resample_index = resample,
                                     this_frame))
       }


    }
  }
}


SUBTELOMERE_LENGTH <- 10000
SUBTELOMERE_NUMBER <- 42

depth_frame <- depth_frame %>%
                group_by(resample_index, percentage_sampled, chr, pos) %>%
                summarise(depth = sum(depth)) %>%
                separate(chr, into = c("chr"), extra = "drop")

depth_frame$chr <- factor(depth_frame$chr,
                          levels = c("1ptel", "1qtel", "2ptel", "2qtel", "3ptel",
                                     "3qtel", "4ptel", "4qtel", "5ptel", "5qtel",
                                     "6ptel", "6qtel", "7qtel", "7ptel", "8ptel", "8qtel",
                                     "9ptel", "9qtel", "10ptel", "10qtel", "11ptel", "11qtel",
                                     "12ptel", "12qtel", "13qtel", "14qtel", "15qtel", "16ptel",
                                     "16qtel", "17ptel", "17qtel", "18ptel", "18qtel", "19ptel",
                                     "19qtel", "20ptel", "20qtel", "21qtel", "22qtel",
                                     "XpYptel", "Xqtel", "Yqtel"),
                          ordered = TRUE)


depth_frame_summary <- depth_frame %>%
                        group_by(resample_index, percentage_sampled) %>%
                        summarise(number_subtelomeres_detected = n_distinct(chr),
                                  avg_depth_in_detected_subtelomeres = sum(depth)/(number_subtelomeres_detected * SUBTELOMERE_LENGTH),
                                  avg_depth_in_all_subtelomeres = sum(depth) / (SUBTELOMERE_NUMBER * SUBTELOMERE_LENGTH))


depth_frame_summary %>%
   ggplot(aes(x = percentage_sampled, y = number_subtelomeres_detected)) +
     geom_point(mapping = aes( colour = as.factor(percentage_sampled))) +
     geom_smooth() +
     scale_colour_viridis_d()

depth_frame_summary %>%
   ggplot(aes(x = percentage_sampled, y = avg_depth_in_detected_subtelomeres)) +
     geom_point(mapping = aes( colour = as.factor(percentage_sampled))) +
     geom_smooth() +
     scale_colour_viridis_d()

depth_frame_summary %>%
   ggplot(aes(x = percentage_sampled, y = avg_depth_in_all_subtelomeres)) +
     geom_point(mapping = aes( colour = as.factor(percentage_sampled))) +
     geom_smooth() +
     scale_colour_viridis_d()



depth_frame %>%
   group_by(resample_index, percentage_sampled, chr) %>%
   summarise(avg_depth_in_subtelomere = sum(depth)/SUBTELOMERE_LENGTH) %>%
   ggplot(aes(x = percentage_sampled, y = avg_depth_in_subtelomere, colour = as.factor(percentage_sampled))) +
     geom_jitter(height = 0) +
     scale_colour_viridis_d() +
     facet_wrap(~ chr)

dev.off()

