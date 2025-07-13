library(ggplot2)
library(ggrepel)
library(reshape2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(corrplot)
library(tidyverse)
library(ggfortify)
library(ggpubr)
library(viridisLite)
library(pheatmap)

clean_colnames <- function(df) {
  df %>% rename_with(~str_replace_all(., "-", "_"))
}

df <- read_tsv("all_metrics_mincov10.tsv") %>%
  clean_colnames()

# Process count data
df_count <- df %>%
  mutate(merged_column = paste(chr, start, end, feature, score, strand, sep = ";")) %>%
  select(merged_column, matches("_count$")) %>%
  mutate(
    SKOV3_Tad_FTO_BT1_count = rowSums(select(., matches("SKOV3_Tad_FTO_BT1.*_count$")) %>% 
                                     mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Tad_FTO_UT_count = rowSums(select(., matches("SKOV3_Tad_FTO_UT.*_count$")) %>% 
                                   mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Apo_RIPK1_BT1_count = rowSums(select(., matches("SKOV3_Apo_RIPK1_BT1.*_count$")) %>% 
                                      mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Apo_RIPK1_UT_count = rowSums(select(., matches("SKOV3_Apo_RIPK1_UT.*_count$")) %>% 
                                     mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Apo_FTO_BT1_count = rowSums(select(., matches("SKOV3_Apo_FTO_BT1.*_count$")) %>% 
                                    mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Apo_FTO_UT_count = rowSums(select(., matches("SKOV3_Apo_FTO_UT.*_count$")) %>% 
                                   mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE)
  ) %>%
  select(merged_column, SKOV3_Tad_FTO_BT1_count, SKOV3_Tad_FTO_UT_count,
         SKOV3_Apo_RIPK1_BT1_count, SKOV3_Apo_RIPK1_UT_count,
         SKOV3_Apo_FTO_BT1_count, SKOV3_Apo_FTO_UT_count)
head(df_count)
write.csv(df_count, "rep_count.csv", row.names = FALSE)

# Process mean data
df_mean <- df %>%
  mutate(merged_column = paste(chr, start, end, feature, score, strand, sep = ";")) %>%
  select(merged_column, matches("_mean$")) %>%
  mutate(
    SKOV3_Tad_FTO_BT1_mean = rowSums(select(., matches("SKOV3_Tad_FTO_BT1.*_mean$")) %>% 
                                     mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Tad_FTO_UT_mean = rowSums(select(., matches("SKOV3_Tad_FTO_UT.*_mean$")) %>% 
                                   mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Apo_RIPK1_BT1_mean = rowSums(select(., matches("SKOV3_Apo_RIPK1_BT1.*_mean$")) %>% 
                                      mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Apo_RIPK1_UT_mean = rowSums(select(., matches("SKOV3_Apo_RIPK1_UT.*_mean$")) %>% 
                                     mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Apo_FTO_BT1_mean = rowSums(select(., matches("SKOV3_Apo_FTO_BT1.*_mean$")) %>% 
                                    mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Apo_FTO_UT_mean = rowSums(select(., matches("SKOV3_Apo_FTO_UT.*_mean$")) %>% 
                                   mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE)
  ) %>%
  select(merged_column, SKOV3_Tad_FTO_BT1_mean, SKOV3_Tad_FTO_UT_mean,
         SKOV3_Apo_RIPK1_BT1_mean, SKOV3_Apo_RIPK1_UT_mean,
         SKOV3_Apo_FTO_BT1_mean, SKOV3_Apo_FTO_UT_mean)
write.csv(df_mean, "rep_mean.csv", row.names = FALSE)

# Process sum data
df_sum <- df %>%
  mutate(merged_column = paste(chr, start, end, feature, score, strand, sep = ";")) %>%
  select(merged_column, matches("_sum$")) %>%
  mutate(
    SKOV3_Tad_FTO_BT1_sum = rowSums(select(., matches("SKOV3_Tad_FTO_BT1.*_sum$")) %>% 
                                     mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Tad_FTO_UT_sum = rowSums(select(., matches("SKOV3_Tad_FTO_UT.*_sum$")) %>% 
                                   mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Apo_RIPK1_BT1_sum = rowSums(select(., matches("SKOV3_Apo_RIPK1_BT1.*_sum$")) %>% 
                                      mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Apo_RIPK1_UT_sum = rowSums(select(., matches("SKOV3_Apo_RIPK1_UT.*_sum$")) %>% 
                                     mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Apo_FTO_BT1_sum = rowSums(select(., matches("SKOV3_Apo_FTO_BT1.*_sum$")) %>% 
                                    mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE),
    SKOV3_Apo_FTO_UT_sum = rowSums(select(., matches("SKOV3_Apo_FTO_UT.*_sum$")) %>% 
                                   mutate(across(everything(), ~as.numeric(as.character(.)))), na.rm = TRUE)
  ) %>%
  select(merged_column, SKOV3_Tad_FTO_BT1_sum, SKOV3_Tad_FTO_UT_sum,
         SKOV3_Apo_RIPK1_BT1_sum, SKOV3_Apo_RIPK1_UT_sum,
         SKOV3_Apo_FTO_BT1_sum, SKOV3_Apo_FTO_UT_sum)
write.csv(df_sum, "rep_sum.csv", row.names = FALSE)

# Score1 for count values
df_count_score1 <- df_count %>%
  separate(merged_column, 
          into = c("chr", "start", "end", "feature", "score", "strand"),
          sep = ";",
          remove = FALSE) %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    region_size = end - start
  ) %>%
  select(-c(chr, start, end, feature, score, strand)) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~as.numeric(as.character(.)))) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~replace_na(., 0))) %>%
  mutate(across(matches("UT"), ~pmax(., 0.1))) %>%
  mutate(
    ratio_RIPK1 = (SKOV3_Apo_RIPK1_BT1_count^2) / SKOV3_Apo_RIPK1_UT_count,
    ratio_ApoFTO = (SKOV3_Apo_FTO_BT1_count^2) / SKOV3_Apo_FTO_UT_count,
    ratio_TadFTO = (SKOV3_Tad_FTO_BT1_count^2) / SKOV3_Tad_FTO_UT_count,
    Combined_Score = if_else(ratio_RIPK1 > 0 & ratio_ApoFTO > 0 & ratio_TadFTO > 0,
                             ratio_RIPK1 * ratio_ApoFTO * ratio_TadFTO,
                             NA_real_)
  ) %>%
  select(-region_size, everything(), region_size) %>%
  arrange(desc(Combined_Score))
write_csv(df_count_score1, "score1_count_v1.csv")

# Score1 for mean values
df_mean_score1 <- df_mean %>%
  separate(merged_column, 
          into = c("chr", "start", "end", "feature", "score", "strand"),
          sep = ";",
          remove = FALSE) %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    region_size = end - start
  ) %>%
  select(-c(chr, start, end, feature, score, strand)) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~as.numeric(as.character(.)))) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~replace_na(., 0))) %>%
  mutate(across(matches("UT"), ~pmax(., 0.1))) %>%
  mutate(
    ratio_RIPK1 = (SKOV3_Apo_RIPK1_BT1_mean^2) / SKOV3_Apo_RIPK1_UT_mean,
    ratio_ApoFTO = (SKOV3_Apo_FTO_BT1_mean^2) / SKOV3_Apo_FTO_UT_mean,
    ratio_TadFTO = (SKOV3_Tad_FTO_BT1_mean^2) / SKOV3_Tad_FTO_UT_mean,
    Combined_Score = if_else(ratio_RIPK1 > 0 & ratio_ApoFTO > 0 & ratio_TadFTO > 0,
                             ratio_RIPK1 * ratio_ApoFTO * ratio_TadFTO,
                             NA_real_)
  ) %>%
  select(-region_size, everything(), region_size) %>%
  arrange(desc(Combined_Score))

write_csv(df_mean_score1, "score1_mean_v1.csv")

# Score1 for sum values
df_sum_score1 <- df_sum %>%
  separate(merged_column, 
          into = c("chr", "start", "end", "feature", "score", "strand"),
          sep = ";",
          remove = FALSE) %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    region_size = end - start
  ) %>%
  select(-c(chr, start, end, feature, score, strand)) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~as.numeric(as.character(.)))) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~replace_na(., 0))) %>%
  mutate(across(matches("UT"), ~pmax(., 0.1))) %>%
  mutate(
    ratio_RIPK1 = (SKOV3_Apo_RIPK1_BT1_sum^2) / SKOV3_Apo_RIPK1_UT_sum,
    ratio_ApoFTO = (SKOV3_Apo_FTO_BT1_sum^2) / SKOV3_Apo_FTO_UT_sum,
    ratio_TadFTO = (SKOV3_Tad_FTO_BT1_sum^2) / SKOV3_Tad_FTO_UT_sum,
    Combined_Score = if_else(ratio_RIPK1 > 0 & ratio_ApoFTO > 0 & ratio_TadFTO > 0,
                             ratio_RIPK1 * ratio_ApoFTO * ratio_TadFTO,
                             NA_real_)
  ) %>%
  select(-region_size, everything(), region_size) %>%
  arrange(desc(Combined_Score))

write_csv(df_sum_score1, "score1_sum_v1.csv")

# Score2 for count values
df_count_score2_log <- df_count %>%
  separate(merged_column, 
          into = c("chr", "start", "end", "feature", "score", "strand"),
          sep = ";",
          remove = FALSE) %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    region_size = end - start
  ) %>%
  select(-c(chr, start, end, feature, score, strand)) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~as.numeric(as.character(.)))) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~replace_na(., 0))) %>%
  mutate(
    ratio_RIPK1 = log2((SKOV3_Apo_RIPK1_BT1_count + 1) / (SKOV3_Apo_RIPK1_UT_count + 1)),
    ratio_ApoFTO = log2((SKOV3_Apo_FTO_BT1_count + 1) / (SKOV3_Apo_FTO_UT_count + 1)),
    ratio_TadFTO = log2((SKOV3_Tad_FTO_BT1_count + 1) / (SKOV3_Tad_FTO_UT_count + 1))
  ) %>%
  mutate(
    Combined_Score = if_else(ratio_RIPK1 > 0 & ratio_ApoFTO > 0 & ratio_TadFTO > 0,
                             ratio_RIPK1 * ratio_ApoFTO * ratio_TadFTO,
                             NA_real_)
  ) %>%
  select(-region_size, everything(), region_size) %>%
  arrange(desc(Combined_Score))

write_csv(df_count_score2_log, "score2_count_v1.csv")

# Score2 for mean values
df_mean_score2_log <- df_mean %>%
  separate(merged_column, 
          into = c("chr", "start", "end", "feature", "score", "strand"),
          sep = ";",
          remove = FALSE) %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    region_size = end - start
  ) %>%
  select(-c(chr, start, end, feature, score, strand)) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~as.numeric(as.character(.)))) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~replace_na(., 0))) %>%
  mutate(
    ratio_RIPK1 = log2((SKOV3_Apo_RIPK1_BT1_mean + 1) / (SKOV3_Apo_RIPK1_UT_mean + 1)),
    ratio_ApoFTO = log2((SKOV3_Apo_FTO_BT1_mean + 1) / (SKOV3_Apo_FTO_UT_mean + 1)),
    ratio_TadFTO = log2((SKOV3_Tad_FTO_BT1_mean + 1) / (SKOV3_Tad_FTO_UT_mean + 1))
  ) %>%
  mutate(
    Combined_Score = if_else(ratio_RIPK1 > 0 & ratio_ApoFTO > 0 & ratio_TadFTO > 0,
                             ratio_RIPK1 * ratio_ApoFTO * ratio_TadFTO,
                             NA_real_)
  ) %>%
  select(-region_size, everything(), region_size) %>%
  arrange(desc(Combined_Score))

write_csv(df_mean_score2_log, "score2_mean_v1.csv")

# Score2 for sum values
df_sum_score2_log <- df_sum %>%
  separate(merged_column, 
          into = c("chr", "start", "end", "feature", "score", "strand"),
          sep = ";",
          remove = FALSE) %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    region_size = end - start
  ) %>%
  select(-c(chr, start, end, feature, score, strand)) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~as.numeric(as.character(.)))) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~replace_na(., 0))) %>%
  mutate(
    ratio_RIPK1 = log2((SKOV3_Apo_RIPK1_BT1_sum + 1) / (SKOV3_Apo_RIPK1_UT_sum + 1)),
    ratio_ApoFTO = log2((SKOV3_Apo_FTO_BT1_sum + 1) / (SKOV3_Apo_FTO_UT_sum + 1)),
    ratio_TadFTO = log2((SKOV3_Tad_FTO_BT1_sum + 1) / (SKOV3_Tad_FTO_UT_sum + 1))
  ) %>%
  mutate(
    Combined_Score = if_else(ratio_RIPK1 > 0 & ratio_ApoFTO > 0 & ratio_TadFTO > 0,
                             ratio_RIPK1 * ratio_ApoFTO * ratio_TadFTO,
                             NA_real_)
  ) %>%
  select(-region_size, everything(), region_size) %>%
  arrange(desc(Combined_Score))

write_csv(df_sum_score2_log, "score2_sum_v1.csv")

# Score3 for count values
df_count_score3_zscore <- df_count %>%
  separate(merged_column, 
          into = c("chr", "start", "end", "feature", "score", "strand"),
          sep = ";",
          remove = FALSE) %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    region_size = end - start
  ) %>%
  select(-c(chr, start, end, feature, score, strand)) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~as.numeric(as.character(.)))) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~replace_na(., 0))) %>%
  mutate(
    ratio_RIPK1 = log2((SKOV3_Apo_RIPK1_BT1_count + 1) / (SKOV3_Apo_RIPK1_UT_count + 1)),
    ratio_ApoFTO = log2((SKOV3_Apo_FTO_BT1_count + 1) / (SKOV3_Apo_FTO_UT_count + 1)),
    ratio_TadFTO = log2((SKOV3_Tad_FTO_BT1_count + 1) / (SKOV3_Tad_FTO_UT_count + 1))
  ) %>%
  mutate(
    RIPK1_zscore = (ratio_RIPK1 - mean(ratio_RIPK1, na.rm = TRUE)) / sd(ratio_RIPK1, na.rm = TRUE),
    ApoFTO_zscore = (ratio_ApoFTO - mean(ratio_ApoFTO, na.rm = TRUE)) / sd(ratio_ApoFTO, na.rm = TRUE),
    TadFTO_zscore = (ratio_TadFTO - mean(ratio_TadFTO, na.rm = TRUE)) / sd(ratio_TadFTO, na.rm = TRUE)
  ) %>%
  mutate(
    Combined_Score = if_else(RIPK1_zscore > 0 & ApoFTO_zscore > 0 & TadFTO_zscore > 0,
                             RIPK1_zscore * ApoFTO_zscore * TadFTO_zscore,
                             NA_real_) 
  ) %>%
  select(-region_size, everything(), region_size) %>%
  arrange(desc(Combined_Score))

write_csv(df_count_score3_zscore, "score3_count_v1.csv")

# Score3 for mean values
df_mean_score3_zscore <- df_mean %>%
  separate(merged_column, 
          into = c("chr", "start", "end", "feature", "score", "strand"),
          sep = ";",
          remove = FALSE) %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    region_size = end - start
  ) %>%
  select(-c(chr, start, end, feature, score, strand)) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~as.numeric(as.character(.)))) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~replace_na(., 0))) %>%
  mutate(
    ratio_RIPK1 = log2((SKOV3_Apo_RIPK1_BT1_mean + 1) / (SKOV3_Apo_RIPK1_UT_mean + 1)),
    ratio_ApoFTO = log2((SKOV3_Apo_FTO_BT1_mean + 1) / (SKOV3_Apo_FTO_UT_mean + 1)),
    ratio_TadFTO = log2((SKOV3_Tad_FTO_BT1_mean + 1) / (SKOV3_Tad_FTO_UT_mean + 1))
  ) %>%
  mutate(
    RIPK1_zscore = (ratio_RIPK1 - mean(ratio_RIPK1, na.rm = TRUE)) / sd(ratio_RIPK1, na.rm = TRUE),
    ApoFTO_zscore = (ratio_ApoFTO - mean(ratio_ApoFTO, na.rm = TRUE)) / sd(ratio_ApoFTO, na.rm = TRUE),
    TadFTO_zscore = (ratio_TadFTO - mean(ratio_TadFTO, na.rm = TRUE)) / sd(ratio_TadFTO, na.rm = TRUE)
  ) %>%
  mutate(
    Combined_Score = if_else(RIPK1_zscore > 0 & ApoFTO_zscore > 0 & TadFTO_zscore > 0,
                             RIPK1_zscore * ApoFTO_zscore * TadFTO_zscore,
                             NA_real_) 
  ) %>%
  select(-region_size, everything(), region_size) %>%
  arrange(desc(Combined_Score))

write_csv(df_mean_score3_zscore, "score3_mean_v1.csv")

# Score3 for sum values
df_sum_score3_zscore <- df_sum %>%
  separate(merged_column, 
          into = c("chr", "start", "end", "feature", "score", "strand"),
          sep = ";",
          remove = FALSE) %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    region_size = end - start
  ) %>%
  select(-c(chr, start, end, feature, score, strand)) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~as.numeric(as.character(.)))) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~replace_na(., 0))) %>%
  mutate(
    ratio_RIPK1 = log2((SKOV3_Apo_RIPK1_BT1_sum + 1) / (SKOV3_Apo_RIPK1_UT_sum + 1)),
    ratio_ApoFTO = log2((SKOV3_Apo_FTO_BT1_sum + 1) / (SKOV3_Apo_FTO_UT_sum + 1)),
    ratio_TadFTO = log2((SKOV3_Tad_FTO_BT1_sum + 1) / (SKOV3_Tad_FTO_UT_sum + 1))
  ) %>%
  mutate(
    RIPK1_zscore = (ratio_RIPK1 - mean(ratio_RIPK1, na.rm = TRUE)) / sd(ratio_RIPK1, na.rm = TRUE),
    ApoFTO_zscore = (ratio_ApoFTO - mean(ratio_ApoFTO, na.rm = TRUE)) / sd(ratio_ApoFTO, na.rm = TRUE),
    TadFTO_zscore = (ratio_TadFTO - mean(ratio_TadFTO, na.rm = TRUE)) / sd(ratio_TadFTO, na.rm = TRUE)
  ) %>%
  mutate(
    Combined_Score = if_else(RIPK1_zscore > 0 & ApoFTO_zscore > 0 & TadFTO_zscore > 0,
                             RIPK1_zscore * ApoFTO_zscore * TadFTO_zscore,
                             NA_real_)
  ) %>%
  select(-region_size, everything(), region_size) %>%
  arrange(desc(Combined_Score))

write_csv(df_sum_score3_zscore, "score3_sum_v1.csv")

# Score4 for count values
df_count_score4_zscore <- df_count %>%
  separate(merged_column, 
          into = c("chr", "start", "end", "feature", "score", "strand"),
          sep = ";",
          remove = FALSE) %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    region_size = end - start
  ) %>%
  select(-c(chr, start, end, feature, score, strand)) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~as.numeric(as.character(.)))) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~replace_na(., 0))) %>%
  mutate(across(matches("UT"), ~pmax(., 0.1))) %>%  # More efficient than ifelse
  mutate(
    ratio_RIPK1 = (SKOV3_Apo_RIPK1_BT1_count^2) / SKOV3_Apo_RIPK1_UT_count,
    ratio_ApoFTO = (SKOV3_Apo_FTO_BT1_count^2) / SKOV3_Apo_FTO_UT_count,
    ratio_TadFTO = (SKOV3_Tad_FTO_BT1_count^2) / SKOV3_Tad_FTO_UT_count
  ) %>%
  mutate(
    RIPK1_zscore = (ratio_RIPK1 - mean(ratio_RIPK1, na.rm = TRUE)) / sd(ratio_RIPK1, na.rm = TRUE),
    ApoFTO_zscore = (ratio_ApoFTO - mean(ratio_ApoFTO, na.rm = TRUE)) / sd(ratio_ApoFTO, na.rm = TRUE),
    TadFTO_zscore = (ratio_TadFTO - mean(ratio_TadFTO, na.rm = TRUE)) / sd(ratio_TadFTO, na.rm = TRUE)
  ) %>%
  mutate(
    Combined_Score = if_else(RIPK1_zscore > 0 & ApoFTO_zscore > 0 & TadFTO_zscore > 0,
                             RIPK1_zscore * ApoFTO_zscore * TadFTO_zscore,
                             NA_real_) 
  ) %>%
  select(-region_size, everything(), region_size) %>%
  arrange(desc(Combined_Score))

write_csv(df_count_score4_zscore, "score4_count_v1.csv")

# Score4 for mean values
df_mean_score4_zscore <- df_mean %>%
  separate(merged_column, 
          into = c("chr", "start", "end", "feature", "score", "strand"),
          sep = ";",
          remove = FALSE) %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    region_size = end - start
  ) %>%
  select(-c(chr, start, end, feature, score, strand)) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~as.numeric(as.character(.)))) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~replace_na(., 0))) %>%
  mutate(across(matches("UT"), ~pmax(., 0.1))) %>%
  mutate(
    ratio_RIPK1 = (SKOV3_Apo_RIPK1_BT1_mean^2) / SKOV3_Apo_RIPK1_UT_mean,
    ratio_ApoFTO = (SKOV3_Apo_FTO_BT1_mean^2) / SKOV3_Apo_FTO_UT_mean,
    ratio_TadFTO = (SKOV3_Tad_FTO_BT1_mean^2) / SKOV3_Tad_FTO_UT_mean
  ) %>%
  mutate(
    RIPK1_zscore = (ratio_RIPK1 - mean(ratio_RIPK1, na.rm = TRUE)) / sd(ratio_RIPK1, na.rm = TRUE),
    ApoFTO_zscore = (ratio_ApoFTO - mean(ratio_ApoFTO, na.rm = TRUE)) / sd(ratio_ApoFTO, na.rm = TRUE),
    TadFTO_zscore = (ratio_TadFTO - mean(ratio_TadFTO, na.rm = TRUE)) / sd(ratio_TadFTO, na.rm = TRUE)
  ) %>%
  mutate(
    Combined_Score = if_else(RIPK1_zscore > 0 & ApoFTO_zscore > 0 & TadFTO_zscore > 0,
                             RIPK1_zscore * ApoFTO_zscore * TadFTO_zscore,
                             NA_real_) 
  ) %>%
  select(-region_size, everything(), region_size) %>%
  arrange(desc(Combined_Score))

write_csv(df_mean_score4_zscore, "score4_mean_v1.csv")

# Score4 for sum values
df_sum_score4_zscore <- df_sum %>%
  separate(merged_column, 
          into = c("chr", "start", "end", "feature", "score", "strand"),
          sep = ";",
          remove = FALSE) %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    region_size = end - start
  ) %>%
  select(-c(chr, start, end, feature, score, strand)) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~as.numeric(as.character(.)))) %>%
  mutate(across(starts_with(c("SKOV3_Apo_", "SKOV3_Tad_")), ~replace_na(., 0))) %>%
  mutate(across(matches("UT"), ~pmax(., 0.1))) %>%
  mutate(
    ratio_RIPK1 = (SKOV3_Apo_RIPK1_BT1_sum^2) / SKOV3_Apo_RIPK1_UT_sum,
    ratio_ApoFTO = (SKOV3_Apo_FTO_BT1_sum^2) / SKOV3_Apo_FTO_UT_sum,
    ratio_TadFTO = (SKOV3_Tad_FTO_BT1_sum^2) / SKOV3_Tad_FTO_UT_sum
  ) %>%
  mutate(
    RIPK1_zscore = (ratio_RIPK1 - mean(ratio_RIPK1, na.rm = TRUE)) / sd(ratio_RIPK1, na.rm = TRUE),
    ApoFTO_zscore = (ratio_ApoFTO - mean(ratio_ApoFTO, na.rm = TRUE)) / sd(ratio_ApoFTO, na.rm = TRUE),
    TadFTO_zscore = (ratio_TadFTO - mean(ratio_TadFTO, na.rm = TRUE)) / sd(ratio_TadFTO, na.rm = TRUE)
  ) %>%
  mutate(
    Combined_Score = if_else(RIPK1_zscore > 0 & ApoFTO_zscore > 0 & TadFTO_zscore > 0,
                             RIPK1_zscore * ApoFTO_zscore * TadFTO_zscore,
                             NA_real_) 
  ) %>%
  select(-region_size, everything(), region_size) %>%
  arrange(desc(Combined_Score))

write_csv(df_sum_score4_zscore, "score4_sum_v1.csv")
