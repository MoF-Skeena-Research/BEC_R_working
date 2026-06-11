find_edaphic_similarity <- function(edatopic_table, df_raw) {
  library(dplyr)
  
  # 1. Group-level medians
  group_meds <- edatopic_table %>%
    group_by(group) %>%
    summarise(
      group_SNR_median  = median(SNR_median, na.rm = TRUE),
      group_aSMR_median = median(aSMR_median, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 2. Identify ungrouped SiteUnits
  ungrouped <- df_raw %>%
    filter(is.na(groups)) %>%
    distinct(SiteUnit)
  
  # 3. Compute distances
  out <- data.frame()
  
  for (su in ungrouped$SiteUnit) {
    
    su_vals <- df_raw %>%
      filter(SiteUnit == su) %>%
      summarise(
        SNR_median  = median(as.numeric(SNR), na.rm = TRUE),
        aSMR_median = median(aSMR, na.rm = TRUE)
      )
    
    for (i in seq_len(nrow(group_meds))) {
      g <- group_meds$group[i]
      
      d_snr  <- abs(su_vals$SNR_median  - group_meds$group_SNR_median[i])
      d_asmr <- abs(su_vals$aSMR_median - group_meds$group_aSMR_median[i])
      
      dist_score <- d_snr + d_asmr
      
      out <- rbind(out, data.frame(
        SiteUnit = su,
        group = g,
        d_SNR = d_snr,
        d_aSMR = d_asmr,
        dist_score = dist_score
      ))
    }
  }
  
  # 4. Rank groups for each ungrouped SU
  out_ranked <- out %>%
    group_by(SiteUnit) %>%
    arrange(dist_score) %>%
    mutate(rank_edatope = row_number()) %>%
    ungroup()
  
  out_ranked
}
