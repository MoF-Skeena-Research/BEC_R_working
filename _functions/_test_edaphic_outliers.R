test_edaphic_outliers <- function(df, group_name) {
  library(dplyr)
  library(tidyr)
  library(purrr)
  
  #-------------------------------
  # Filter to selected group
  #-------------------------------
  df_g <- df %>%
    filter(groups == group_name) %>%
    filter(is.finite(SNR_num), is.finite(aSMR))
  
  #-------------------------------
  # Compute actual edatopic ranges + medians
  #-------------------------------
  ranges <- df_g %>%
    group_by(SiteUnit) %>%
    summarise(
      SNR_min     = min(SNR_num),
      SNR_median  = median(SNR_num),
      SNR_max     = max(SNR_num),
      aSMR_min    = min(aSMR),
      aSMR_median = median(aSMR),
      aSMR_max    = max(aSMR),
      .groups = "drop"
    )
  
  su_names <- ranges$SiteUnit
  
  #-------------------------------
  # Pairwise median comparison table
  #-------------------------------
  range_pairs <- expand.grid(
    SU1 = su_names,
    SU2 = su_names,
    stringsAsFactors = FALSE
  ) %>%
    filter(SU1 < SU2) %>%
    rowwise() %>%
    mutate(
      SU1_SNR_median  = ranges$SNR_median[ranges$SiteUnit == SU1],
      SU2_SNR_median  = ranges$SNR_median[ranges$SiteUnit == SU2],
      SU1_aSMR_median = ranges$aSMR_median[ranges$SiteUnit == SU1],
      SU2_aSMR_median = ranges$aSMR_median[ranges$SiteUnit == SU2]
    ) %>%
    ungroup()
  
  #-------------------------------
  # Group-level medians
  #-------------------------------
  group_medians <- ranges %>%
    summarise(
      group_SNR_median  = median(SNR_median),
      group_aSMR_median = median(aSMR_median)
    )
  
  #-------------------------------
  # Per-SiteUnit median flags
  #-------------------------------
  median_flags <- ranges %>%
    mutate(
      group_SNR_median  = group_medians$group_SNR_median,
      group_aSMR_median = group_medians$group_aSMR_median,
      flag_SNR  = abs(SNR_median  - group_SNR_median)  > 1,
      flag_aSMR = abs(aSMR_median - group_aSMR_median) > 1,
      flag_any  = flag_SNR | flag_aSMR
    )
  
  #-------------------------------
  # Return results
  #-------------------------------
  list(
    range_overlap = range_pairs,
    median_flags  = median_flags
  )
}