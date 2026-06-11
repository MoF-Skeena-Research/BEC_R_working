climate_comparison <- function(df,
                                             group_col = "groups",
                                             unit_col = "SiteUnit",
                                             vars = c("DD5_an", "MCMT", "MSP", "PAS_an"),
                                             low_overlap_threshold = 0.10) {
  
  
  group_col <- rlang::sym(group_col)
  unit_col  <- rlang::sym(unit_col)
  
  # -----------------------------
  # 1. Compute overlap table
  # -----------------------------
  results <- df %>%
    group_by(!!group_col) %>%
    group_modify(~{
      df_group <- .x
      su_vals  <- unique(df_group[[rlang::as_string(unit_col)]])
      
      purrr::map_dfr(su_vals, function(su) {
        purrr::map_dfr(vars, function(v) {
          
          x <- df_group %>% filter(!!unit_col == su) %>% pull(!!sym(v)) %>% na.omit()
          y <- df_group %>% filter(!!unit_col != su) %>% pull(!!sym(v)) %>% na.omit()
          
          ov <- if (length(x) < 2 || length(y) < 2) NA 
          else overlapping::overlap(list(x, y))$OV
          
          tibble(
            group     = unique(df_group[[rlang::as_string(group_col)]]),
            SiteUnit  = su,
            variable  = v,
            overlap   = ov,
            mean_su   = mean(x),
            sd_su     = sd(x),
            mean_grp  = mean(y),
            sd_grp    = sd(y)
          )
        })
      })
    }) %>%
    ungroup()
  
  # -----------------------------
  # 2. Flag outliers
  # -----------------------------
  results_flagged <- results %>%
    group_by(groups, variable) %>%
    mutate(
      group_mean = mean(overlap, na.rm = TRUE),
      group_sd   = sd(overlap, na.rm = TRUE),
      flag_low_overlap = overlap < low_overlap_threshold,
      flag_sd = overlap < (group_mean - group_sd)
    ) %>%
    ungroup()#%>% 
  # group_by(SiteUnit, variable) %>% 
  #     mutate(su_mean = mean(overlap, na.rm = TRUE),
  #   su_sd   = sd(overlap, na.rm = TRUE))
  # 
  
  # -----------------------------
  # 3. Summary per SiteUnit (UPDATED)
  # -----------------------------
  summary_table <- results_flagged %>%
    group_by(groups, SiteUnit) %>%
    summarise(
      n_vars_flagged = sum(flag_low_overlap, na.rm = TRUE),
      variables_flagged = paste(variable[flag_low_overlap], collapse = ", "),
      mean_overlap = mean(overlap, na.rm = TRUE),
      sd_overlap   = sd(overlap, na.rm = TRUE),
      min_overlap  = min(overlap, na.rm = TRUE),
      max_overlap  = max(overlap, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    overlap_table = results_flagged,
    summary = summary_table
  )
}