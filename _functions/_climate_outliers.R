climate_outliers <- function(
    plot_climate,
    su_anal,
    climate_vars,
    group_col = "tempAssoc_Code",
    unit_col = "SiteUnit",
    low_overlap_threshold = 0.10,
    alpha_global = .99,
    alpha_local = 0.90
) {
  
  # 1. Build association table
  assoc_climate <- plot_climate %>%
    left_join(su_anal, by = "PlotNumber") %>%
    select({{unit_col}}, all_of(climate_vars), groups = {{group_col}}) %>%
    filter(!is.na({{unit_col}}))
  
  # 2. Run climate comparison
  climate.compare <- climate_comparison(
    assoc_climate,
    group_col = "groups",
    unit_col = unit_col,
    vars = climate_vars,
    low_overlap_threshold = low_overlap_threshold
  )
  
  climate.assoc <- climate.compare$overlap_table
  
  # 3. Extract flagged SiteUnit–variable pairs
  flagged.units <- climate.compare$summary %>%
    filter(n_vars_flagged > 0) %>%
    drop_na(groups) %>%
    separate_rows(variables_flagged, sep = ",\\s*") %>%
    rename(variable = variables_flagged) %>%
    select(SiteUnit = {{unit_col}}, variable)
  
  # 4. Extract climate outliers
  climate_outliers <- climate.assoc %>%
    semi_join(flagged.units, by = c("SiteUnit", "variable")) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    select(groups, SiteUnit, variable, mean_grp, mean_su, everything())
  
  # 5. Build SiteUnit climate matrix
  su_matrix <- climate.assoc %>%
    group_by(SiteUnit, variable) %>%
    summarise(su_mean = mean(mean_su, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = variable, values_from = su_mean)
  
  # 6. Build group climate matrix
  group_matrix <- climate.assoc %>%
    group_by(groups, variable) %>%
    summarise(grp_mean = mean(mean_grp, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = variable, values_from = grp_mean)
  
  # 7. Covariance matrix for Mahalanobis
  cov_mat <- assoc_climate %>%
    select(all_of(climate_vars)) %>%
    stats::cov(use = "pairwise.complete.obs")
  
  # 8. Identify ungrouped SiteUnits
  ungrouped_su <- climate.assoc %>%
    filter(is.na(groups) | groups == "ungrouped") %>%
    pull(SiteUnit) %>%
    unique()
  
  # 9. Compute Mahalanobis distances
  library(dplyr)
  library(tidyr)
  library(FNN)
  library(purrr)
  
  # ---------------------------------------------------------
  # 1. BETWEEN-GROUP MAHALANOBIS DISTANCES
  # ---------------------------------------------------------
  
  # Expand all SU × group combinations
  mahal_results <- expand.grid(
    SiteUnit = ungrouped_su,
    groups   = group_matrix$groups,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      distance = purrr::pmap_dbl(
        list(SiteUnit, groups),
        function(su, g) {
          su_vec <- su_matrix %>%
            filter(SiteUnit == su) %>%
            select(-SiteUnit) %>%
            as.numeric()
          
          g_vec <- group_matrix %>%
            filter(groups == g) %>%
            select(-groups) %>%
            as.numeric()
          
          stats::mahalanobis(su_vec, g_vec, cov_mat)
        }
      )
    )
  
  climate_similarity <- mahal_results %>%
    group_by(SiteUnit) %>%
    arrange(distance) %>%
    mutate(rank_climate = row_number()) %>%
    ungroup() %>%
    rename(group = groups)
  
  flagged_su <- climate_outliers %>%
    pull(SiteUnit) %>%
    unique()
  
  mahal_flagged <- climate_similarity %>%
    filter(SiteUnit %in% flagged_su) %>%
    arrange(SiteUnit, distance)
  
  # ---------------------------------------------------------
  # 2. WITHIN-GROUP GLOBAL MAHALANOBIS
  # ---------------------------------------------------------
  
  # SU climate means
  su_climate <- climate.assoc %>%
    group_by(SiteUnit, groups, variable) %>%
    summarise(su_mean = mean(mean_su, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = variable, values_from = su_mean)
  
  # Group means + covariances
  group_stats <- assoc_climate %>%
    group_by(groups) %>%
    summarise(
      mean_vec = list(colMeans(across(all_of(climate_vars)), na.rm = TRUE)),
      cov_mat  = list(cov(across(all_of(climate_vars)),
                          use = "pairwise.complete.obs")),
      .groups = "drop"
    )
  
  # Join SU climate vectors with group stats
  su_with_group <- su_climate %>%
    left_join(group_stats, by = "groups") %>%
    drop_na(all_of(climate_vars))
  
  
  # ---------------------------------------------------------
  # 3. WITHIN-GROUP GLOBAL + LOCAL MAHALANOBIS
  # ---------------------------------------------------------
  
  # Precompute matrix for kNN
  clim_mat <- as.matrix(su_with_group[climate_vars])
  
  k <- 20
  nn_index <- get.knn(clim_mat, k = k)$nn.index
  
  su_within_mahal <- su_with_group %>%
    mutate(row_id = row_number()) %>%
    rowwise() %>%
    mutate(
      # -----------------------------
      # GLOBAL MAHALANOBIS
      # -----------------------------
      cov_reg = list({
        m <- cov_mat
        diag(m) <- diag(m) + 1e-6
        m
      }),
      distance_global = mahalanobis(
        x = c_across(all_of(climate_vars)),
        center = mean_vec,
        cov = cov_reg
      ),
      cutoff_global = qchisq(alpha_global, df = length(climate_vars)),
      is_outlier_global = distance_global > cutoff_global,
      
      # -----------------------------
      # LOCAL MAHALANOBIS
      # -----------------------------
      neigh_idx = list(nn_index[row_id, ]),
      neigh_mat = list(clim_mat[neigh_idx, , drop = FALSE]),
      
      local_mean = list(colMeans(neigh_mat)),
      local_cov = list({
        m <- cov(neigh_mat)
        diag(m) <- diag(m) + 1e-6
        m
      }),
      
      distance_local = mahalanobis(
        x = c_across(all_of(climate_vars)),
        center = local_mean,
        cov = local_cov
      ),
      cutoff_local = qchisq(alpha_local, df = length(climate_vars)),
      is_outlier_local = distance_local > cutoff_local
    ) %>%
    ungroup() %>%
    drop_na(groups)
  
  # Outliers (global or local)
  within_group_outliers <- su_within_mahal %>%
    filter(is_outlier_global | is_outlier_local) %>%
    arrange(groups, desc(distance_global))
  
  
  # Return everything cleanly
  list(
    assoc_climate = assoc_climate,
    climate_association = climate.assoc,
    flagged_units = flagged.units,
    climate_outliers = climate_outliers,
    su_matrix = su_matrix,
    group_matrix = group_matrix,
    covariance_matrix = cov_mat,
    mahal_all = climate_similarity,
    mahal_flagged = mahal_flagged,
    within_group_mahal = su_within_mahal,
    within_group_outliers = within_group_outliers
  )
}