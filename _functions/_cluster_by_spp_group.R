cluster_by_spp_group <- function(unit.compare,
                                 analysis_tab = NULL,
                                 group_var = NULL,
                                 group = NULL,
                                 cut.level = 0.7,
                                 combine.level = 0.5,
                                 include.low.diag = FALSE,
                                 low.diag.threshold = 25,
                                 include.few.plots = FALSE,
                                 few.plots.threshold = 5,
                                 output.file = "cluster_memberships.csv") {
  cluster_results <- list()
  all_memberships <- list()
  
  process_unit_compare <- function(unit_subset) {
    unit.compare2 <- unit_subset
    
    if (include.few.plots) {
      unit.few <- unit.compare2 %>%
        filter(nplots.x < few.plots.threshold) %>%
        select(Unit1, nplots.x) %>%
        distinct()
      
      unit.compare2 <- unit.compare2 %>%
        mutate(Unit1 = ifelse(Unit1 %in% unit.few$Unit1, paste0(Unit1, "-nplot"), Unit1),
               Unit2 = ifelse(Unit2 %in% unit.few$Unit1, paste0(Unit2, "-nplot"), Unit2))
    }
    
    if (include.low.diag) {
      unit.simple <- unit.compare2 %>%
        filter(unit.diag.sum.x < low.diag.threshold) %>%
        select(Unit1, unit.diag.sum.x) %>%
        distinct()
      
      unit.compare2 <- unit.compare2 %>%
        mutate(Unit1 = ifelse(Unit1 %in% unit.simple$Unit1, paste0(Unit1, "-lowdiag"), Unit1),
               Unit2 = ifelse(Unit2 %in% unit.simple$Unit1, paste0(Unit2, "-lowdiag"), Unit2))
    }
    
    unit.compare2 <- unit.compare2 %>%
      select(Unit1, Unit2, BEC.sim.min)
    
    return(unit.compare2)
  }
  
  compute_membership <- function(comp_data, group_label) {
    unit_labels <- unique(c(comp_data$Unit1, comp_data$Unit2))
    sim_matrix <- matrix(1, nrow = length(unit_labels), ncol = length(unit_labels),
                         dimnames = list(unit_labels, unit_labels))
    
    for (i in seq_len(nrow(comp_data))) {
      u1 <- comp_data$Unit1[i]
      u2 <- comp_data$Unit2[i]
      val <- comp_data$BEC.sim.min[i]
      sim_matrix[u1, u2] <- val
      sim_matrix[u2, u1] <- val
    }
    
    dissimilarity <- as.dist(1 - sim_matrix)
    hc <- hclust(dissimilarity, method = "average")
    groups <- cutree(hc, h = combine.level)
    
    return(data.frame(
      SiteUnit = names(groups),
      Cluster = groups,
      Spp.group = group_label,
      stringsAsFactors = FALSE
    ))
  }
  
  if (is.null(analysis_tab)) {
    message("🔄 No anal.groups provided — clustering all available units together...")
    unit_subset <- unit.compare
    unit.compare2 <- process_unit_compare(unit_subset)
    remaining_units <- unique(c(unit.compare2$Unit1, unit.compare2$Unit2))
    if (length(remaining_units) < 2) return(NULL)
    
    cluster_results[["All_Units"]] <- draw_dendro_split(unit.compare2, cut.level = cut.level)
    all_memberships[[1]] <- compute_membership(unit.compare2, "All_Units")
  } else {
    if (is.null(spp.groups)) {
      spp.groups <- unique(analysis_tab$group_var)
    }
    
    for (grp in spp.groups) {
      unit.choose <- analysis_tab %>%
        filter(group_var == group) %>%
        pull(SiteUnit) %>%
        unique()
      
      if (length(unit.choose) < 2) next
      
      unit_subset <- unit.compare %>%
        filter(Unit1 %in% unit.choose, Unit2 %in% unit.choose)
      
      unit.compare2 <- process_unit_compare(unit_subset)
      
      remaining_units <- unique(c(unit.compare2$Unit1, unit.compare2$Unit2))
      if (length(remaining_units) < 2) next
      
      cluster_results[[grp]] <- draw_dendro_split(unit.compare2, cut.level = cut.level)
      all_memberships[[grp]] <- compute_membership(unit.compare2, grp)
    }
  }
  
  membership_df <- do.call(rbind, all_memberships)# %>% mutate(Step1_Group = NA)
  write.csv(membership_df, file = output.file, row.names = FALSE)
  message("✅ Cluster membership file saved to: ", output.file)
  
  return(cluster_results)
}
