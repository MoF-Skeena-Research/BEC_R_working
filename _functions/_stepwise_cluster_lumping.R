stepwise_cluster_lumping <- function(unit.compare,
                                     similarity_col = "BEC.sim.min",
                                     unit1_col = "Unit1",
                                     unit2_col = "Unit2",
                                     include.low.diag = FALSE,
                                     low.diag.threshold = 25,
                                     include.few.plots = FALSE,
                                     few.plots.threshold = 5,
                                     combine.level = 0.7,
                                     max.steps = 10,
                                     stop.no.change = TRUE,
                                     verbose = TRUE,
                                     initial.groups = NULL) {
  library(dplyr)
  library(tidyr)
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
  
  # Initialize unit history
  original_units <- unique(c(unit.compare[[unit1_col]], unit.compare[[unit2_col]]))
  
  if (!is.null(initial.groups)) {
    if (is.vector(initial.groups)) {
      unit_history <- data.frame(
        OriginalUnit = names(initial.groups),
        Step0_Group = as.character(initial.groups),
        stringsAsFactors = FALSE
      )
    } else if (is.data.frame(initial.groups) && all(c("SiteUnit", "Association") %in% names(initial.groups))) {
      unit_history <- initial.groups %>%
        select(OriginalUnit = SiteUnit, Step0_Group = Association) %>%
        mutate(Step0_Group = as.character(Step0_Group))
    } else {
      stop("initial.groups must be a named vector or a data frame with columns 'SiteUnit' and 'Association'")
    }
    
    # Add any missing units as their own group
    missing_units <- setdiff(original_units, unit_history$OriginalUnit)
    if (length(missing_units) > 0) {
      unit_history <- bind_rows(
        unit_history,
        data.frame(OriginalUnit = missing_units, Step0_Group = missing_units, stringsAsFactors = FALSE)
      )
    }
  } else {
    unit_history <- data.frame(
      OriginalUnit = original_units,
      Step0_Group = original_units,
      stringsAsFactors = FALSE
    )
  }
  
  prev_groups <- unit_history$Step0_Group
  
  for (step in 1:max.steps) {
    if (verbose) message("🔁 Step ", step)
    
    current_col <- paste0("Step", step - 1, "_Group")
    next_col <- paste0("Step", step, "_Group")
    
    # Map original units to current groups
    group_map <- unit_history %>%
      select(OriginalUnit, CurrentGroup = !!sym(current_col))
    
    # Join group labels to both Unit1 and Unit2
    grouped_compare <- unit.compare %>%
      left_join(group_map, by = setNames("OriginalUnit", unit1_col)) %>%
      rename(Group1 = CurrentGroup) %>%
      left_join(group_map, by = setNames("OriginalUnit", unit2_col)) %>%
      rename(Group2 = CurrentGroup)
    
    grouped_compare <- grouped_compare %>%
      filter(!is.na(Group1), !is.na(Group2))
    
    # Aggregate similarity between groups
    group_sim <- grouped_compare %>%
      group_by(Group1, Group2) %>%
      summarise(Sim = mean(.data[[similarity_col]], na.rm = TRUE), .groups = "drop")
    
    # Create symmetric similarity matrix
    group_labels <- unique(c(group_sim$Group1, group_sim$Group2))
    sim_matrix <- matrix(1, nrow = length(group_labels), ncol = length(group_labels),
                         dimnames = list(group_labels, group_labels))
    
    for (i in seq_len(nrow(group_sim))) {
      g1 <- group_sim$Group1[i]
      g2 <- group_sim$Group2[i]
      sim_matrix[g1, g2] <- group_sim$Sim[i]
      sim_matrix[g2, g1] <- group_sim$Sim[i]
    }
    
    # Convert to dissimilarity and cluster
    dissimilarity <- as.dist(1 - sim_matrix)
    hc <- hclust(dissimilarity, method = "average")
    new_groups <- cutree(hc, h = combine.level)
    
    # Create mapping of current group → new group
    group_map_new <- tibble::tibble(
      !!current_col := names(new_groups),
      !!next_col := paste0("Step", step, "_Grp", new_groups)
    )
    
    # Join to unit_history
    unit_history <- unit_history %>%
      left_join(group_map_new, by = current_col)
    
    current_groups <- unit_history[[next_col]]
    n_current <- length(unique(current_groups))
    
    if (stop.no.change && step > 1) {
      if (n_current == n_prev) {
        message("✅ Number of groups unchanged — stopping early at step ", step)
        break
      }
    }
    n_prev <- n_current
    
    if (verbose) {
      n_unique <- length(unique(current_groups))
      message("   → ", n_unique, " groups after step ", step)
    }
  }
  
  return(unit_history)
}