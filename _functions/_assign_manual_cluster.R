assign_manual_cluster <- function(su.tab,
                                  cluster.file,
                                  group_col = "Association",
                                  use_last_step = FALSE) {
  # Load cluster membership data
  if (is.character(cluster.file) && file.exists(cluster.file)) {
    membership <- data.table::fread(cluster.file)
  } else if (is.data.frame(cluster.file)) {
    membership <- cluster.file
  } else {
    stop("cluster.file must be either a valid file path or a data frame.")
  }
  
  # If using the last step from a stepwise lineage table
  if (use_last_step) {
    step_cols <- grep("^Step\\d+_Group$", names(membership), value = TRUE)
    if (length(step_cols) == 0) {
      stop("No stepwise group columns (e.g., Step1_Group) found in cluster data.")
    }
    last_col <- step_cols[which.max(as.numeric(gsub("\\D", "", step_cols)))]
    new_name_col <- last_col
  }
  
  # Ensure required columns exist
  if (!all(c("SiteUnit", group_col) %in% colnames(membership))) {
    stop("Required columns not found in cluster data.")
  }
  
  # Rename for clarity
  # membership <- membership %>%
  #   dplyr::rename(SiteUnit = all_of(siteunit_col),
  #                 Association = all_of(group_col))
  
  # Join and preserve original SiteUnit
  updated_data <- su.tab %>%
    dplyr::left_join(membership, by = "SiteUnit") %>% select(PlotNumber, "SiteUnit", all_of(group_col)) %>% rename(SiteUnit = 2, Step1_Group = 3) %>%  mutate(Step1_Group = ifelse(is.na(Step1_Group), SiteUnit, Step1_Group))
  
  return(updated_data)
}
