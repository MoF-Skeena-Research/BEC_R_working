combine_SU_databases <- function(folder_path) {
  db_files <- list.files(folder_path, pattern = "\\.accdb$", full.names = TRUE)  # Get all Access files
  
  SU_combined <- bind_rows(lapply(db_files, buildSU), .id = "source_db")  # Merge results
  
  return(SU_combined)
}