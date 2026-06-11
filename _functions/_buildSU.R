buildSU <- function(db) {
  correlation <- dbConnect(odbc::odbc(), .connection_string = paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)}; DBQ=", db, ";"))
  
  su_tables <- dbListTables(correlation) %>% str_subset("_SU$")  # Get tables ending in "_SU"
  
  all_su <- lapply(setNames(nm = su_tables), dbReadTable, conn = correlation)  # Read tables
  
  dbDisconnect(correlation)  # Close connection
  
  if (length(all_su) > 0) {
    SU <- do.call(rbind.data.frame, all_su) %>%
      mutate(bgc = substr(SiteUnit, 1, 9)) %>%
      drop_na() %>%
      distinct(PlotNumber, .keep_all = TRUE) %>%
      filter(!grepl('poor|low|[$]|add|nudum|_[[:alpha:]]|X|omit|unplaced|moved|support', SiteUnit)) %>%
      arrange(desc(PlotNumber)) %>%
      select(PlotNumber, SiteUnit)  # Keep relevant columns
    
    return(SU)
  } else {
    return(NULL)  # Return NULL if no _SU tables exist
  }
}