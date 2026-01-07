create_sci_names <- function(SS_tbl, su.column = "SiteSeriesLongName", spp.master){
  setDT(SS_tbl)
  su_var <- sym(su.column)
  taxon.info <- spp.master
  SiteSeriesLongName <- SS_tbl$su_var
  nms <- SS_tbl[,.(SiteSeriesLongName)]##only need long name
  nms[SiteSeriesLongName != "",]
  nms$SiteSeriesLongName <- str_replace(nms$SiteSeriesLongName, " - ", " – ")
  #nms <- as.data.table(nms)
  test <- as.data.table(tstrsplit(nms$SiteSeriesLongName, split = " – ", fixed = T))
  V4 <- "V4"
  if (!V4 %in% colnames(test)) {
    test[, (V4) := NA]
  }
  nms[, c("S1","S2","S3","S4") := test]  
  
 # nms[,c("S1","S2","S3","S4") := tstrsplit(SiteSeriesLongName," – ", fixed = T)]##need to split into components
  nms[,Trees := ifelse(nchar(S1) < 9,S1,NA)] ### move tree spps to new columns
  cols <- c("S1","S2","S3","S4")
  nms[,(cols) := lapply(.SD,tolower), .SDcols = cols] ##convert non-trees to lower
  #xx= nms
  lookup <- taxon.info[,.(TreeCode,ScientificName,EnglishName, ShortenedGuideName)]
  lookup <- lookup  %>% mutate(EnglishName = ifelse(!is.na(ShortenedGuideName), ShortenedGuideName, EnglishName))
  setDT(lookup)
  lookup[,EnglishName.lc := tolower(EnglishName)] ##also convert lookup to lower

  ###merge each column
  nms[lookup, Sc1 := i.ScientificName, on = c(S1 = "EnglishName.lc")]
  nms[lookup, Sc2 := i.ScientificName, on = c(S2 = "EnglishName.lc")]
  nms[lookup, Sc3 := i.ScientificName, on = c(S3 = "EnglishName.lc")]
  nms[lookup, Sc4 := i.ScientificName, on = c(S4 = "EnglishName.lc")]
  
  ##now deal with trees 
  lookup2 <- lookup[!is.na(TreeCode)]
  nms[,Trees := gsub("\\$","",Trees)]##remove $
  test <- as.data.table(tstrsplit(nms$Trees, split = "(?<=.)(?=[[:upper:]])", perl=TRUE))
  V3 <- "V3"
    if (!V3 %in% colnames(test)) {
      test[, (V3) := NA]
    }
  nms[, c("T1", "T2", "T3") := test]  
  #test <- as.data.table(tstrsplit(nms$Trees,split = "(?<=.)(?=[[:upper:]])"))
  #nms[,c("T1","T2") := tstrsplit(Trees,split = "(?<=.)(?=[[:upper:]])", perl = T)]##split before capital letters
  #xx <- nms
  ##merge tree codes
  nms[lookup2, Ts1 := i.ScientificName, on = c(T1 = "TreeCode")]
  nms[lookup2, Ts2 := i.ScientificName, on = c(T2 = "TreeCode")]
  nms[lookup2, Ts3 := i.ScientificName, on = c(T3 = "TreeCode")]
  
  ### check to see if there are any missing values
  nms[, missed_spp := FALSE]
  nms[is.na(Ts1) & !is.na(T1), missed_spp := TRUE]
  nms[is.na(Ts2) & !is.na(T2), missed_spp := TRUE]
  nms[is.na(Ts3) & !is.na(T3), missed_spp := TRUE]
  #nms[is.na(Sc1) & !is.na(S1), missed_spp := TRUE]
  nms[is.na(Sc2) & !is.na(S2), missed_spp := TRUE]
  nms[is.na(Sc3) & !is.na(S3), missed_spp := TRUE]
  nms[is.na(Sc4) & !is.na(S4), missed_spp := TRUE]
  
  # Check if any flag is TRUE and return the SiteUnit names
  if (any(nms$missed_spp)) {
    site_units <- nms[missed_spp == TRUE, ..su.column]
    merged_table <- merge(SS_tbl, site_units,  all.x = FALSE)
    merged_table <- merged_table[order(SS_NoSpace)]
    message("The following units have missed scientific names ", paste(merged_table$SS_NoSpace, collapse = ", "))
  }
  ###finally,
  # nms <- SS_tbl[,.(SiteSeriesLongName)]##only need long name
  # nms <- nms[SiteSeriesLongName != "",]
  # nms$SiteSeriesLongName <- str_replace(nms$SiteSeriesLongName, " - ", " – ")
  # nms <- as.data.table(nms)
  # nms[,c("S1","S2","S3","S4") := tstrsplit(SiteSeriesLongName," – ", fixed = T)]##need to split into components
  # nms[,Trees := ifelse(nchar(S1) < 7,S1,NA)] ### move tree spps to new columns
  # cols <- c("S1","S2","S3","S4")
  # nms[,(cols) := lapply(.SD,tolower), .SDcols = cols] ##convert non-trees to lower
  # 
  # lookup <- taxon.info[,.(TreeCode,ScientificName,EnglishName)] ##lookup table
  # lookup[,EnglishName := tolower(EnglishName)] ##also convert lookup to lower
  # ###merge each column
  # nms[lookup, Sc1 := i.ScientificName, on = c(S1 = "EnglishName")]
  # nms[lookup, Sc2 := i.ScientificName, on = c(S2 = "EnglishName")]
  # nms[lookup, Sc3 := i.ScientificName, on = c(S4 = "EnglishName")]
  # nms[lookup, Sc4 := i.ScientificName, on = c(S4 = "EnglishName")]
  # 
  # ##now deal with trees 
  # nms[,Trees := gsub("\\$","",Trees)]##remove $
  # nms[,c("T1","T2","T3") := tstrsplit(Trees,split = "(?<=.)(?=[[:upper:]])", perl = T)]##split before capital letters
  # ##merge tree codes
  # nms[lookup, Ts1 := i.ScientificName, on = c(T1 = "TreeCode")]
  # nms[lookup, Ts2 := i.ScientificName, on = c(T2 = "TreeCode")]
  # nms[lookup, Ts3 := i.ScientificName, on = c(T3 = "TreeCode")]
  ###finally, paste together
  nms <- unite(nms,col = "FinalScName", Ts1,Ts2,Ts3,Sc1,Sc2,Sc3,Sc4,na.rm = T, sep = " – ")
  SS_tbl$SiteSeriesScientificName <- nms$FinalScName
  return(SS_tbl)
}
