##formats a vegetation summary table
#vsum = vegSum; spp =taxon.lifeform

format_veg_table_nostrata <- function(
    vsum = vegSum,
    spp = species,
    reorder_units = list(move = "_101$", before = "_110$")
) {
  
  #------------------------------------------------------------
  # 0. Lifeform lookup table (editable, explicit)
  #------------------------------------------------------------
  lifeform_map <- c(
    "1" = "Tree",
    "2" = "Tree",
    "3" = "Shrub",
    "4" = "Shrub",
    "12" = "Herb",
    "5" = "Herb",
    "6" = "Herb",
    "7" = "Herb",
    "8" = "Herb",
    "9" = "Moss",
    "10" = "Moss",
    "11" = "Moss"
    
  )
  
  lifeform_names <- c(
    "1" = "Conifer Tree",
    "2" = "Deciduous Tree",
    "3" = "Evergreen Shrub",
    "4" = "Deciduous Shrub",
    "12" = "Dwarf woody plants", 
    "6" = "Graminoids",
    "7" = "Forbs",
    "5" = "Ferns and Allies",   
    "8" = "Parasitic",
    "9" = "Mosses",
    "10" = "Liverworts",
    "11" = "Lichens",
    "13" = "Macro algae"
  )  
  
  strata_order <- c("Tree", "Shrub", "Herb", "Moss")
  #------------------------------------------------------------
  # Helper: reorder wide-format SiteUnit columns in data.table
  #------------------------------------------------------------
  .move_before_dt <- function(DT, move_suffix, before_suffix) {
    
    cols <- colnames(DT)
    
    move_cols   <- grep(move_suffix, cols, value = TRUE)
    before_cols <- grep(before_suffix, cols, value = TRUE)
    
    # If either group is missing, return unchanged
    if (length(move_cols) == 0 || length(before_cols) == 0)
      return(DT)
    
    cols_no_move <- cols[!cols %in% move_cols]
    pos_before   <- match(before_cols[1], cols_no_move)
    
    new_order <- append(cols_no_move, move_cols, after = pos_before - 1)
    
    DT[, ..new_order]
  }
  
  #------------------------------------------------------------
  # 1. Encode vegetation summary
  #------------------------------------------------------------
  # vsum[, Species2 := paste0(Species, "_", Layer)]
  # 
  # vsum[, code := encode_veg_sum(MeanCov, Constancy),
  #      by = .(Species2, SiteUnit)]
  
  vsum[, LifeformGroup := lifeform_map[as.character(Lifeform)]]
  
  vsum[, code := encode_veg_sum(MeanCov, Constancy)]
  
  vsum[, LifeformGroup := factor(LifeformGroup, levels = strata_order)]
  #------------------------------------------------------------
  # 2. Join species metadata
  #------------------------------------------------------------
  vsum <- merge(
    vsum,
    spp[, c("Code", "ScientificName", "EnglishName")],
    by.x = "Species",
    by.y = "Code",
    all.x = TRUE
  )
  
  
  #------------------------------------------------------------
  # 3. Count plots per SiteUnit
  #------------------------------------------------------------
  nPlots <- unique(vsum[, .(SiteUnit, nplots)])[order(nplots, decreasing = TRUE)]
  nPlots <- nPlots[order(SiteUnit)]
  
  
  #------------------------------------------------------------
  # 3.5 Detect duplicates BEFORE casting
  #------------------------------------------------------------
  dup_keys <- c("LifeformGroup", "ScientificName", "EnglishName", "SiteUnit")
  
  # summary of duplicated combinations
  dup_summary <- vsum[, .N, by = dup_keys][N > 1]
  
  # full rows for inspection
  dup_rows <- vsum[vsum[, .I[.N > 1], by = dup_keys]$V1]
  
  # (optional) message if duplicates exist
  if (nrow(dup_summary) > 0) {
    message("⚠️ Duplicate combinations detected before casting. Inspect `dup_summary` and `dup_rows`.")
  }
  
  
  #------------------------------------------------------------
  # 4. Cast to wide format
  #------------------------------------------------------------
  vsum <- data.table::dcast(
    vsum,
    LifeformGroup + ScientificName + EnglishName ~ SiteUnit,
    value.var = "code",
    fill = ""
  )[order(LifeformGroup, ScientificName)]
  
  #------------------------------------------------------------
  # 5. Rename columns
  #------------------------------------------------------------
  data.table::setnames(
    vsum,
    old = c("LifeformGroup", "ScientificName", "EnglishName"),
    new = c("Layer", "Scientific name", "Common name")
  )
  
  #------------------------------------------------------------
  # 6. Clean codes
  #------------------------------------------------------------
  vsum <- vsum %>% mutate_all(str_replace_all, "-remove", "")
  
  #------------------------------------------------------------
  # 7. Reorder rows by layer
  #------------------------------------------------------------
  vsum <- vsum %>%
    select(order(colnames(vsum))) %>%
    select(Layer, `Scientific name`, everything()) %>%
    relocate(`Common name`, .after = last_col()) %>%
    arrange(match(Layer, c("Tree", "Shrub", "Herb", "Moss")), Layer)
  
  #------------------------------------------------------------
  # 8. Add nPlots header row
  #------------------------------------------------------------
  nPlotRow <- c("", "n Plots", nPlots$nplots, "") %>%
    matrix(nrow = 1) %>%
    data.frame() %>%
    stats::setNames(names(vsum)) %>%
    data.table::as.data.table()
  
  vsum <- rbind(nPlotRow, vsum)
  
  #------------------------------------------------------------
  # 9. Optional: reorder SiteUnit columns (wide format)
  #------------------------------------------------------------
  if (!is.null(reorder_units)) {
    vsum <- .move_before_dt(
      vsum,
      move_suffix   = reorder_units$move,
      before_suffix = reorder_units$before
    )
  }
  return(vsum)
}

#############  OLD  #######################
# format_veg_table_strata <- function(vsum = vegSum, spp = species){
#   #create new variable Species2 in vegSum from Species and Layer
#   #has_layer    <- "Layer" %in% names(vsum)
#   # has_lifeform <- "Lifeform" %in% names(vsum)
#   # 
#   # # Join spp for Lifeform if missing
#   # if (!has_lifeform) {
#   #   message("Lifeform not found in input; deriving from spp lookup.")
# 
#   # vsum <- vsum %>% mutate(Layer = as.character(Lifeform))
#   # vsum$Layer <-  case_match(vsum$Layer, "1" ~ "Tree",
#   #                            "2" ~ "Tree",
#   #                            "3" ~ "Shrub",
#   #                            "4" ~ "Shrub",
#   #                            "5" ~ "Herb",
#   #                            "6" ~ "Herb" ,
#   #                            "7" ~ "Herb",
#   #                            "8" ~ "Herb",
#   #                            "9" ~ "Moss",
#   #                           "10" ~ "Moss",
#   #                           "11" ~ "Moss",
#   #                           "12" ~ "Herb")
#   vsum$Species2 <- paste0(vsum$Species, "_", vsum$Layer)
#   vsum[ , code := encode_veg_sum(MeanCov, Constancy), by = .(Species2, SiteUnit)]
#   
#   vsum <- vsum %>%
#     merge(
#       spp %>% select(Code, ScientificName, Lifeform, EnglishName),
#       by.x = "Species",
#       by.y = "Code",
#       all.x = TRUE
#     )
#   
#   #merge(spp, by.x = 'Species', by.y = 'Code') |>
#   #mutate(SiteUnit = str_replace(SiteUnit, "101", "109"))
#   nPlots <- unique(vsum[ ,.(SiteUnit, nplots)])[order(nplots, decreasing = TRUE), ] %>% 
#     arrange(SiteUnit)
#   vsum <- data.table::dcast(vsum, 
#                             Layer + ScientificName + EnglishName ~ SiteUnit, 
#                             value.var = 'code',
#                             fill = '')[order(Layer, ScientificName), ]
#   #vsum[duplicated(ReportName, fromLast = TRUE), ReportName := '']
#   vsum[ , c('Layer','ScientificName', nPlots$SiteUnit, 'EnglishName'), with = FALSE]
#   data.table::setnames(vsum, old = c('Layer', 'ScientificName', 
#                                      'EnglishName'), new = c('Layer', 'Scientific name', 'Common name'))
#   vsum <- vsum %>% 
#     mutate_all(str_replace_all,"-remove", "")
#   # vsum$Layer <-  case_match(vsum$Layer, "1" ~ "A",
#   #                            "2" ~ "A",
#   #                            "3" ~ "B",
#   #                            "4" ~ "B",
#   #                            "5" ~ "C",
#   #                            "6" ~ "C" ,
#   #                            "7" ~ "C",
#   #                            "8" ~ "C",
#   #                            "9" ~ "D",
#   #                           "10" ~ "D",
#   #                           "11" ~ "D",
#   #                           "12" ~ "C")
#   #vsum <- vsum[ order(match(vsum$Scientific, indic.order$ScientificName)), ]
#   if (exists("indic.order") && "ScientificName" %in% names(indic.order)) {
#     vsum <- vsum[ order(match(vsum$Scientific, indic.order$ScientificName)), ]
#   }
#   vsum2 <-   vsum %>% select(order(colnames(vsum))) %>%  
#    select(Layer, `Scientific name`, everything()) %>% 
#     relocate(`Common name`, .after = last_col()) %>%
#     arrange(match(Layer, c("Tree", "Regen", "Shrub", "Herb", "Moss")), Layer)
#   
#   #colnames(vsum2) <- gsub(paste0(bgc.choose,"?"), "", colnames(vsum2))
#   nPlotRow <- c('', 'n Plots', nPlots$nplots, '') |> 
#     matrix(nrow = 1) |>
#     data.frame() |> 
#     stats::setNames(names(vsum2)) |> 
#     data.table::as.data.table()
#   vsum3 <- rbind(nPlotRow, vsum2)
# }
