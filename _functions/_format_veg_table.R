##formats a vegetation summary table
##types of summary tables: BGC =  site series within a BGC; zonal = comparison of zonal vegtation between BGCs
#type = "Zonal"
vsum = vegSum; spp=taxon.lifeform; cons.1 = 70; cons.2 = 50; strata.by = "Auto"; type = NULL

format_veg_table <- function(
    vsum = vegSum,
    spp = taxon.all,
    type = "BGC",
    cons.1 = 70,
    cons.2 = 50,
    strata.by = c("Layer", "Lifeform", "Auto"),
    indic.order = NULL
){
  
  strata.by <- match.arg(strata.by)
  
  #------------------------------------------------------------
  # 0. Lifeform lookup table (editable, explicit)
  #------------------------------------------------------------
  lifeform_map <- c(
    "1" = "Tree",
    "2" = "Regen",
    "3" = "Shrub",
    "4" = "Shrub",
    "5" = "Herb",
    "6" = "Herb",
    "7" = "Herb",
    "8" = "Herb",
    "9" = "Moss",
    "10" = "Moss",
    "11" = "Moss",
    "12" = "Herb"
  )

  lifeform_names <- c(
    "1" = "Conifer Tree",
    "2" = "Deciduous Tree",
    "3" = "Evergreen Shrub",
    "4" = "Deciduous Shrub",
    "5" = "Ferns and Allies",
    "6" = "Graminoids",
    "7" = "Forbs",
    "8" = "Parasitic",
    "9" = "Mosses",
    "10" = "Liverworts",
    "11" = "Lichens",
    "12" = "Dwarf woody plants",
    "13" = "Macro algae"
  )  
  
  #------------------------------------------------------------
  # 1. Encode coverage + constancy into symbol codes
  #------------------------------------------------------------
 
   encode_veg_sum <- function(coverage, constancy) {
    black <- "n"; grey <- "l"; star <- "v"
    char <- black
    if (constancy < cons.1) char <- grey
    if (constancy < cons.2) char <- star
    
    color <- "remove"
    
    data.table::fcase(
      coverage <   1, sprintf("%s-%s", strrep(char, 1), color),
      coverage <   3, sprintf("%s-%s", strrep(char, 2), color),
      coverage <  10, sprintf("%s-%s", strrep(char, 3), color),
      coverage <  25, sprintf("%s-%s", strrep(char, 4), color),
      coverage < 100, sprintf("%s-%s", strrep(char, 5), color),
      default = strrep(char, 6)
    )
   }
  #------------------------------------------------------------
  # 1b. Grouping logic (Layer / Lifeform / Auto)
  #------------------------------------------------------------
  group_species <- function(df, mode) {
    if (mode == "Layer") {
      df$Group <- df$Layer
    } else if (mode == "Lifeform") {
      df$Group <- df$LifeformName
    } else if (mode == "Auto") {
      df$Group <- ifelse(df$Layer == "Shrub" & df$Lifeform %in% c(1,2),
                         "Regen",
                         df$Layer)
    } else {
      stop("Invalid strata.by value.")
    }
    df
  }
  
  #------------------------------------------------------------
  # 2. Prepare vsum and harmonize Layer / Lifeform inputs
  #------------------------------------------------------------
  has_layer    <- "Layer" %in% names(vsum)
  has_lifeform <- "Lifeform" %in% names(vsum)
  
  vsum <- vsum %>%
    merge(spp %>% select(ScientificName, EnglishName, Code), by.x = 'Species', by.y = 'Code')
  # Ensure Layer is character if present
  if (has_layer) {
    vsum <- vsum %>% mutate(Layer = as.character(Layer))
  }
  
  
  
  # Join spp for Lifeform if missing
  if (has_layer && !has_lifeform) {
    message("Lifeform not found in input; deriving from spp lookup.")
    # vsum <- vsum %>%
    #   merge(spp %>% select(ScientificName, Lifeform),
    #         by = "ScientificName",
    #         all.x = TRUE)
  }
  
  # Derive Layer from Lifeform if missing
  if (has_lifeform && !has_layer) {
    message("Layer not found in input; deriving from Lifeform mapping.")
    vsum <- vsum %>%
      mutate(Layer = lifeform_map[as.character(Lifeform)])
  }
  
  # Add LifeformName
  vsum <- vsum %>%
    mutate(LifeformName = lifeform_map[as.character(Lifeform)])
  
  if (any(is.na(vsum$LifeformName))) {
    warning("Some species have missing or unmapped Lifeform values. See saved missing_lifeform.csv for details.")
    missing_lf <- vsum %>%
      filter(is.na(LifeformName)) %>%
      select(Species, ScientificName, Lifeform) %>%
      distinct()
    fwrite(missing_lf, "missing_lifeform.csv")
    
  }
  
  # Apply strata.by BEFORE Species2 is created
  vsum <- group_species(vsum, strata.by)
  
  # Species2 reflects chosen grouping variable
  vsum <- vsum %>%
    mutate(Species2 = paste0(Species, "_", Group))
  
  
  #------------------------------------------------------------
  # 3. Apply encoding
  #------------------------------------------------------------
  vsum <- as.data.table(vsum)
  # vsum[, code := encode_veg_sum(MeanCov, Constancy),
  #      by = .(Species, SiteUnit)]
  vsum[ , code := encode_veg_sum(MeanCov, Constancy), by = .(Species2, SiteUnit)]
  #vsum <- vsum %>% merge(spp, by.x = 'Species', by.y = 'Code') 
  # |>    mutate(SiteUnit = str_replace(SiteUnit, "101", "109"))
  
  # nPlot row
  nPlots <- unique(vsum[, .(SiteUnit, nplots)])
  # nPlotRow <- c("", "n Plots", nPlots$nplots, "") |>
  #   matrix(nrow = 1) |>
  #   data.frame() |>
  #   stats::setNames(names(vsum)) |>
  #   data.table::as.data.table()

  #------------------------------------------------------------
  # 5. Cast to wide format
  #------------------------------------------------------------
  cast_table <- function(df) {
    out <- data.table::dcast(
      df,
      Group + ScientificName + EnglishName ~ SiteUnit,
      value.var = "code",
      fill = ""
    )
    out[order(Group, ScientificName), ]
  }
  
  vsum <- cast_table(vsum)
  
  #------------------------------------------------------------
  # 6. Finalize table (ordering, renaming, nPlot row, BGC cleanup)
  #------------------------------------------------------------
  finalize_table <- function(vsum, spp, type) {
    
    # Replace "-remove"
    vsum <- vsum %>% mutate_all(str_replace_all, "-remove", "")
    
    # Optional ordering
    if (exists("indic.order")) {
      vsum <- vsum[order(match(vsum$ScientificName, indic.order$ScientificName)), ]
    } else {
      message("No vegetation ordering applied (indic.order not found).")
    }
    
    # Add lifeform for Auto logic
    lifeform <- spp %>% select(ScientificName, Lifeform)
    
    vsum2 <- vsum %>%
      left_join(lifeform, by = "ScientificName") %>%
      relocate(EnglishName, .after = last_col()) %>%
      select(-Lifeform)
    
    # BGC cleanup
    # if (type == "BGC") {
    #   colnames(vsum2) <- gsub("109", "101", colnames(vsum2))
    #   colnames(vsum2) <- gsub(paste0(bgc.choose, "?"), "", colnames(vsum2))
    #   colnames(vsum2) <- gsub("_", "", colnames(vsum2))
    # }
    
    # wide table: vsum_wide
    # long table: nplots_long (columns: Unit, nplots)
    
    # 1. Create an empty row with the correct number of columns
    plot_row <- as.list(rep("", ncol(vsum2)))
    names(plot_row) <- colnames(vsum2)
    
    # 2. Fill in only the unit columns
    unit_cols <- intersect(colnames(vsum2), nPlots$SiteUnit)
    
    plot_row[unit_cols] <- nPlots$nplots[match(unit_cols, nPlots$SiteUnit)]
    
    # 3. Convert to data.frame
    plot_row <- as.data.frame(plot_row, check.names = FALSE)
    
    # # Save the correct names
    # col_names <- names(vsum2)
    # 
    # # Remove names from both objects
    # names(plot_row) <- NULL
    # names(vsum2)    <- NULL
    
    # Bind by position
    vsum <- rbind(plot_row, vsum2)
    
    # Restore names
    # names(vsum_with_plots) <- col_names
  }
  
  vsum <- finalize_table(vsum, spp, type)
  
  return(vsum)
}
