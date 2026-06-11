##formats a vegetation summary table
##types of summary tables: BGC =  site series within a BGC; zonal = comparison of zonal vegtation between BGCs
#type = "Zonal"
#vsum = vegSum; spp=taxon.lifeform; cons.1 = 70; cons.2 = 50; strata.by = "Lifeform"; type = NULL; indic.order = indic.order

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
  
  strata_order <- c("Tree", "Regen", "Shrub", "Herb", "Moss")

  
  #------------------------------------------------------------
  # 1. Encode coverage + constancy into symbol codes
  #------------------------------------------------------------
 ## use sourced funtion
   # encode_veg_sum <- function(coverage, constancy) {
   #  black <- "n"; grey <- "v"; star <- "l"
   #  char <- black
   #  if (constancy < cons.1) char <- grey
   #  if (constancy < cons.2) char <- star
   #  
   #  color <- "remove"
   #  
   #  data.table::fcase(
   #    coverage <   1, sprintf("%s-%s", strrep(char, 1), color),
   #    coverage <   3, sprintf("%s-%s", strrep(char, 2), color),
   #    coverage <  10, sprintf("%s-%s", strrep(char, 3), color),
   #    coverage <  25, sprintf("%s-%s", strrep(char, 4), color),
   #    coverage < 100, sprintf("%s-%s", strrep(char, 5), color),
   #    default = strrep(char, 6)
   #  )
   # }
  #------------------------------------------------------------
  # 1b. Grouping logic (Layer / Lifeform / Auto)
  #------------------------------------------------------------
  order_groups <- function(df, strata.by) {
    
    if (strata.by == "Layer") {
      # Layer-based ordering
      df$Group <- factor(df$Group, levels = strata_order)
      df <- df[order(df$Group), ]
      
    } else if (strata.by == "Lifeform") {
      # LifeformName-based ordering
      df$Group <- factor(df$Group, levels = lifeform_names)
      df <- df[order(df$Group), ]
      
    } else if (strata.by == "Auto") {
      # Auto uses Layer logic (Tree → Regen → Shrub → Herb → Moss)
      df$Group <- factor(df$Group, levels = strata_order)
      df <- df[order(df$Group), ]
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
    vsum <- vsum %>%
      merge(spp %>% select(ScientificName, Lifeform),
            by = "ScientificName",
            all.x = TRUE)
  }
  
  # Derive Layer from Lifeform if missing
  if (has_lifeform && !has_layer) {
    message("Layer not found in input; deriving from Lifeform mapping.")
    vsum <- vsum %>%
      mutate(Layer = lifeform_map[as.character(Lifeform)])
  }
  
  # Add LifeformName
  vsum <- vsum %>%
    mutate(LifeformName = lifeform_names[as.character(Lifeform)])
  
  if (any(is.na(vsum$LifeformName))) {
    warning("Some species have missing or unmapped Lifeform values. See saved missing_lifeform.csv for details.")
    missing_lf <- vsum %>%
      filter(is.na(LifeformName)) %>%
      select(Species, ScientificName, Lifeform) %>%
      distinct()
    fwrite(missing_lf, "missing_lifeform.csv")
    
  }
  
  # Apply strata.by BEFORE Species2 is created
  #vsum <- group_species(vsum, strata.by)
  
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
  vsum2 = vsum
  
  #------------------------------------------------------------
  # 6. Finalize table (ordering, renaming, nPlot row, BGC cleanup)
  #------------------------------------------------------------
  
  finalize_table <- function(vsum, spp, type, strata.by = NULL) {
    
    # Replace "-remove"
    vsum <- vsum %>% mutate_all(str_replace_all, "-remove", "")
    
    # ------------------------------------------------------------
    # CASE 1: strata.by == "Layer" → skip LifeformName logic
    # ------------------------------------------------------------
    if (!is.null(strata.by) && strata.by == "Layer") {
      
      strata_order <- c("Tree", "Regen", "Shrub", "Herb", "Moss")
      
      # Order by strata
      if ("Layer" %in% names(vsum)) {
        vsum <- vsum %>%
          mutate(Layer = factor(Layer, levels = strata_order)) %>%
          arrange(Layer)
      }
      
      # No LifeformName → skip lifeform_map and skip spp join
      vsum2 <- vsum
      
    } else {
      
      # ------------------------------------------------------------
      # CASE 2: Normal logic (LifeformName → Lifeform)
      # ------------------------------------------------------------
      # Optional species ordering
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
      
    }
    
    # ------------------------------------------------------------
    # Add nPlots row
    # ------------------------------------------------------------
    
    plot_row <- as.list(rep("", ncol(vsum2)))
    names(plot_row) <- colnames(vsum2)
    
    unit_cols <- intersect(colnames(vsum2), nPlots$SiteUnit)
    plot_row[unit_cols] <- nPlots$nplots[match(unit_cols, nPlots$SiteUnit)]
    
    plot_row <- as.data.frame(plot_row, check.names = FALSE)
    
    # Bind
    vsum <- rbind(plot_row, vsum2)
        return(vsum)
  }
  
  vsum.final <- finalize_table(vsum, spp, type)
  
  return(vsum.final)
}
