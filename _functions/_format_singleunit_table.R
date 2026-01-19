##formats a vegetation summary table
##types of summary tables: BGC =  site series within a BGC; zonal = comparison of zonal vegtation between BGCs
#type = "Zonal"
# vsum = vegSum; spp=taxon.lifeform; min.cons = 50 ; strata.by = "Layer"; type = "numeric"
# file_name = "AA_D_upland_units2"; file_path = "./vegsum.tables/"

format_singleunit_table <- function(
    vsum = vegSum,
    spp = taxon.all,
    type = "BGC",
    min.cons = 50,
    strata.by = c("Layer", "Lifeform"),
    output.type = c("numeric", "graphic", "both"),
    file_path = "./vegsum.tables/",
    file_name = "working_unit_tables"
){
  
  strata.by   <- match.arg(strata.by)
  output.type <- match.arg(output.type)
  
  # Add lifeform for Auto logic
  taxa.info <- spp %>% select(ScientificName, EnglishName, Code)
  #------------------------------------------------------------
  # 0. Lifeform lookup table (editable, explicit)
  #------------------------------------------------------------
  lifeform_layer <- c(
    "1" = "Tree",
    "2" = "Tree",
    "3" = "Shrub",
    "4" = "Shrub",
    "5" = "Herb",
    "6" = "Herb",
    "7" = "Herb",
    "8" = "Herb",
    "12" = "Herb", 
    "9" = "Moss",
    "10" = "Moss",
    "11" = "Moss"
  )
  # Convert to joinable table
  layer_tbl <- tibble::tibble(
    Lifeform = names(lifeform_layer),
    LayerName = unname(lifeform_layer)
  )
  
  
  lifeform_names <- c(
    "1" = "Conifer Tree",
    "2" = "Deciduous Tree",
    "3" = "Evergreen Shrub",
    "4" = "Deciduous Shrub",
    "12" = "Dwarf woody plant", 
    "5" = "Fern",
    "6" = "Graminoid",
    "7" = "Forb",
    "8" = "Parasitic",
    "9" = "Moss",
    "10" = "Liverwort",
    "11" = "Lichen",
    "13" = "Macro algae"
  )
  
  # Convert to joinable table
  lifeform_tbl <- tibble::tibble(
    Lifeform = names(lifeform_names),
    LifeformName = unname(lifeform_names)
  )
  #------------------------------------------------------------
  # 1. Prepare long table
  #------------------------------------------------------------
  
  vsum2 <- vsum %>%
    mutate(
      Constancy    = round(Constancy, 1),
      MeanCov      = round(MeanCov, 1),
      LifeformName = lifeform_names[as.character(Lifeform)],
      Layer        = lifeform_layer[as.character(Lifeform)],   # <-- NEW
      cons.cov     = paste0(Constancy, "-", MeanCov)
    ) %>%
    left_join(taxa.info, by = c("Species" = "Code")) %>%
    tidyr::drop_na()
  
  
  
  dups <- vsum2 %>%
    count(LifeformName, ScientificName, EnglishName, SiteUnit) %>%
    filter(n > 1)
  
  dups
  
  #------------------------------------------------------------
  # 2. Cast to wide format
  #------------------------------------------------------------
  cast_table <- function(df) {
    
    if (strata.by == "Lifeform") {
      out <- data.table::dcast(
        df,
        Lifeform + LifeformName + ScientificName + EnglishName ~ SiteUnit,
        value.var = "cons.cov",
        fill = ""
      )
      out <- out[order(Lifeform, ScientificName), ]
    }
    
    if (strata.by == "Layer") {
      out <- data.table::dcast(
        df,
        Layer + ScientificName + EnglishName ~ SiteUnit,
        value.var = "cons.cov",
        fill = ""
      )
      out <- out[order(Layer, ScientificName), ]
    }
    
    out
  }
  
  vsum2 <- cast_table(vsum2)
  
  #------------------------------------------------------------
  # 3. Identify unit columns correctly
  #------------------------------------------------------------
  if (strata.by == "Lifeform") {
    unit_cols <- setdiff(colnames(vsum2),
                         c("Lifeform", "LifeformName", "ScientificName", "EnglishName"))
  }
  
  if (strata.by == "Layer") {
    unit_cols <- setdiff(colnames(vsum2),
                         c("Layer", "ScientificName", "EnglishName"))
  }
  
  #------------------------------------------------------------
  # 4. Split cons.cov
  #------------------------------------------------------------
  split_cons_cov <- function(x) {
    parts <- strsplit(x, "-", fixed = TRUE)
    const <- as.numeric(sapply(parts, `[`, 1))
    mean  <- as.numeric(sapply(parts, `[`, 2))
    tibble(Constancy = const, MeanCover = mean)
  }
  
  #------------------------------------------------------------
  # 5. Build unit tables
  #------------------------------------------------------------
  unit_tables <- lapply(unit_cols, function(unit) {
    
    if (strata.by == "Lifeform") {
      
      df_unit <- vsum2 %>%
        select(Lifeform, LifeformName, ScientificName, EnglishName, all_of(unit)) %>%
        rename(cons.cov = all_of(unit)) %>%
        filter(cons.cov != "") %>%
        bind_cols(split_cons_cov(.$cons.cov)) %>%
        arrange(Lifeform, desc(MeanCover)) %>%
        select(LifeformName, ScientificName, Constancy, MeanCover, EnglishName)
      
      present_levels <- unique(df_unit$LifeformName)
      
      df_unit$LifeformName <- factor(df_unit$LifeformName,
                                     levels = present_levels,
                                     ordered = TRUE)
      
      df_unit <- df_unit %>% arrange(LifeformName, desc(MeanCover))
    }
    
    if (strata.by == "Layer") {
      
      df_unit <- vsum2 %>%
        select(Layer, ScientificName, EnglishName, all_of(unit)) %>%
        rename(cons.cov = all_of(unit)) %>%
        filter(cons.cov != "") %>%
        bind_cols(split_cons_cov(.$cons.cov)) %>%
        arrange(Layer, desc(MeanCover)) %>%
        select(Layer, ScientificName, Constancy, MeanCover, EnglishName)
      
      desired_layers <- c("Tree", "Shrub", "Herb", "Moss")
      present_layers <- intersect(desired_layers, unique(df_unit$Layer))
      
      df_unit$Layer <- factor(df_unit$Layer,
                              levels = present_layers,
                              ordered = TRUE)
      
      df_unit <- df_unit %>% arrange(Layer, desc(MeanCover))
    }
    
    df_unit
  })
  
  names(unit_tables) <- unit_cols
  
  #------------------------------------------------------------
  # 6. Save workbook
  #------------------------------------------------------------
  file_path <- ifelse(grepl("/$", file_path), file_path, paste0(file_path, "/"))
  full_file <- paste0(file_path, file_name, ".xlsx")
  
  openxlsx::write.xlsx(unit_tables, file = full_file)
}

