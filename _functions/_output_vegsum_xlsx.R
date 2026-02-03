##formats a vegetation summary table and exports to excel for Unit Reporting functions
## will want to combine all exports into one function eventually
##types of summary tables: BGC =  site series within a BGC; zonal = comparison of zonal vegtation between BGCs
#type = "Zonal"
# vsum = vegSum; spp=taxon.lifeform; min.cons = 50 ; strata.by = "Layer"; type = "numeric"
# file_name = "AA_D_upland_units2"; file_path = "./vegsum.tables/"
output_vegsum_xlsx <- function(
    vsum = vegSum,
    spp = taxon.all,
    type = "BGC",
    min.cons = 50,
    strata.by = c("Layer", "Lifeform"),
    output.type = c("numeric", "graphic", "both"),
    file_path = "./vegsum.tables/",
    file_name = "working_unit_tables",
    separate_files = FALSE
){
  
  strata.by   <- match.arg(strata.by)
  output.type <- match.arg(output.type)
  
  # Detect the unit variable from the first column
  unit_var <- names(vsum)[1]
  
  taxa.info <- spp %>% select(ScientificName, EnglishName, Code)
  
  lifeform_layer <- c(
    "1" = "Tree","2" = "Tree",
    "3" = "Shrub","4" = "Shrub",
    "5" = "Herb","6" = "Herb","7" = "Herb","8" = "Herb","12" = "Herb",
    "9" = "Moss","10" = "Moss","11" = "Moss"
  )
  
  lifeform_names <- c(
    "1" = "Conifer Tree","2" = "Deciduous Tree",
    "3" = "Evergreen Shrub","4" = "Deciduous Shrub",
    "12" = "Dwarf woody plant",
    "5" = "Fern","6" = "Graminoid","7" = "Forb","8" = "Parasitic",
    "9" = "Moss","10" = "Liverwort","11" = "Lichen","13" = "Macro algae"
  )
  lifeform_order <- unname(lifeform_names)
  
  vsum2 <- vsum %>%
    mutate(
      Constancy    = round(Constancy, 1),
      MeanCov      = round(MeanCov, 1),
      LifeformName = lifeform_names[as.character(Lifeform)],
      Layer        = lifeform_layer[as.character(Lifeform)],
      cons.cov     = paste0(Constancy, "-", MeanCov)
    ) %>%
    left_join(taxa.info, by = c("Species" = "Code")) %>%
    tidyr::drop_na()
  
  # Duplicate detection using dynamic unit variable
  dups <- vsum2 %>%
    count(LifeformName, ScientificName, EnglishName, !!sym(unit_var)) %>%
    filter(n > 1)
  
  # Cast to wide format
  cast_table <- function(df) {
    
    if (strata.by == "Lifeform") {
      out <- data.table::dcast(
        df,
        as.formula(
          paste("Lifeform + LifeformName + ScientificName + EnglishName ~", unit_var)
        ),
        value.var = "cons.cov",
        fill = ""
      )
      out <- out[order(Lifeform, ScientificName), ]
    }
    
    if (strata.by == "Layer") {
      out <- data.table::dcast(
        df,
        as.formula(
          paste("Layer + ScientificName + EnglishName ~", unit_var)
        ),
        value.var = "cons.cov",
        fill = ""
      )
      out <- out[order(Layer, ScientificName), ]
    }
    
    out
  }
  
  vsum2 <- cast_table(vsum2)
  
  # Identify unit columns
  if (strata.by == "Lifeform") {
    unit_cols <- setdiff(colnames(vsum2),
                         c("Lifeform", "LifeformName", "ScientificName", "EnglishName"))
  }
  
  if (strata.by == "Layer") {
    unit_cols <- setdiff(colnames(vsum2),
                         c("Layer", "ScientificName", "EnglishName"))
  }
  
  split_cons_cov <- function(x) {
    parts <- strsplit(x, "-", fixed = TRUE)
    tibble(
      Constancy = as.numeric(sapply(parts, `[`, 1)),
      MeanCover = as.numeric(sapply(parts, `[`, 2))
    )
  }
  
  unit_tables <- lapply(unit_cols, function(unit) {
    
    if (strata.by == "Lifeform") {
      df_unit <- vsum2 %>%
        select(Lifeform, LifeformName, ScientificName, EnglishName, all_of(unit)) %>%
        rename(cons.cov = all_of(unit)) %>%
        filter(cons.cov != "") %>%
        bind_cols(split_cons_cov(.$cons.cov)) %>%
        filter(Constancy >= min.cons) %>%
        arrange(Lifeform, desc(MeanCover)) %>%
        select(LifeformName, ScientificName, Constancy, MeanCover, EnglishName)
      
      df_unit$LifeformName <- factor(df_unit$LifeformName,
                                     levels = lifeform_order,
                                     ordered = TRUE)
      
      df_unit <- df_unit %>% arrange(LifeformName, desc(MeanCover))
    }
    
    if (strata.by == "Layer") {
      df_unit <- vsum2 %>%
        select(Layer, ScientificName, EnglishName, all_of(unit)) %>%
        rename(cons.cov = all_of(unit)) %>%
        filter(cons.cov != "") %>%
        bind_cols(split_cons_cov(.$cons.cov)) %>%
        filter(Constancy >= min.cons) %>%
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
  
  
  write_scientific_txt <- function(scientific_names, out_file) {
    
    # Ensure directory exists
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    
    # Collapse names into a simple text block
    txt <- paste(scientific_names, collapse = "\n")
    
    # Write the text file
    writeLines(txt, con = out_file)
  }
  
  

  # write_scientific_docx <- function(scientific_names, out_file, use_style = TRUE) {
  #   
  #   # Ensure directory exists
  #   dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  #   
  #    doc <- read_docx("BEC_Concept_Template.docx")
  #   
  #   for (nm in scientific_names) {
  #     
  #     if (use_style) {
  #       # Scientific must be a PARAGRAPH style
  #       doc <- body_add_par(
  #         doc,
  #         value = nm,
  #         style = "Scientific"
  #       )
  #       
  #     } else {
  #       doc <- body_add_par(
  #         doc,
  #         value = ftext(nm, prop = fp_text(italic = TRUE)),
  #         style = "Normal"
  #       )
  #     }
  #   }
  #   dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  #   
  #   print(doc, target = out_file)
  # }
  # 
  
  # Save workbook(s)
  file_path <- ifelse(grepl("/$", file_path), file_path, paste0(file_path, "/"))
  
  if (!separate_files) {
    full_file <- paste0(file_path, file_name, ".xlsx")
    openxlsx::write.xlsx(unit_tables, file = full_file)
  } else {
    for (unit in names(unit_tables)) {
      
      # --- Excel output ---
      subdir <- file.path(file_path, unit, "vegsum_raw")
      dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
      out_file <- file.path(subdir, paste0(unit, "_vegsum.xlsx"))
      openxlsx::write.xlsx(list(unit_tables[[unit]]), file = out_file)
      
      # --- Word output ---
      concept_dir <- file.path(file_path, unit, "concept_description")
      dir.create(concept_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Extract scientific names for this unit
      sci_names <- unit_tables[[unit]]$ScientificName |> unique()
      
      # Output file
      txt_file <- file.path(concept_dir, paste0(unit, "_concept.txt"))
      
      write_scientific_txt(
        scientific_names = sci_names,
        out_file = txt_file
      )
      
    }
    
  }
}
   
