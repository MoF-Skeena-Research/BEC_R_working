combine_species <- function(
    vegdata,
    lumpfile,
    strata.by = c("Auto", "Lifeform", "Layer"),
    use.subtaxa = FALSE
) {
  
  strata.by <- match.arg(strata.by)
  
  # --- Auto-detect grouping variable ---
  if (strata.by == "Auto") {
    if ("Lifeform" %in% names(vegdata)) {
      strata.by <- "Lifeform"
    } else if ("Layer" %in% names(vegdata)) {
      strata.by <- "Layer"
    } else {
      stop("Auto mode could not find either 'Lifeform' or 'Layer' in vegdata.")
    }
  }
  
  # --- Apply lumping only if lumpfile is provided ---
  if (!is.null(lumpfile)) {
    setDT(vegdata)[setDT(lumpfile),
                   "Species" := LumpCode,
                   on = c("Species" = "SppCode")
    ]
  }
  
  # --- Remove subtaxa if requested ---
  if (isFALSE(use.subtaxa)) {
    vegdata[, Species := gsub("[0-9]+", "", Species)]
  }
  
  # --- Dynamic grouping ---
  group_vars <- c("PlotNumber", "Species", strata.by)
  
  vegdata2 <- vegdata[, .(Cover = sum(Cover)), by = group_vars]
  
  # --- Return tidy tibble ---
  vegdata_out <- vegdata2 %>%
    dplyr::select(PlotNumber, Species, Cover, all_of(strata.by)) %>%
    dplyr::arrange(PlotNumber, Species)
  
  return(vegdata_out)
}
