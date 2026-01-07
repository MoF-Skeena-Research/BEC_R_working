create_veg_sum <- function(
    vegdata,
    siteUnits,
    minconstancy = 60,
    noiseconstancy = 10,
    strata.by = c("Auto","Layer", "Lifeform"),
    minimportance = 0,
    lumpfile = NULL
) {
  
  strata.by <- match.arg(strata.by)
  
  # --- Auto-detect strata variable ---
  if (strata.by == "Auto") {
    if ("Layer" %in% names(vegdata)) {
      strata.by <- "Layer"
    } else if ("Lifeform" %in% names(vegdata)) {
      strata.by <- "Lifeform"
    } else {
      stop("Auto mode could not find either 'Layer' or 'Lifeform' in vdat.")
    }
  }
  
  # --- Apply species lumping depending on strata ---
  if (strata.by == "Layer") {
    vegdata <- combine_species(vegdata, lumpfile = lumpfile)
  } else if (strata.by == "Lifeform") {
    vegdata <- combine_species(vegdata, lumpfile = lumpfile)
  }
  
  # --- Merge site units ---
  setDT(vegdata)
  vdat <- merge(vegdata, siteUnits, by = "PlotNumber")
  vdat <- vdat[PlotNumber %in% siteUnits$PlotNumber]
  
  # --- Keep only species present in each SiteUnit ---
  vdat <- vdat[, if (.N >= 1) .SD, by = .(SiteUnit, Species)]
  
  # --- Count plots per SiteUnit ---
  vdat[, nplots := uniqueN(PlotNumber), by = .(SiteUnit)]
  
  # --- Dynamic grouping variable ---
  group_vars <- c("SiteUnit", "Species", strata.by)
  
  # --- Summaries by strata ---
  vdat <- vdat[, .(
    MeanCov = sum(Cover, na.rm = TRUE) / unique(nplots),
    Constancy = (.N / unique(nplots)) * 100,
    nplots = unique(nplots),
    importance = (sum(Cover, na.rm = TRUE) / unique(nplots))^(1/2) *
      (.N / unique(nplots))
  ), by = group_vars]
  
  # --- Species-level maxima for filtering ---
  vdat[, maxcons := max(Constancy), by = .(Species)]
  vdat[, maximportance := max(importance), by = .(Species)]
  
  # --- Apply filters ---
  vdat <- vdat[maxcons > minconstancy]
  vdat <- vdat[Constancy > noiseconstancy]
  vdat <- vdat[importance > minimportance]
  
  return(vdat)
}
