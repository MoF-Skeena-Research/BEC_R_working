create_veg_sum <- function(
    vegdata,
    siteunit.tbl,
    siteunit.var = "SiteUnit",
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
      stop("Auto mode could not find either 'Layer' or 'Lifeform' in vegdata.")
    }
  }
  
  # --- Apply species lumping ---
  if (strata.by %in% c("Layer", "Lifeform")) {
    vegdata <- combine_taxa(vegdata, lumpfile = lumpfile)
  }
  
  # --- Merge site units using user-specified table and variable ---
  setDT(vegdata)
  setDT(siteunit.tbl)
  
  # Ensure siteunit.tbl has PlotNumber
  if (!"PlotNumber" %in% names(siteunit.tbl)) {
    stop("siteunit.tbl must contain a 'PlotNumber' column.")
  }
  
  # Ensure siteunit.var exists
  if (!siteunit.var %in% names(siteunit.tbl)) {
    stop(paste0("Column '", siteunit.var, "' not found in siteunit.tbl."))
  }
  
  vdat <- merge(vegdata, siteunit.tbl, by = "PlotNumber")
  
  # Keep only plots that appear in the siteunit table
  vdat <- vdat[PlotNumber %in% siteunit.tbl$PlotNumber]
  
  # --- Keep only species present in each site unit ---
  vdat <- vdat[, if (.N >= 1) .SD, by = c(siteunit.var, "Species")]
  
  # --- Count plots per site unit ---
  vdat[, nplots := uniqueN(PlotNumber), by = siteunit.var]
  
  # --- Dynamic grouping variable ---
  group_vars <- c(siteunit.var, "Species", strata.by)
  
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
