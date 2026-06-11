## summarizes veg dat using su table
create_veg_sum_bgc <- function(vegdata, siteunit.tbl, BGC = bgc.choose, minconstancy = 50, noiseconstancy = 10, minimportance = 0, strata.by = "Layer") {
  if (strata.by == "Layer") {
    vegdata <- combine_taxa_strata(vegdata, lump = lump)
  } else if (strata.by == "Lifeform") {
    vegdata <- vegdata <- lump_species_lifeform(vegdata, lump = lump)
  }

  vegdata <- merge(vegdata, siteunit.tbl, by = "PlotNumber")
  setDT(vegdata)
  vegdata <- vegdata[PlotNumber %in% siteunit.tbl$PlotNumber, ]
  # vegdata <- vegdata %>% filter(bgc %in% BGC)
  vegdata <- vegdata[bgc %in% BGC, ]
## remove trees in moss layer
  # vegdata <-  vegdata  %>% filter(!Species %in% tree_seedlings)
  # vegdata <-  vegdata  %>% filter(!(Species %in% trees & Layer == "Moss"))

#veg.dat2 <- lump_species2(vegdata = vegdata, lump, use.subtaxa = FALSE)
  
  vegdata <- vegdata[, if (.N > 0) .SD, by = .(SiteUnit, Species)]
  vegdata[, nplots := length(unique(PlotNumber)), by = .(SiteUnit)]
  if (strata.by == "Layer") {
    vegdata <- vegdata[, .(
      MeanCov = sum(Cover, na.rm = TRUE) / unique(nplots), # should this just be mean, is NA assumed to be 0?
      Constancy = (.N / unique(nplots)) * 100,
      nplots = unique(nplots),
      importance = (sum(Cover, na.rm = TRUE) / unique(nplots))^(1/2) * (.N / unique(nplots))
    ), by = .(SiteUnit, Species, Layer)]
  } else if (strata.by == "Lifeform") {
    vegdata <- vegdata[, .(
      MeanCov = sum(Cover, na.rm = TRUE) / unique(nplots), # should this just be mean, is NA assumed to be 0?
      Constancy = (.N / unique(nplots)) * 100,
      nplots = unique(nplots),
      importance = (sum(Cover, na.rm = TRUE) / unique(nplots))^(1/2) * (.N / unique(nplots))
    ), by = .(SiteUnit, Species, Lifeform)]
  }
  vegdata[, maxcons := max(Constancy), by = .(Species)]
  vegdata[, maximportance := max(importance), by = .(Species)]
  vegdata <- vegdata[maxcons >= minconstancy, ]
  vegdata <- vegdata[Constancy >= noiseconstancy, ]
  vegdata <- vegdata[importance >= minimportance, ]
}

