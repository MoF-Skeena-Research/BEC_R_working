#vegdata = vegdata; lumpfile = lump; use.subtaxa = FALSE
combine_taxa_strata <- function(vegdata, lumpfile, use.subtaxa = FALSE){
  setDT(vegdata)[setDT(lumpfile), "Species" := LumpCode, on = c("Species" = "SppCode")]
  if (isFALSE(use.subtaxa)){
    vegdata$Species <-   gsub('[0-9]+', '', vegdata$Species)
  }
  
  vegdata <- vegdata[
    ,
    .(
      Cover = sum(Cover),
      Lifeform = first(Lifeform)
    ),
    by = .(PlotNumber, Species, Layer)
  ]
  
  vegdata <- vegdata  %>% dplyr::select(PlotNumber,Species, Cover, Layer, Lifeform)
  
  return(vegdata)
}
