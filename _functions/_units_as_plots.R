## converts meancover + constancy values from a veg summary table into cover values for comparison to plots

units_as_plots <- function(vegsum, mincover = 0){
  vegsum <- as.data.frame(vegsum)
  vegsum$MeanCov[vegsum$MeanCov >100] <- 100
  vegsum$cover <- vegsum$MeanCov * (vegsum$Constancy *.01)
  vegsum <- vegsum %>% dplyr::filter(cover >= mincover)
  ### remove species where the maximum constancy for species is < contancy
  return(vegsum)
}
