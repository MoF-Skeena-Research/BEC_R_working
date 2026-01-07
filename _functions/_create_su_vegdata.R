 # su = su
 # vegdata = veg.dat2
 # unit = "SiteUnit"
create_su_vegdata <- function(vegdata, su, constancycut = 0){
  su.choice <- su %>% select(SiteUnit)
  vegdat <- as.data.table(vegdata)
  vegdat[su, SiteUnit := i.SiteUnit, on = "PlotNumber"] ## limit data to those listed in SiteUnit
  vegdat <- vegdat[!is.na(SiteUnit) & SiteUnit != "",]
  vegdat <- unique(vegdat[!is.na(SiteUnit) & SiteUnit != "",])
  vegdat3 <- vegdat[,if(.N > 1) .SD, by = .(SiteUnit,Species)]
  vegdat3[,nplots := length(unique(PlotNumber)), by = .(SiteUnit)]
  vegsum <- vegdat3[,.(MeanCov = sum(Cover, na.rm = TRUE)/nplots[1], Constancy = (.N/nplots[1])*100, nplots = nplots[1]), by = .(SiteUnit,Species)]
  return(vegsum)
}
create_hier_vegdata <- function(vegdata, su, unit = "SiteUnit", constancycut = 0){
  su.choice <- su %>% select(SiteUnit, all_of(unit)) %>% rename(unit = 2)
  vegdat <- as.data.table(vegdata)
  vegdat[su, SiteUnit := i.SiteUnit, on = "PlotNumber"] ## limit data to those listed in SiteUnit
  vegdat <- vegdat[!is.na(SiteUnit) & SiteUnit != "",]
  vegdat <- unique(vegdat[!is.na(SiteUnit) & SiteUnit != "",])
  vegdat3 <- vegdat[,if(.N > 1) .SD, by = .(SiteUnit,Species)]
  vegdat4 <- vegdat3[su.choice, "unit" := unit, on = "SiteUnit"] ## swith SiteUnit to level selected
  vegdat4[,nplots := length(unique(PlotNumber)), by = .(unit)]
  vegsum <- vegdat4[,.(MeanCov = sum(Cover, na.rm = TRUE)/nplots[1], Constancy = (.N/nplots[1])*100, nplots = nplots[1]), by = .(unit,Species)]
  return(vegsum)
}
