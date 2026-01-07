##convert cover + constancy into importance value for analysis
## apply to table built from VegdatSUsummary function
#vegsum = veg_anal.tree; minimportance = 1; minconstancy = 60; noiseconstancy = 40; minplots = 0; covadj = .75
create_analysis_vegsum <- function(vegsum, minimportance, minconstancy, noiseconstancy, minplots, covadj){
  vegsum <- as.data.frame(vegsum)
  vegsum$MeanCov[vegsum$MeanCov >100] <- 100
  vegsum$spp_importance <- vegsum$MeanCov^covadj
  vegsum$spp_importance[vegsum$spp_importance < minimportance] <- NA 
  #vegsum$spp_importance[vegsum$spp_importance < 1.1] <- 0.1 
  #vegsum$spp_importance[is.na(vegsum$spp_importance)] <- 0 
  vegsum <- vegsum %>% dplyr::mutate(spp_importance = spp_importance * (Constancy/100)) %>% dplyr::filter(spp_importance > minimportance)
  ### remove species where the maximum constancy for species is < contancy
  speciesmax <- vegsum  %>% dplyr::group_by(Species) %>% 
    dplyr::summarise(maxcons = max(Constancy)) %>% 
    dplyr::filter(maxcons >= minconstancy)
  vegsum <- vegsum %>%  dplyr::filter(Species %in% speciesmax$Species) %>%
    dplyr::filter(nplots >= minplots) %>% filter(Constancy > noiseconstancy)
  return(vegsum)
}
