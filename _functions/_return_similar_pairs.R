##convert cover + constancy into importance value for analysis
## apply to table built from VegdatSUsummary function
#vegsum = veg_anal.tree; minimportance = 1; minconstancy = 60; noiseconstancy = 40; minplots = 0; covadj = .75
return_similar_pairs <- function(veg.data, similarity = 1, neighbours = 1){
  su.similar <- vegsum.pairs %>% select(Unit1, Unit2, diag.ratio) %>% dplyr::filter(diag.ratio >similarity) %>% filter(!Unit1 == Unit2)  %>%  distinct
  su.similar <- su.similar[!duplicated(apply(su.similar,1,function(x) paste(sort(x),collapse=''))),] %>% rowid_to_column("Working_ID") %>% mutate(Working_ID = paste0("working-", Working_ID)) 
  su.similar1 <- su.similar %>%  select(Unit1, Working_ID, diag.ratio) %>% rename(SiteUnit = 1)
  su.similar2 <- su.similar %>%  select(Unit2, Working_ID, diag.ratio) %>% rename(SiteUnit = 1)
  su.similar.tot <- rbind(su.similar1, su.similar2) %>% distinct 
  su.similar <- su.similar.tot %>% group_by(SiteUnit) %>% mutate(n = n()) %>% slice_min(diag.ratio, n =neighbours) %>% ungroup %>%
    left_join(su.similar.tot, by = "Working_ID") %>% filter(!SiteUnit.x == SiteUnit.y)
  return(su.similar)
}
