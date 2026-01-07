process_unit_compare <- function(unit_subset) {
  unit.compare2 <- unit_subset
  
  if (include.few.plots) {
    unit.few <- unit.compare2 %>%
      filter(nplots.x < few.plots.threshold) %>%
      select(Unit1, nplots.x) %>%
      distinct()
    
    unit.compare2 <- unit.compare2 %>%
      mutate(Unit1 = ifelse(Unit1 %in% unit.few$Unit1, paste0(Unit1, "-nplot"), Unit1),
             Unit2 = ifelse(Unit2 %in% unit.few$Unit1, paste0(Unit2, "-nplot"), Unit2))
  }
  
  if (include.low.diag) {
    unit.simple <- unit.compare2 %>%
      filter(unit.diag.sum.x < low.diag.threshold) %>%
      select(Unit1, unit.diag.sum.x) %>%
      distinct()
    
    unit.compare2 <- unit.compare2 %>%
      mutate(Unit1 = ifelse(Unit1 %in% unit.simple$Unit1, paste0(Unit1, "-lowdiag"), Unit1),
             Unit2 = ifelse(Unit2 %in% unit.simple$Unit1, paste0(Unit2, "-lowdiag"), Unit2))
  }
  
  unit.compare2 <- unit.compare2 %>%
    select(Unit1, Unit2, BEC.sim.min)
  
  return(unit.compare2)
}