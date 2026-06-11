find_nearest_neighbours <- function(unit.compare, single.units) {
  library(data.table)
  
  dt <- as.data.table(unit.compare)[, .(Unit1, Unit2, BEC.sim.min, BEC.sim.max)]
  dt <- unique(dt)
  
  su.units <- unique(su2$SiteUnit)
  single.units2 <- single.units %>% filter(SiteUnit %in% su.units) %>% distinct %>% pull()
  
  results <- list()
  
  for (SiteUnit in single.units2) {
    
    su <- as.character(SiteUnit)
    
    filtered <- dt[(get("Unit1") == su | get("Unit2") == su)]
    filtered <- filtered[complete.cases(filtered)]
    
    if (nrow(filtered) == 0) next
    
    # 🔥 Normalize so Unit1 is always the SiteUnit
    filtered <- filtered[, {
      swap <- Unit2 == su
      list(
        Unit1 = su,
        Unit2 = ifelse(swap, Unit1, Unit2),
        BEC.sim.min = BEC.sim.min,
        BEC.sim.max = BEC.sim.max
      )
    }]
    filtered <- filtered %>% rename(Singleton = Unit1, Neighbor = Unit2)
    
    take_n <- min(2, nrow(filtered))
    
    nearest <- filtered[
      order(BEC.sim.min, decreasing = TRUE)
    ][1:take_n]
    
    results[[su]] <- nearest
    
  }
  
  final_results <- rbindlist(results, use.names = TRUE, fill = TRUE)
  return(final_results)
}