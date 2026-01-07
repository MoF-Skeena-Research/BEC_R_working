### Converts Vpro veg table into long form analysis set

create_veg_analysis_dataset <- function(plot.veg, taxon.lifeform){
  plot.veg <- as.data.table(plot.veg)
  fields = c("PlotNumber", "Species", "TotalA", "TotalB", "Cover6", "Cover7")
  vegdat <- plot.veg %>% dplyr::select(all_of(fields)) 
  vegdat[, Cover := rowSums(.SD, na.rm = TRUE), .SDcols = 3:6]
  vegdat <- vegdat %>% dplyr::select(PlotNumber, Species,Cover)
  ### add in lifeform
  vegdat[setDT(taxon.lifeform), "Lifeform" := Lifeform, on = c("Species" = "Code")]
  vegdat <- vegdat %>% filter(Cover > 0 )
  tree_seedlings <- taxon.lifeform %>% filter(Lifeform %in% c("1", "2")) %>% mutate(Code = paste0(Code, "D")) %>% pull(Code)
  trees <- taxon.lifeform %>% filter(Lifeform %in% c("1", "2")) %>% pull(Code)
  vegdat <- vegdat %>% filter(!Species %in% tree_seedlings)
  return(vegdat)
}
  
create_veg_reports_dataset <- function(plot.veg, taxon.lifeform){
  plot.veg <- as.data.table(plot.veg)
  fields = c("PlotNumber", "Species", "TotalA", "TotalB", "Cover6", "Cover7")
  vegdat <- plot.veg %>% dplyr::select(fields) 
  vegdat[setDT(taxon.lifeform), "Lifeform" := Lifeform, on = c("Species" = "Code")]  
  vegdat <- vegdat %>% pivot_longer(cols = c("TotalA", "TotalB", "Cover6", "Cover7"), 
                                    names_to = "Layer", values_to = "Cover")
  
  tree_seedlings <- taxon.lifeform %>% filter(Lifeform %in% c("1", "2")) %>% mutate(Code = paste0(Code, "D")) %>% pull(Code)
  trees <- taxon.lifeform %>% filter(Lifeform %in% c("1", "2")) %>% pull(Code)
  vegdat <- vegdat  %>% filter(!Species %in% tree_seedlings)
  vegdat <- vegdat  %>% filter(!(Species %in% trees & Layer == "Moss"))
  # change value"Total A" in Layer to "Tree" using dplyr
  setDT(vegdat)[, Layer := fcase(
    Layer == "TotalA", "Tree",
    Layer == "TotalB", "Shrub",
    Layer == "Cover6", "Herb",
    Layer == "Cover7", "Moss",
    default = "Layer"
  )]
  # vegdat2 <- vegdat %>% mutate(Layer = ifelse(Layer == "TotalA", "Tree",
  #                                            ifelse(Layer == "TotalB", "Shrub",
  #                                                   ifelse(Layer == "Cover6", "Herb",
  #                                                          ifelse(Layer == "Cover7", "Moss", Layer)))))
  

  vegdat <- vegdat %>% filter(Cover > 0 )
  return(vegdat)
}
