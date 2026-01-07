
build_species_ordering_hier <- function(vdat, vsum = vegSum, code.lump = lump, siteUnits, hier.level = "Working") {
  vegdata <- lump_species(vdat, code.lump)
  #hier.grp <- siteUnits %>% dplyr::select({{hier.level}}) %>% distinct()
  #cluster.choose <- unique(hier.grp[[hier.level]])
  #su.choose <- siteUnits %>% filter(hier.grp[[hier.level]] %in% cluster.choose)
  
  filter_site_units <- function(siteUnits, hier.level) {
    # Capture the column symbol
    col_sym <- ensym(hier.level)
    
    # Get distinct group values from siteUnits
    cluster.choose <- siteUnits %>%
      distinct(!!col_sym) %>%
      pull(!!col_sym) %>%
      unique()
    
    # Filter rows based on cluster.choose
    su.choose <- siteUnits %>%
      filter((!!col_sym) %in% cluster.choose)
    
    return(su.choose)
  }
  
  # Example usage:
  su.choose <- filter_site_units(siteUnits, hier.level = "Working")
  
  
  setDT(vegdata)
  vegdata <- merge(su.choose, vegdata, by = "PlotNumber")
  
  SS <- vegdata %>%
    select(PlotNumber, SiteUnit) %>%
    distinct() %>%
    pull(SiteUnit)
  
  veg_anal <- vegdata %>%
    filter(Species %in% vsum$Species) %>%
    select(PlotNumber, Species, Cover) %>%
    group_by(PlotNumber, Species) %>%
    summarise(Cover = sum(Cover, na.rm = TRUE)) %>%
    ungroup() %>%
    data.frame()
  
  veg_anal <- matrify(veg_anal)
  
  n_units <- length(unique(SS))
  
  if (n_units < 2) {
    message("Only one site unit in group. Returning species sorted by total cover.")
    spp_summary <- veg_anal %>%
      as.data.frame() %>%
      rownames_to_column("spp") %>%
      pivot_longer(-spp, names_to = "PlotNumber", values_to = "Cover") %>%
      group_by(spp) %>%
      summarise(indic = sum(Cover, na.rm = TRUE)) %>%
      left_join(taxon.lifeform, by = c("spp" = "Code")) %>%
      arrange(desc(indic))
    return(spp_summary)
  }
  
  library(indicspecies)
  library(permute)
  library(doParallel)
  
  # Set number of cores
  num_cores <- parallel::detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Set up permutation control using permute's `how()` function
  ctrl <- how(nperm = 9)  # More perms now that we’re speeding things up
  
  # Run multipatt with parallel backend
  indval <- multipatt(veg_anal, SS, control = ctrl)
  
  # Stop the cluster when done
  stopCluster(cl)
  
  #indval <- indicspecies::multipatt(veg_anal, SS, control = how(nperm = 9))
  
  indic.order <- indval$str %>%
    data.frame() %>%
    select(1:(all_of(n_units))) %>%
    rownames_to_column("spp") %>%
    pivot_longer(-spp, names_to = "siteunit", values_to = "indic") %>%
    group_by(spp) %>%
    mutate(max_indic = max(indic)) %>%
    filter(indic == max_indic) %>%
    ungroup() %>%
    left_join(taxon.lifeform, by = c("spp" = "Code")) %>%
    arrange(desc(indic))
  
  return(indic.order)
}
