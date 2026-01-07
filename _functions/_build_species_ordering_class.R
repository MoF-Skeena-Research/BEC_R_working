## order species in table using multipatt from the indicspecies package
build_species_ordering <- function(vdat, vsum = vegSum,  siteUnits, BGC) { #code.lump = NULL,
  # if (!is.null(code.lump)) {
  #   vegdata <- lump_species(vdat, code.lump)
  #   
  # } else {
  vegdata <- vdat#}
  #vegdata <- lump_species2(vdat, code.lump)
  vegdata <- merge(vegdata, siteUnits, by = "PlotNumber")
  vegdata <- vegdata[PlotNumber %in% siteUnits$PlotNumber, ]
  # vegData <- vegData %>% filter(bgc %in% BGC)
  vegdata <- vegdata[bgc %in% BGC, ] %>% data.frame()
  SS <- vegdata %>%
    select(PlotNumber, SiteUnit) %>%
    distinct()
  ss <- SS$SiteUnit
  veg_anal <- vegdata %>%
    filter(Species %in% vsum$Species) %>%
    select(PlotNumber, Species, Cover) %>%
    group_by(PlotNumber, Species) %>%
    summarise(Cover = sum(Cover, na.rm = TRUE)) %>%
    ungroup() %>%
    data.frame()
   veg_anal <- labdsv::matrify(veg_anal)
   ## capture conditions where no species in a plot are in the selected vsum species
   # veg.anal.plots <- rownames(veg_anal)%>% data.frame() %>% rename( PlotNumber= 1)
   # veg.anal.plots$present <- "exists"
   # present <- left_join(SS, vv) %>% filter(is.na(present)) %>% select(PlotNumber)
   # SS <- SS %>% filter(PlotNumber %in% present$PlotNumber)
   # ss <- SS$SiteUnit
   n_units <- length(unique(vsum$SiteUnit))

  indval <- indicspecies::multipatt(veg_anal, ss,
    control = how(nperm = 9)
  )
  # summary(indval, alpha=1)
  indic.order <- indval$str %>%
    data.frame() %>%
    select(1:(all_of(n_units))) %>%
    rownames_to_column("spp") %>%
    pivot_longer(-spp, names_to = "siteunit", values_to = "indic") %>%
    group_by(spp) %>%
    mutate(max_indic = max(indic)) %>%
    filter(indic == max_indic) %>%
    ungroup() %>%
    left_join(taxon.lifeform, by = c(spp = "Code")) %>%
    arrange(desc(indic)) %>%
    mutate(siteunit = str_replace(siteunit, "101", "109")) %>%
    arrange(siteunit)
}


build_species_ordering_all <- function(vdat, vsum = vegSum, code.lump = NULL, siteUnits) {
  if (!is.null(code.lump)) {
    vegdata <- lump_species(vdat, code.lump)
    
  } else {
    vegdata <- vdat}
  setDT(vegdata)
  vegdata <- merge(vegdata, siteUnits, by = "PlotNumber")
  vegdata <- vegdata[PlotNumber %in% siteUnits$PlotNumber, ]
  # vegData <- vegData %>% filter(bgc %in% BGC)
  #vegdata <- vegdata[clusters %in% cluster.choose, ] %>% data.frame()
  SS <- vegdata %>%
    select(PlotNumber, SiteUnit) %>%
    distinct()
  SS <- SS$SiteUnit
  veg_anal <- vegdata %>%
    filter(Species %in% vsum$Species) %>%
    select(PlotNumber, Species, Cover) %>%
    group_by(PlotNumber, Species) %>%
    summarise(Cover = sum(Cover, na.rm = TRUE)) %>%
    ungroup() %>%
    data.frame()
  veg_anal <- matrify(veg_anal)
  
  
  n_units <- length(unique(vsum$SiteUnit))
  indval <- indicspecies::multipatt(veg_anal, SS,
                                    control = how(nperm = 9)
  )
  # summary(indval, alpha=1)
  indic.order <- indval$str %>%
    data.frame() %>%
    select(1:(all_of(n_units))) %>%
    rownames_to_column("spp") %>%
    pivot_longer(-spp, names_to = "siteunit", values_to = "indic") %>%
    group_by(spp) %>%
    mutate(max_indic = max(indic)) %>%
    filter(indic == max_indic) %>%
    ungroup() %>%
    left_join(taxon.lifeform, by = c(spp = "Code")) %>%
    arrange(desc(indic)) %>% arrange(siteunit)
}
