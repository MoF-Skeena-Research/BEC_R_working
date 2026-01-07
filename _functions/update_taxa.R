

# sppmaster = taxon.all
# vegdata = vtab
update_taxa <- function(vegdata, sppmaster){
  lookup <- sppmaster[sppmaster$OldCode != sppmaster$Code,c("Code", "OldCode")] #%>% dplyr::rename("Species" = OldCode) ##creates code update crosswalk
  setDT(vegdata)[setDT(lookup), "Species" := Code, on = c("Species" = "OldCode")]
  return(vegdata)
}


