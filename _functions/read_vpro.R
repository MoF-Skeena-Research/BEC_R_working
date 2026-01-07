## require(Hmisc)
read_vpro <- function(database = "D:/BECMaster64/BECMaster64.accdb", project = "BECMaster64"){
  becmaster <- dbConnect(odbc::odbc(), .connection_string = paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)}; DBQ=", database,";"))
   plot.env <- dbReadTable(becmaster, paste0(project,"_ENV")) %>% mutate(Longitude = ifelse(Longitude<0, Longitude, 0-Longitude))
  plot.veg <- dbReadTable(becmaster, paste0(project,"_Veg"))  
  plot.info <- dbReadTable(becmaster, paste0(project,"_Admin"))
  plot.min   <- dbReadTable(becmaster, paste0(project,"_Mineral"))  
  plot.ff   <- dbReadTable(becmaster, paste0(project,"_Humus")) 
  dbDisconnect(becmaster)
  return(list(plot.env, plot.veg, plot.info, plot.min, plot.ff))
}

