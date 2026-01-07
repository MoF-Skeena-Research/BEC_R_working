
# db = "E:/Arctic/Vpro64_Data_Arctic/AVACanada2024.accdb"
# proj = "AVACanada2024"
require(DBI)
import_vpro <- function(db, proj){
  becmaster <- dbConnect(odbc::odbc(), .connection_string = paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)}; DBQ=", db, ";"))
  plot.env <- dbReadTable(becmaster, paste0(proj,"_Env")) %>% 
    mutate(Longitude = ifelse(Longitude<0, Longitude, 0-Longitude))
  plot.admin <- dbReadTable(becmaster, paste0(proj,"_Admin"))
  plot.veg <- dbReadTable(becmaster, paste0(proj,"_Veg"))
  plot.mineral <- dbReadTable(becmaster, paste0(proj,"_Mineral"))
    plot.humus <- dbReadTable(becmaster, paste0(proj,"_Humus"))
  plot.meta <- dbReadTable(becmaster, paste0(proj,"_Metadata"))
  dbDisconnect(becmaster)
     return(list(env = plot.env, veg = plot.veg, info = plot.admin, soil = plot.mineral, humus = plot.humus, metadata = plot.meta))
 }

