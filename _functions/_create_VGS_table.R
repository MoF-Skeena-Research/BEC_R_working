create_VGS_table <- function(veg.sum.table, table.name = ""){
  veg.sum.gt <- gt::gt(veg.sum.table) |> gt::tab_options(table.font.size = 8)|> gt::tab_style(
    style = gt::cell_text(font = gt::google_font(name = "wingdings")),
    locations = gt::cells_body(columns = 4:(ncol(veg.sum.table))-1, rows = 2:nrow(veg.sum.table))) |>  gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(rows =1)) |>
    gt::tab_header(
      title = paste0(table.name),
      subtitle = paste("Summary Vegetation Table")
    )
  return(veg.sum.gt)
}