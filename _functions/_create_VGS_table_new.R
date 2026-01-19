##df = veg.sum.table; title = "test";  table.name = "test"; save_path = "./figures/"; file_type = "png"; table_type = "all"

create_VGS_table <- function(
    veg.sum.table,
    table.name = "",
    save_path = NULL,
    file_type = c("html", "png", "pdf", "rtf", "docx", "tex"),
    #table_type = c("all", "by_unit"),
    font_family = "Helvetica"   # <-- NEW PARAMETER
){
  
  file_type  <- match.arg(file_type)
  #table_type <- match.arg(table_type)
  
  # Helper to build a gt table
  build_gt <- function(df, title, unit_col = NULL) {
    
    # Determine wingding columns
    #if (is.null(unit_col)) {
      wing_cols <- 4:(ncol(df) - 1)
    #} else {
    #  wing_cols <- unit_col
    #}
    
    gt::gt(df) |>
      gt::opt_table_font(font = gt::google_font(font_family)) |>   # <-- NEW
      gt::tab_options(table.font.size = 8) |>
      gt::tab_style(
        style = gt::cell_text(font = gt::google_font("wingdings")),
        locations = gt::cells_body(
          columns = wing_cols,
          rows = 2:nrow(df)
        )
      ) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_body(rows = 1)
      ) |>
      gt::tab_style(                                              # <-- NEW
        style = gt::cell_text(style = "italic"),
        locations = gt::cells_body(columns = ScientificName)
      ) |>
      gt::tab_header(
        title = title,
        subtitle = "Summary Vegetation Table"
      )
  }
  
  
  # Helper to save a gt table
  save_gt <- function(gt_obj, name) {
    
    if (is.null(save_path)) return(invisible(NULL))
    
    # If save_path is a directory, append filename
    if (dir.exists(save_path) || !grepl("\\.[A-Za-z0-9]+$", save_path)) {
      save_dir <- save_path
      if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
      out_path <- file.path(save_dir, paste0(name, ".", file_type))
    } else {
      out_path <- save_path
      if (!grepl(paste0("\\.", file_type, "$"), out_path, ignore.case = TRUE)) {
        out_path <- paste0(out_path, ".", file_type)
      }
      dir_path <- dirname(out_path)
      if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
    }
    
    gt::gtsave(gt_obj, filename = out_path)
  }
  
  # -------------------------------------------------------------------
  # CASE 1: table_type = "all"  → current behavior
  # -------------------------------------------------------------------
  #if (table_type == "all") {
    
    gt_obj <- build_gt(veg.sum.table, title = table.name, unit_col = unit)
    
    save_gt(gt_obj, table.name)
    
    return(gt_obj)
  }
#}
  # -------------------------------------------------------------------
  # CASE 2: table_type = "by_unit"
  # -------------------------------------------------------------------
#   tables <- list()
#   
#   # identify useful column names
#   all_cols    <- colnames(veg.sum.table)
#   layer_col   <- all_cols[1]
#   sci_col     <- all_cols[2]
#   english_col <- tail(all_cols, 1)
#   
#   # adjust these if your ns/ls/vs columns are named differently
#   ns_col <- "ns"
#   ls_col <- "ls"
#   vs_col <- "vs"
#   
#   for (unit in unit_cols) {
#     
#     # Start from full table so we can use ns/ls/vs for sorting
#     sub_df_full <- veg.sum.table
#     
#     # --- 1. Split plot row from species rows ---
#     plot_row   <- sub_df_full[1, , drop = FALSE]
#     species_df <- sub_df_full[-1, , drop = FALSE]
#     
#     # --- 2. Keep only species with non-empty values in this unit ---
#     species_df <- species_df[
#       species_df[[unit]] != "" & !is.na(species_df[[unit]]),
#       , drop = FALSE]
#     
#     # If no species remain, skip this unit
#     if (nrow(species_df) == 0) next
#     
#     # --- 3. Order species by Layer, then ns, ls, vs ---
#     desired_levels <- c("Tree", "Shrub", "Herb", "Moss")
#     
#     # keep only levels that actually appear in the data
#     present_levels <- intersect(desired_levels, unique(species_df[[layer_col]]))
#     
#     species_df[[layer_col]] <- factor(
#       species_df[[layer_col]],
#       levels = present_levels,
#       ordered = TRUE
#     )
#     
#     
#     species_df <- species_df[
#       order(
#         species_df[[layer_col]],
#         -species_df[[n_col]],
#         -species_df[[l_col]],
#         -species_df[[v_col]]
#       ),
#       , drop = FALSE]
#     
#     # --- 4. Build the display table: plot row + ordered species, with only required columns ---
#     
#     # plot row label: "nplots = <number>"
#     nplots_value <- plot_row[[unit]]
#     plot_row[[sci_col]]   <- paste0("nplots = ", nplots_value)
#     plot_row[[unit]]      <- ""  # don't show the number again in the unit column
#     plot_row[[english_col]] <- ""  # usually blank for the plot row, adjust if needed
#     
#     # now restrict to the 4 desired columns
#     plot_row_disp <- plot_row[, c(layer_col, sci_col, unit, english_col), drop = FALSE]
#     species_disp  <- species_df[, c(layer_col, sci_col, unit, english_col), drop = FALSE]
#     
#     # combine back together
#     sub_df <- rbind(plot_row_disp, species_disp)
#     
#     # --- 5. Build and save gt table ---
#     gt_obj <- build_gt(sub_df, title = unit, unit_col = unit)
#     
#     save_gt(gt_obj, name = unit)
#     
#     tables[[unit]] <- gt_obj
#   }
#   
#   return(tables)
# }
