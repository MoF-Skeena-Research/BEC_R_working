##table graphic that shows plot count by edaphic position for each site series 
### modify this script to group site series by BGC orderly layout.
## see this link for some ideas https://stackoverflow.com/questions/65835639/arrange-gt-tables-side-by-side-or-in-a-grid-or-table-of-tables
# plot.env; su=su2; bgc.choose = "ICHmk2"; su.choose = "ICHmk2_101"

draw_ss_edatope <- function(plot.env, su, bgc.choose) {
  
  su2 <- su %>%
    dplyr::filter(bgc == bgc.choose)
  
  edatopic <- plot.env %>%
    dplyr::select(PlotNumber, rSMR = MoistureRegime, SNR = NutrientRegime) %>%
    dplyr::distinct()
  
  su.edatopic <- su2 %>%
    dplyr::left_join(edatopic) %>%
    dplyr::mutate(
      rSMR = gsub("[+-]", "", as.character(rSMR)),
      SNR  = gsub("[+-]", "", as.character(SNR))
    )
  
  bgc_working <- su.edatopic %>%
    dplyr::filter(!is.na(rSMR), !is.na(SNR)) %>%
    dplyr::select(-PlotNumber) %>%
    dplyr::mutate(edatope = paste0(rSMR, SNR)) %>%
    dplyr::group_by(SiteUnit, edatope) %>%
    dplyr::mutate(edatopic = length(edatope)) %>%
    dplyr::distinct() %>%
    dplyr::ungroup()
  
  su.unique <- unique(bgc_working$SiteUnit)
  
  good_snr <- c("A","B","C","D","E")
  good_smr <- c(0:8)
  
  su_tables <- list()
  
  for (su.choose in su.unique) {
    
    edatopic_default <- expand.grid(SNR = good_snr, rSMR = good_smr) %>%
      as.data.frame() %>%
      dplyr::mutate(
        SiteUnit = su.choose,
        bgc = bgc.choose,
        edatope = paste0(rSMR, SNR),
        edatopic = 0
      ) %>%
      dplyr::mutate(rSMR = as.character(rSMR))
    
    su_working <- bgc_working %>%
      dplyr::filter(SiteUnit == su.choose) %>%
      dplyr::filter(SNR %in% good_snr, rSMR %in% good_smr) %>%
      dplyr::mutate(edatopic = as.integer(edatopic)) %>%
      dplyr::select(SiteUnit, bgc, SNR, rSMR, edatope, edatopic) %>%
      rbind(edatopic_default) %>%
      dplyr::group_by(edatope) %>%
      dplyr::slice_max(edatopic, n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::select(rSMR, SNR, edatopic) %>%
      tidyr::spread(SNR, edatopic) %>%
      dplyr::arrange(rSMR) %>%
      tibble::column_to_rownames("rSMR")
    
    su_working[su_working == 0] <- ""
    
    
    p <- ggplot(
      su_working %>%
        tibble::rownames_to_column("rSMR") %>%
        tidyr::pivot_longer(-rSMR, names_to = "SNR", values_to = "value"),
      aes(x = SNR, y = rSMR, fill = as.numeric(value))
    ) +
      geom_tile(color = "black") +
      geom_text(aes(label = ifelse(value == "", "", value)), size = 2.5) +
      scale_fill_gradient(low = "white", high = "steelblue", na.value = "white") +
      scale_x_discrete(drop = FALSE, limits = c("A","B","C","D","E")) +
      scale_y_discrete(drop = FALSE, limits = rev(as.character(0:8))) +
      coord_fixed() +
      labs(title = su.choose) +
      theme_minimal(base_size = 8) +
      theme(
        legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_blank(),
        plot.margin = margin(2, 2, 2, 2),
        axis.text = element_text(size = 6)
      )
    
    su_tables[[length(su_tables) + 1]] <- p
    

  }
  
  # layout
  if (length(su_tables) == 0) return(NULL)
  library(patchwork)
  
  n <- length(su_tables)
  ncol <- 4
  nrow <- ceiling(n/ncol)
  
  wrap_plots(su_tables, ncol = ncol) +
    plot_layout(
      widths = rep(1, ncol),
      heights = rep(1, nrow)
    ) +
    plot_annotation(title = paste("BGC:", bgc.choose))
  
}

# 
# library(gridExtra)
# draw_ss_edatope <- function(plot.env, su, bgc.choose) {
#   su2 <- su %>%
#     filter(bgc == bgc.choose)
#   bgc.unique <- unique(su2$bgc)
# 
#   edatopic <- plot.env %>%
#     select(PlotNumber, rSMR = MoistureRegime, SNR = NutrientRegime) %>%
#     distinct()
# 
#   su.edatopic <- su2 %>%
#     left_join(edatopic)# %>%
#     #mutate(across(where(~ is.character(.x) || is.factor(.x)),
#                  # ~ gsub("[+-]", "", as.character(.x))))
# 
#   su.edatopic <- su.edatopic %>%
#     mutate(
#       rSMR = gsub("[+-]", "", as.character(rSMR)),
#       SNR  = gsub("[+-]", "", as.character(SNR))
#     )
# 
# for (bgc.choose in bgc.unique) {
#   
#   bgc_working <- su.edatopic %>%
#     filter(bgc == bgc.choose) %>%
#     filter(!is.na(rSMR) & !is.na(SNR)) %>%
#     select(-PlotNumber) %>%
#     mutate(edatope = paste0(rSMR, SNR)) %>%
#     group_by(SiteUnit, edatope) %>%
#     mutate(edatopic = length(edatope)) %>%
#     distinct() %>%
#     ungroup()
#   
#   su.unique <- unique(bgc_working$SiteUnit)
#   
#   good_snr <- c("A", "B", "C", "D", "E")
#   good_smr <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)
#   
#   # ✅ collect tables here
#   su_tables <- list()
#   
#   for (su.choose in su.unique) {
#     
#     edatopic_default <- expand.grid(SNR = good_snr, rSMR = good_smr) %>%
#       as.data.frame() %>%
#       mutate(
#         SiteUnit = su.choose,
#         bgc = bgc.choose,
#         edatope = paste0(rSMR, SNR),
#         edatopic = 0
#       ) %>%
#       arrange(SNR, rSMR) %>%
#       mutate(rSMR = as.character(rSMR)) %>%
#       select(SiteUnit, bgc, SNR, rSMR, edatope, edatopic)
#     
#     su_working <- bgc_working %>%
#       filter(SiteUnit == su.choose) %>%
#       filter(SNR %in% good_snr, rSMR %in% good_smr) %>%
#       as.data.frame() %>%
#       mutate(edatopic = as.integer(edatopic)) %>%
#       arrange(SNR) %>%
#       select(SiteUnit, bgc, SNR, rSMR, edatope, edatopic) %>%
#       rbind(edatopic_default) %>%
#       group_by(edatope) %>%
#       slice_max(edatopic, n = 1) %>%
#       ungroup() %>%
#       select(rSMR, SNR, edatopic) %>%
#       tidyr::spread(SNR, edatopic) %>%
#       arrange(rSMR) %>%
#       column_to_rownames("rSMR") %>%
#       as.data.frame()
#     
#     su_working[su_working == 0] <- ""
#     
#     su_gt <- gt::as_gtable(
#       gt::gt(su_working, rownames_to_stub = TRUE) %>%
#         gt::fmt_number(decimals = 0) |>
#         gt::tab_options(
#           table.font.size = 10,
#           table_body.hlines.color = "gray25",
#           table_body.hlines.width = 1,
#           table_body.vlines.color = "gray25",
#           table_body.vlines.width = 1
#         ) |>
#         gt::tab_header(title = paste0(su.choose)) |>
#         gt::cols_width(everything() ~ gt::px(30)) %>%
#         gt::tab_style(
#           style = gt::cell_borders(
#             sides = "all",
#             color = "#000000",
#             style = "solid",
#             weight = gt::px(1)
#           ),
#           locations = gt::cells_body()
#         ),
#       plot = TRUE
#     )
#     
#     # ✅ store it
#     su_tables[[length(su_tables) + 1]] <- su_gt
#   }
#   
#   n <- length(su_tables)
#   ncol <- ceiling(sqrt(n * (8.5 / 11)))
#   nrow <- ceiling(n / ncol)
#   
#   # ✅ create combined grob (DO NOT draw yet)
#   combined_plot <- arrangeGrob(
#     grobs = su_tables,
#     ncol = ncol,
#     nrow = nrow,
#     top = textGrob(
#       paste("BGC:", bgc.choose),
#       gp = gpar(fontsize = 14, fontface = "bold")
#     )
#   )
#   
#   # ✅ print → Quarto renders it inline
# combined_plot
#   cat("\\newpage")
#   # ✅ combine & draw all SU tables for this BGC
#   # grid::grid.newpage()
#   # gridExtra::grid.arrange(
#   #   grobs = su_tables,
#   #   ncol = 3,                      # adjust layout
#   #   top = paste("BGC:", bgc.choose)
#   # )
# }
# }


# draw_ss_edatope <- function(plot.env, su, bgc.choose) {
#   su2 <- su %>%
#     filter(bgc == bgc) 
#   bgc.unique <- unique(su2$bgc)
#   
#   edatopic <- plot.env %>%
#     select(PlotNumber, rSMR = MoistureRegime, SNR = NutrientRegime) %>%
#     distinct()
# 
#   su.edatopic <- su2 %>%
#     left_join(edatopic)# %>%
#     #mutate(across(where(~ is.character(.x) || is.factor(.x)),
#                  # ~ gsub("[+-]", "", as.character(.x))))
#   
#   su.edatopic <- su.edatopic %>%
#     mutate(
#       rSMR = gsub("[+-]", "", as.character(rSMR)),
#       SNR  = gsub("[+-]", "", as.character(SNR))
#     )
#   
#   
#   for (bgc.choose in bgc.unique) {
#     bgc_working <- su.edatopic %>%
#       filter(bgc == bgc.choose) %>%
#       #dplyr::rename(rSMR = MoistureRegime, SNR = NutrientRegime) %>%
#       #mutate(rSMR = gsub(remove, '', rSMR), SNR = gsub(remove, '', SNR)) %>%
#       filter(!is.na(rSMR) & !is.na(SNR)) %>%
#       select(-PlotNumber) %>%
#       mutate(edatope = paste0(rSMR, SNR)) %>%
#       group_by(SiteUnit, edatope) %>%
#       mutate(edatopic = length(edatope)) %>%
#       distinct() %>%
#       ungroup()
#     su.unique <- unique(bgc_working$SiteUnit)
#     good_snr <- c("A", "B", "C", "D", "E")
#     good_smr <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)
#     
#     for (su.choose in su.unique) {
#       edatopic_default <- expand.grid(SNR = good_snr, rSMR = good_smr) %>%
#         as.data.frame() %>%
#         mutate(SiteUnit = su.choose, bgc = bgc.choose, edatope = paste0(rSMR, SNR), edatopic = 0) %>%
#         arrange(SNR, rSMR) %>% mutate(rSMR=as.character(rSMR)) %>% 
#         select(SiteUnit, bgc, SNR, rSMR, edatope, edatopic)
#       
#       su_working <- bgc_working %>%
#         filter(SiteUnit == su.choose) %>%
#         filter(SNR %in% good_snr, rSMR %in% good_smr) %>%
#         as.data.frame() %>% mutate(edatopic = as.integer(edatopic)) %>% 
#         arrange(SNR) %>%
#         select(SiteUnit, bgc, SNR, rSMR, edatope, edatopic) %>%
#         rbind(edatopic_default) %>%
#         group_by(edatope) %>%
#         slice_max(edatopic, n = 1) %>%
#         ungroup() %>%
#         select(rSMR, SNR, edatopic) %>%
#         spread(SNR, edatopic) %>%
#         arrange(rSMR) %>%
#         column_to_rownames("rSMR") %>%
#         as.data.frame()
#       su_working[su_working == 0] <- ""
#       su_gt <-  gt::as_gtable(gt::gt(su_working, rownames_to_stub = TRUE) %>%
#                                 gt::fmt_number(decimals = 0) |>
#                                 gt::tab_options(table.font.size = 10, table_body.hlines.color = "gray25", table_body.hlines.width = 1, table_body.vlines.color = "gray25", table_body.vlines.width = 1) |>
#                                 gt::tab_header(title = paste0(su.choose)) |>
#                                 gt::cols_width(everything() ~ px(30)) |>
#                                 gt::tab_options() %>%
#                                 gt::tab_style(style = gt::cell_borders(sides = "all", color = "#000000", style = "solid", weight = gt::px(1)), locations = gt::cells_body()), plot = TRUE)# add grid lines to gtable
#       # add grid lines to gtable
#     }
#     
#   }
# }
