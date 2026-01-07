##table graphic that shows plot count by edaphic position for each site series 
### modify this script to group site series by BGC orderly layout.
## see this link for some ideas https://stackoverflow.com/questions/65835639/arrange-gt-tables-side-by-side-or-in-a-grid-or-table-of-tables
draw_ss_edatope <- function(plot.env, su, bgc.choose) {
  su2 <- su %>%
    filter(bgc == bgc) 
  bgc.unique <- unique(su2$bgc)
  
  edatopic <- plot.env %>%
    select(PlotNumber, MoistureRegime, NutrientRegime) %>%
    distinct()
  # join edatopic to su
  su.edatopic <- su2 %>%
    left_join(edatopic)
  remove <- c("+|-")
  for (bgc.choose in bgc.unique) {
    bgc_working <- su.edatopic %>%
      filter(bgc == bgc.choose) %>%
      dplyr::rename(rSMR = MoistureRegime, SNR = NutrientRegime) %>%
      mutate(rSMR = gsub(remove, '', rSMR), SNR = gsub(remove, '', SNR)) %>%
      filter(!is.na(rSMR) & !is.na(SNR)) %>%
      select(-PlotNumber) %>%
      mutate(edatope = paste0(rSMR, SNR)) %>%
      group_by(SiteUnit, edatope) %>%
      mutate(edatopic = length(edatope)) %>%
      distinct() %>%
      ungroup()
    su.unique <- unique(bgc_working$SiteUnit)
    good_snr <- c("A", "B", "C", "D", "E")
    good_smr <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)
    
    for (su.choose in su.unique) {
      edatopic_default <- expand.grid(SNR = good_snr, rSMR = good_smr) %>%
        as.data.frame() %>%
        mutate(SiteUnit = su.choose, bgc = bgc.choose, edatope = paste0(rSMR, SNR), edatopic = 0) %>%
        arrange(SNR, rSMR) %>%
        select(SiteUnit, bgc, SNR, rSMR, edatope, edatopic)
      
      su_working <- bgc_working %>%
        filter(SiteUnit == su.choose) %>%
        filter(SNR %in% good_snr, rSMR %in% good_smr) %>%
        as.data.frame() %>%
        arrange(SNR) %>%
        rbind(edatopic_default) %>%
        group_by(edatope) %>%
        slice_max(edatopic, n = 1) %>%
        ungroup() %>%
        select(rSMR, SNR, edatopic) %>%
        spread(SNR, edatopic) %>%
        arrange(rSMR) %>%
        column_to_rownames("rSMR") %>%
        as.data.frame()
      su_working[su_working == 0] <- ""
      su_gt <-  gt::as_gtable(gt::gt(su_working, rownames_to_stub = TRUE) %>%
                                gt::fmt_number(decimals = 0) |>
                                gt::tab_options(table.font.size = 10, table_body.hlines.color = "gray25", table_body.hlines.width = 1, table_body.vlines.color = "gray25", table_body.vlines.width = 1) |>
                                gt::tab_header(title = paste0(su.choose)) |>
                                gt::cols_width(everything() ~ px(30)) |>
                                gt::tab_options() %>%
                                gt::tab_style(style = gt::cell_borders(sides = "all", color = "#000000", style = "solid", weight = gt::px(1)), locations = gt::cells_body()), plot = TRUE)# add grid lines to gtable
      # add grid lines to gtable
    }
    
  }
}
