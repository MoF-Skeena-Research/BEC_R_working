###########################################################################
#                                                                         #   
# Data Manipulation                                                       #
#                                                                         #
###########################################################################
make_veg_sum_code <- function(coverage, constancy) {
  fullBlock <- r'(\blacksquare)' #'▊' #'&#x2588;'
  star <- r'(\textasteriskcentered)'
  if (constancy <= 60) {
    char <- star
    spec <- r'(\textcolor{gray}{%s})'
    wrap <- ''
  } else if (constancy < 80) {
    char <- fullBlock
    spec <- r'(\textcolor{gray}{%s})'
    wrap <- '$'
  } else {
    char <- fullBlock
    spec <- r'(\textcolor{black}{%s})'
    wrap <- '$'
  }
  code <- data.table::fcase(
    coverage <   1,
    sprintf('%s%s%s', wrap, paste0(rep(char, 1), collapse=''), wrap),
    coverage <   3,
    sprintf('%s%s%s', wrap, paste0(rep(char, 2), collapse=''), wrap),
    coverage <  10,
    sprintf('%s%s%s', wrap, paste0(rep(char, 3), collapse=''), wrap),
    coverage <  25,
    sprintf('%s%s%s', wrap, paste0(rep(char, 4), collapse=''), wrap),
    coverage < 100,
    sprintf('%s%s%s', wrap, paste0(rep(char, 5), collapse=''), wrap),
    default = paste0(rep(char, 6), collapse='')
  )
  sprintf(spec, code)
}
encode_veg_sum <- function(coverage, constancy) {
  fullBlock <- '▊'
  star <- '*'
  color <- 'black'
  char <- fullBlock
  if (constancy <= 60) {
    char <- star
  }
  if (constancy < 80) {
    color <- 'gray'
  }
  code <- data.table::fcase(
    coverage <   1, 
    sprintf('%s-%s', paste0(rep(char, 1), collapse=''), color),
    coverage <   3, 
    sprintf('%s-%s', paste0(rep(char, 2), collapse=''), color),
    coverage <  10, 
    sprintf('%s-%s', paste0(rep(char, 3), collapse=''), color),
    coverage <  25, 
    sprintf('%s-%s', paste0(rep(char, 4), collapse=''), color),
    coverage < 100, 
    sprintf('%s-%s', paste0(rep(char, 5), collapse=''), color),
    default = paste0(rep(char, 6), collapse='')
  )
  code
}
convert_veg_sum_code <- function(vegSumCode) {
  fullBlock <- r'(\blacksquare)'
  star <- r'(\textasteriskcentered)'
  vegSumCode |> 
    gsub(pattern = '([▊\\*]+)-(gray|black)', replacement = '\\\\textcolor{\\2}{\\1}') |> 
    gsub(pattern = '*', replacement = star, fixed = TRUE) |> 
    sub(pattern = '{▊', replacement = '{$▊', fixed = TRUE) |> 
    sub(pattern = '▊}', replacement = '▊$}', fixed = TRUE) |> 
    gsub(pattern = '▊', replacement = fullBlock, fixed = TRUE)
}
summarize_vegdata <- function(vegData, siteUnits) {
  vegData <- merge(vegData, siteUnits, by = 'PlotNumber')
  vegData <- vegData[,if(.N > 1) .SD, by = .(SiteUnit,Species)]
  vegData[ , nplots := length(unique(PlotNumber)), by = .(SiteUnit)]
  vegData[,.(
    MeanCov = sum(Cover, na.rm = TRUE) / unique(nplots), # should this just be mean, is NA assumed to be 0?
    Constancy = (.N / unique(nplots)) * 100, 
    nplots = unique(nplots)
  ), by = .(SiteUnit, Species, Lifeform)]
}
lump_species <- function(vegData, lump, useSubtaxa = FALSE) {
  lumpSpecies <- merge(vegData, lump, by.x = 'Species', by.y = 'SppCode', all.x = TRUE)
  if (isFALSE(useSubtaxa)){
    lumpSpecies[ , Species := gsub('[0-9]+', '', Species)]
  }
  lumpSpecies[ , .(Cover = sum(Cover)), by = .(PlotNumber, Species, Lifeform)]
}
format_dir_name <- function(siteUnit) {
  gsub('\\s|\\/', '', siteUnit)
}
###########################################################################
#                                                                         #   
# Figures/Tables                                                          #
#                                                                         #
###########################################################################
convert_docx_to_tex <- function(siteUnitDir, fileName = 'unit-description.tex') {
  filePath <- list.files(siteUnitDir, full.names = TRUE, 
    pattern = '(\\.doc|\\.docx)$')
  if (length(filePath) == 0) {
    stop('No .doc or .docx file found')
  } else if (length(filePath) > 1) {
    stop('More than one .doc or docx file found.')
  } else {
    rmarkdown::pandoc_convert(filePath, to = 'latex', wd = getwd(),
      output = sprintf('%s/%s', siteUnitDir, fileName))
  }
}
copy_docx_template <- function(siteUnit, fileName = 'unit-description.docx') {
  filePath <- list.files(sprintf('reports/%s', format_dir_name(siteUnit)), 
    pattern = '(\\.doc|\\.docx)$')
  if (length(filePath) > 0) {
    message(sprintf('Word file "%s" already exists.', filePath[1]))
  } else {
    filePath <- sprintf('reports/%s/%s', format_dir_name(siteUnit), fileName)
    file.copy(fileName, filePath)
  }
}
save_map_plot <- function(siteUnit, latLon, mapBg, fileName = 'unit-range-map.pdf', 
  overwrite = FALSE) {
  filePath <- sprintf('reports/%s/%s', format_dir_name(siteUnit), fileName)
  if (file.exists(filePath) && isFALSE(overwrite)) {
    message(sprintf('"%s" already exists and overwrite is set to FALSE.', 
      filePath))
    return(NULL)
  }
  g <- ggplot2::ggplot(mapBg) +
    ggplot2::geom_sf(fill='transparent', color = 'black') +
    ggplot2::geom_sf(data = latLon[grep(pattern = siteUnit, latLon$SiteUnit), ], 
      color='red') +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5, 
      family = 'serif'))
  ggplot2::ggsave(filename = filePath, plot = g, width = 3, height = 3)
}
save_veg_sum <- function(siteUnit, vegSum, fileName = 'veg-sum.xlsx', 
  overwrite = FALSE) {
  filePath <- sprintf('reports/%s/%s', format_dir_name(siteUnit), fileName)
  if (file.exists(filePath) && isFALSE(overwrite)) {
    message(sprintf('"%s" already exists and overwrite is set to FALSE.', 
      filePath))
    return(NULL)
  }
  vegSum <- vegSum[grepl(pattern = siteUnit, x = SiteUnit) & Constancy >= 40, ]
  vegSum[ , code := encode_veg_sum(MeanCov, Constancy), by = .(Species, SiteUnit)]
  nPlots <- unique(vegSum[ ,.(SiteUnit, nplots)])[order(nplots, decreasing = TRUE), ]
  vegSum <- data.table::dcast(vegSum, 
    Lifeform + ReportName + ScientificName + EnglishName ~ SiteUnit, 
    value.var = 'code',
    fill = '')[order(Lifeform, ReportName, ScientificName), ]
  vegSum[duplicated(ReportName, fromLast = TRUE), ReportName := '']
  vegSum <- vegSum[ , c('ReportName', 'ScientificName', nPlots$SiteUnit, 'EnglishName'), with = FALSE]
  data.table::setnames(vegSum, old = c('ReportName', 'ScientificName', 
    'EnglishName'), new = c('Lifeform', 'Scientific', 'Common Name'))
  nPlotRow <- c('', 'n Plots', nPlots$nplots, '') |> 
    matrix(nrow = 1) |>
    data.frame() |> 
    stats::setNames(names(vegSum)) |> 
    data.table::as.data.table()
  vegSum <- rbind(nPlotRow, vegSum)
  writexl::write_xlsx(x = vegSum, path = filePath)
}
###########################################################################
#                                                                         #   
# Output                                                                  #
#                                                                         #
###########################################################################
prepare_quarto_materials <- function(siteUnit, latLon, mapBg, vegSum) {
  dirName <- sprintf('reports/%s', format_dir_name(siteUnit))
  if (!dir.exists(dirName)) {
    dir.create(dirName)
  }
  copy_docx_template(siteUnit = siteUnit)
  save_map_plot(siteUnit = siteUnit, latLon = latLon, mapBg = mapBg, overwrite = TRUE)
  save_veg_sum(siteUnit = siteUnit, vegSum = vegSum, overwrite = TRUE)
  print(sprintf('%s processed', siteUnit))
}
build_quarto_report <- function(siteUnitDir, associations) {
  tempQmd <- sprintf('reports/%s/bec-report-template.qmd', siteUnitDir)
  on.exit({file.remove(tempQmd)})
  file.copy('bec-report-template.qmd', tempQmd)
  convert_docx_to_tex(siteUnitDir = sprintf('reports/%s', siteUnitDir))
  tryCatch({
  quarto::quarto_render(
    input = tempQmd, 
    output_file = sprintf('%s.pdf', siteUnitDir),
    execute_params = list(siteUnitDir = siteUnitDir,
       title01 = paste(associations[SiteUnits == siteUnitDir, Class], "Association"),
       title02 = associations[SiteUnits == siteUnitDir, SiteSeriesLongName],
       title03 = associations[SiteUnits == siteUnitDir, SiteSeriesScientificName])
  )
  }, error = function(e) print(e$message))
  file.remove(sprintf('reports/%s/%s.pdf', siteUnitDir, siteUnitDir))
  file.rename(sprintf('%s.pdf', siteUnitDir), sprintf('reports/%s/%s.pdf', siteUnitDir, siteUnitDir))
}
