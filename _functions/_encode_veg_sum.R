## sets characters for the vegetation summary
encode_veg_sum <- function(coverage, constancy) {
  
  black <- 'n'
  grey  <- 'v'
  star  <- 'l'
  
  # Vectorized assignment of char
  char <- data.table::fifelse(
    constancy <= 50, star,
    data.table::fifelse(
      constancy < 70, grey, black
    )
  )
  
  color <- "remove"
  
  code <- data.table::fcase(
    coverage <=   1, paste0(char, "-", color),
    coverage <=   3, paste0(char, char, "-", color),
    coverage <=  10, paste0(strrep(char, 3), "-", color),
    coverage <=  25, paste0(strrep(char, 4), "-", color),
    coverage <= 100, paste0(strrep(char, 5), "-", color),
    default = strrep(char, 6)
  )
  
  return(code)
}


# encode_veg_sum <- function(coverage, constancy) {
#   black <- 'n'
#   grey <- 'v'
#   star <- 'l'
#   char <- black
#   if (constancy < 70) {
#     char <- grey
#   } 
#   if (constancy <= 50) {
#     char <- star
#   }
#   color = "remove"
#   code <- data.table::fcase(
#     coverage <=   1,
#     sprintf('%s-%s', paste0(rep(char, 1), collapse=''), color),
#     coverage <=   3,
#     sprintf('%s-%s', paste0(rep(char, 2), collapse=''), color),
#     coverage <=  10,
#     sprintf('%s-%s', paste0(rep(char, 3), collapse=''), color),
#     coverage <=  25,
#     sprintf('%s-%s', paste0(rep(char, 4), collapse=''), color),
#     coverage <= 100,
#     sprintf('%s-%s', paste0(rep(char, 5), collapse=''), color),
#     default = paste0(rep(char, 6), collapse='')
#   )
# }

