## sets characters for the vegetation summary
encode_veg_sum <- function(coverage, constancy) {
  black <- 'n'
  grey <- 'l'
  star <- 'v'
  char <- black
  if (constancy < 70) {
    char <- grey
  } 
  if (constancy <= 50) {
    char <- star
  }
  color = "remove"
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
}
