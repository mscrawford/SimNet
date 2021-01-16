#' Write succulent model plant parameter file
#'
#' @param x Data framw with plant parameters
#' @param file Name of file to create
#' @export
write_plantpars <- function(x, file) {
  text <- apply(x, 1, function(z) paste(colnames(x), z, sep = " = ", collapse = " "))
  writeLines(text, file)
}
