#' Read succulent model plant parameter file
#'
#' @param file File containing plant parameter values
#' @return Data frame with plant parameter values
#' @export
read_plantpars <- function(file) {
  x <- readLines(file)
  x <- gsub(" = ", " ", x)
  x <- strsplit(x, " ")
  # sapply(x, function(z) z)
  m <- lapply(x, matrix, ncol = 2, byrow = TRUE)
  m <- sapply(m, function(x) {
    v <- as.numeric(x[,2])
    names(v) <- x[,1]
    v
  })
  data.frame(t(m))
}
