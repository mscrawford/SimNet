#' Make simulation pools
#'
#' @param simulations Data frame with simulation parameters
#' @return Named list with species pools
#' @export
make_species_pools <- function(simulations) {
  xx <- simulations[simulations$SeedRain == 0 & simulations$Wdh == 1, c("Ninitial", "Rep")]
  pool <- apply(xx, 1, function(x) sort(sample.int(64, x[["Ninitial"]])))

  # Need to check that there are no duplicates...
  pool_size <- sapply(pool, length)
  pool_size_1 <- which(pool_size == 1)
  for (i in seq(pool_size_1)) {
    pool[[pool_size_1[i]]] <- i
  }


  for (i in unique(pool_size)) {
      if(any(duplicated(unname(sapply(pool[pool_size == i], function(x) x))))) {
        stop("Duplicated species in pool size", i)
        }
  }

  names(pool) <- paste0(formatC(xx$Ninitial, width = 2, flag = "0"),
                                "_",
                        formatC(xx$Rep, width = 2, flag = "0"))
  pool
}
