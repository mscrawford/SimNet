# dissmodelr

<!-- badges: start -->
<!-- badges: end -->

The package provides the succulent simulation model from Reineking, Veste, Wissel, Huth (2006) Environmental variability and allocation trade-offs maintain species diversity in a process-based model of succulent plant communities. Ecological Modelling 199: 486-504.

It further provides code to run the simulations used in Crawford et al. (submitted) The function-dominance correlation drives the direction and strength of biodiversity-ecosystem functioning relationships.

## Installation

You can install dissmodelr on Unix-alikes from a local file with:

``` r
install.packages("dissmodelr.tar.gz", repos = NULL, type = "source")
```

## Simulation runs for simnet
``` r
library("dissmodelr")
run_simnet(4)
# run analysis_simnet_dryland.R, you will have to adjust the path to the simulation output
```
