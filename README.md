# Welcome to the Granger causality in extremes repository!

This repository contains the code to reproduce the results of the manuscript 'Granger Causality in Extremes' by Juraj Bodik and Olivier C. Pasche, available on ArXiv at <https://arxiv.org/abs/2407.09632>. 
For any questions or error reports, please contact [juraj.bodik@unil.ch](mailto:juraj.bodik@unil.ch).

## Overview

### ExtremeGranger R package

The methods introduced in the paper are implemented as an open-source `ExtremeGranger` R package, available at <https://github.com/opasche/ExtremeGranger>. 
For more details, see its documentation website at <https://opasche.github.io/ExtremeGranger/>. 
To reproduce the results from the manuscript, please install version `0.1.0` of the package, using the following command in R:

```r
# install.packages("devtools")
devtools::install_github("opasche/ExtremeGranger@v0.1.0")
```

The main functions of the package are `Extreme_causality_test` and `Extreme_causality_graph` (and its parallelized version `Extreme_causality_graph_parallel`). 

* The `Extreme_causality_test` function tests whether the tails/extremes of time series `X` cause those of time series `Y`, given confounders `Z`. See [its documentation](https://opasche.github.io/ExtremeGranger/reference/Extreme_causality_test.html) for details on the inputs and outputs.
* The `Extreme_causality_graph` function estimates the causal graph (path diagram) between a set of time series. (See [its documentation](https://opasche.github.io/ExtremeGranger/reference/Extreme_causality_graph.html)).
* The `Extreme_causality_graph_parallel` does the same, but using `doFuture` parallelization. (See [its documentation](https://opasche.github.io/ExtremeGranger/reference/Extreme_causality_graph_parallel.html)).

A local copy of those functions is also available at `./R/Main_functions_localcopy.R`, which can be used without installing the package, as a redundancy backup.  

### Repository contents

The rest of the repository is organized as follows:

* The folder `./main/` contains the code to reproduce the results of the manuscript, in the respective subfolders.
* Results from these scripts are output in the respective `./Results/` subfolders.
* The `./R/` folder contains the R functions used in the `./main/` scripts, as well as a local copy of the package functions, as described above.
* The `./data/Crypto_Application/` folder contains the data used in the Cryptocurrency application. 
The river and meteorological data should is expected in the `./data/Meteo_Application/` folder, but is not included in the repository, due to copyright reasons. 
Please follow instructions in the manuscript to obtain them.


## Example usage of the ExtremeGranger package

Here is a usage example with a 4-dimensional toy VAR time series with `lag = 2`.

```r
# (Install and) load v0.1.0 of the ExtremeGranger R package
# # install.packages("devtools")
# devtools::install_github("opasche/ExtremeGranger@v0.1.0")
library(ExtremeGranger)

# Load EnvStats for the example (or any other package that can generate Pareto noise)
library(EnvStats) 

# Example: Generating a 4-dimensional VAR time series with lag = 2
n <- 5000

set.seed(0) # for reproducibility
epsilon_x <- rpareto(n, 1, 1)
epsilon_y <- rpareto(n, 1, 1)
epsilon_z1 <- rpareto(n, 1, 1)
epsilon_z2 <- rpareto(n, 1, 1)
x <- rep(0, n)
y <- rep(0, n)
z1 <- rep(0, n)
z2 <- rep(0, n)

for (i in 3:n) {
  z1[i] <- 0.5 * z1[i - 1] + epsilon_z1[i]
  z2[i] <- 0.5 * z2[i - 1] + epsilon_z2[i]
  x[i] <- 0.5 * x[i - 1] + 0.5 * z1[i - 2] + 0.5 * z2[i - 1] + epsilon_x[i]
  y[i] <- 0.5 * y[i - 1] + 0.5 * z1[i - 1] + 0.5 * z2[i - 2] + 0.2 * x[i - 1] + epsilon_y[i]
}

# Running the extreme causality tests
z <- data.frame(z1, z2)
result_xy <- Extreme_causality_test(x, y, z, max_causal_lag = 2, p_value_computation = FALSE)
result_yx <- Extreme_causality_test(y, x, z, max_causal_lag = 2, p_value_computation = FALSE)

# Estimating the full causality graph
w <- data.frame(z1, z2, x, y)
G <- Extreme_causality_graph(w, max_causal_lag = 2)  # Try it out also with lag = 1. You will see that the lagged edges disappear
```

One can then use, for example, the `igraph` package for plotting the estimated graph.

```r
# Visualizing the final graph estimate using the igraph package
library(igraph)

graph <- graph_from_edgelist(G$G)
V(graph)$name <- names(w)
plot.igraph(
  graph, 
  layout = layout_nicely(graph), 
  vertex.label = V(graph)$name,
  margin = c(0, 0, 0, 0),
  vertex.color = "white",
  vertex.label.color = "black",
  vertex.size = 30,
  edge.color = "black"
)
```


## References

Bodik, J. and Pasche, O. C. (2024). "Granger Causality in Extremes." *ArXiv Preprint* 2407.09632. [doi:10.48550/arXiv.2407.09632](https://doi.org/10.48550/arXiv.2407.09632).

