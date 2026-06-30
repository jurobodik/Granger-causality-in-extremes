# Simulations for our Extreme Causal method, for comparison with other state-of-the-art methods on VAR and GARCH models.
# The competitor PCMCI methods are implemented in Python, using the Tigramite package; the code can be found in the corresponding python script.

# generate_data is the function that generated the data with sample size n, number of variables p (p=m from the manuscript) and structure='VAR' or 'GARCH' and heavy_tailed=TRUE or FALSE.
# For generating random graph we use the function erdos.renyi.game(p, 1/p, directed = TRUE,loops = FALSE) from igraph library
# We generate data with random graph + estimates the graph using our method + compute the distance between true graph and estimated graph + repeat 100 times and return the mean of the distances
# Results are saved for later processing.

# (Install and) load v0.1.0 of the ExtremeGranger R package
# devtools::install_github("opasche/ExtremeGranger@v0.1.0")
library(ExtremeGranger)
# source("./R/Main_functions_localcopy.R")

source("./R/utils.R")

# Other packages used in the analysis
library(EnvStats)
library(ExtremeRisks)
library(foreach)
library(doFuture)
library(igraph)

## ==== PARAMETERS ====

# Scenario parameters
structure <- "VAR" # "VAR" or "GARCH"
heavy_tailed <- TRUE # TRUE or FALSE

nu_x <- 0.3
q_y <- 0.2
q_z <- 0.1
max_causal_lag <- 1
max_confounder_lag <- 0


# Define the ranges for n and p
number_of_simulations <- 100 # 100
sequence_of_p <- c(2, 4, 7, 10, 20) # c(2, 4, 7, 10, 20)
sequence_of_n <- c(500, 5000) # c(500, 5000)

strategy <- "multisession" # "sequential", "multisession", "multicore"
n_workers <- 6 # future::availableCores() - 1

# Folder where results will be stored
results_folder <- "./Results/Simulations_comparison/tables/" # folder where results will be stored

check_directory(results_folder, recursive = TRUE, no_warning = TRUE)

set.seed(0)

start_time <- Sys.time()


generate_data <- function(n = 200, p = 4, structure = "VAR", heavy_tailed = TRUE) {
  generate_data_internal <- function(n, p, structure, heavy_tailed) {

    true_graph <- erdos.renyi.game(p, 1 / p, directed = TRUE, loops = FALSE)
    X <- matrix(0, nrow = n, ncol = p)

    if (structure == "VAR") {
      ar_coef <- 0.1
      cross_coef <- 0.5
      noise_sd <- 1
      normalize_by_indegree <- FALSE
      
      if (heavy_tailed) {
        epsilon <- matrix(rpareto(n * p, 1, 1), nrow = n, ncol = p)
      } else {
        epsilon <- matrix(rnorm(n * p, sd = noise_sd), nrow = n, ncol = p)
      }
      
      B <- matrix(0, p, p) # B[i,k] = effect of source k on target i
      el <- as_edgelist(true_graph) # rows: (source, target)
      if (nrow(el) > 0) {
        for (e in seq_len(nrow(el))) {
          B[el[e, 2], el[e, 1]] <- cross_coef
        }
      }
      if (normalize_by_indegree) {
        indeg <- rowSums(B != 0)
        indeg[indeg == 0] <- 1
        B <- B / indeg
      }
      diag(B) <- ar_coef
      
      for (ti in 2:n) {
        X[ti, ] <- as.vector(B %*% X[ti - 1, ]) + epsilon[ti, ]
      }
      X <- data.frame(X)
    }


    else if (structure == "GARCH") {
      if (heavy_tailed) {
        epsilon <- matrix(rcauchy(n * p), nrow = n, ncol = p)
      } else {
        epsilon <- matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
      }
      
      effect <- 0.5
      
      # Create adjacency matrix from the graph: A[i,k] = 1 if edge from k to i exists
      A <- matrix(0, p, p)
      el <- as_edgelist(true_graph)
      if (nrow(el) > 0) {
        for (e in seq_len(nrow(el))) {
          # A[el[e, 2], el[e, 1]] <- effect
          A[el[e, 2], el[e, 1]] <- 1
        }
      }
      
      # Vectorized loop for GARCH dynamics
      for (ti in 2:n) {
        # X_prev_squared <- X[ti - 1, ]^2
        # addition <- as.vector(A %*% (X[ti - 1, ]^2))
        addition <- as.vector(A %*% (X[ti - 1, ]^2)) * effect
        X[ti, ] <- sqrt(0.1 + addition) * epsilon[ti, ]
      }
      
      X <- data.frame(X)
    }
    
    else {
      stop("Invalid structure type. Must be 'VAR' or 'GARCH'.")
    }
    
    # if (all(is.finite(X)) && max(abs(X[n, ])) < max_abs) {
    #   colnames(X) <- paste0("X", seq_len(p))
    #   return(list(data = as.data.frame(X), true_graph = true_graph))
    # }

    return(list(data = X, true_graph = true_graph))
  }

  X <- generate_data_internal(n = n, p = p, structure = structure, heavy_tailed = heavy_tailed)
  trynm <- 0
  # while (!(all(is.finite(X$data)) && max(abs(X$data[n, ])) < 1e9)) {
  while (sum(is.na(X$data[n, ])) != 0) {
    X <- generate_data_internal(n = n, p = p, structure = structure, heavy_tailed = heavy_tailed)
    # DEBUG CHECK: to be sure it is numerically stationary.
    trynm <- trynm + 1
    if (trynm > 100) {
      cat("Non-stationary DGP (n =", n, ", p =", p, ") -- Failed to generate valid data after 100 attempts.\n")
      # warning(paste0("Non-stationary DGP (n =", n, ", p =", p, ") -- Failed to generate valid data after 100 attempts.\n"))
      stop(paste0("Non-stationary DGP (n =", n, ", p =", p, ") -- Failed to generate valid data after 100 attempts.\n"))
    }
  }
  return(list(data = X$data, true_graph = X$true_graph))
}


distance_between_two_graphs <- function(graph1, graph2) {
  # graph1 <- graph_from_edgelist(graph1, directed = TRUE)
  # graph2 <- graph_from_edgelist(graph2, directed = TRUE)
  suppressWarnings({
    common <- length(E((graph1) %s% (graph2)))
  })
  return(length(E(graph1)) + length(E(graph2)) - 2 * common)
}



change_our_output_into_graph <- function(G, p, vtx_names = NULL) {
  if (is.list(G) && !is.null(G$G) && length(G$G) > 0) {
    G$G <- rbind(G$G, c(p, p))
    graph <- graph_from_edgelist(matrix(G$G, ncol = 2), directed = TRUE)
    # V(graph)$name <- names(X$data)
    if (!is.null(vtx_names)) {
      V(graph)$name <- vtx_names
    }
    graph <- graph - E(graph)[length(E(graph))]
    G$G <- G$G[-nrow(G$G), ]
    # graph <- make_empty_graph(n = p, directed = TRUE)
    # graph <- add_edges(graph, as.vector(t(matrix(G$G, ncol = 2))))
  } else {
    graph <- make_empty_graph(n = p, directed = TRUE)
  }
  return(graph)
}


if (heavy_tailed) {
  ht_str <- "heavy"
} else {
  ht_str <- "light"
}


cat("Starting simulations for structure =", structure, " and heavy_tailed =", heavy_tailed, "\n")

# distance_results <- data.frame(
#   n = integer(),
#   p = integer(),
#   replication = integer(),
#   distance = numeric()
# )

index_pairs <- expand.grid("sim" = seq(number_of_simulations), "p" = sequence_of_p, "n" = sequence_of_n)[, c("n", "p", "sim")]

`%fun%` <- set_doFuture_strategy(strategy, n_workers)


distance_results <- foreach::foreach(i = seq(nrow(index_pairs)), .combine = rbind, .options.future = list(seed = TRUE)) %fun% {
  n <- index_pairs$n[i]
  p <- index_pairs$p[i]
  one_simulation <- index_pairs$sim[i]

  # cat("-- Starting simulation", one_simulation, "for n =", n, "and p =", p, "(at", difftime(Sys.time(), start_time, units = "mins"), "minutes).\n")

  X <- generate_data(n = n, p = p, structure = structure, heavy_tailed = heavy_tailed)
  
  G <- Extreme_causality_graph(X$data,
    max_causal_lag = max_causal_lag, max_confounder_lag = max_confounder_lag,
    nu_x = nu_x, q_y = q_y, q_z = q_z,
    both_tails = TRUE, p_value_based = FALSE
  )
  
  graph <- change_our_output_into_graph(G, p, vtx_names = names(X$data))
  one_result <- distance_between_two_graphs(graph, X$true_graph)

  data.frame(
    n = n,
    p = p,
    replication = one_simulation,
    distance = one_result
  )
}


# add columns for method, structure and heavy_tailed to distance_results
# then, reorder columns in distance_results
distance_results$method <- "ExtremeCausal"
distance_results$structure <- structure
distance_results$heavy_tailed <- heavy_tailed
distance_results <- distance_results[, c("method", "structure", "heavy_tailed", "n", "p", "replication", "distance")]

# Save results on disk
write.csv(distance_results, file = paste0(results_folder, "results_ExtremeCausal_", structure, "_", ht_str, ".csv"), row.names = FALSE)
distance_results

end_doFuture_strategy()



end_time <- Sys.time()
print(end_time - start_time)
# write execution time to a file
total_exec_time <- difftime(end_time, start_time, units = "hours")
write(paste("Total execution time:", as.numeric(total_exec_time), "hours"),
  file = paste0(results_folder, "execution_time_ExtremeCausal_", structure, "_", ht_str, ".txt")
)


# # Time measurement - how long does it take to estimate one graph?
# X <- generate_data(n = 500, p = 20, structure = structure, heavy_tailed = heavy_tailed)
# start_time <- Sys.time()

# G <- Extreme_causality_graph(X$data)

# end_time <- Sys.time()
# print(end_time - start_time)


