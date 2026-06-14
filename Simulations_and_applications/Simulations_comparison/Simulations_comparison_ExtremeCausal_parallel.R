# Simulations for our Extreme Causal method, for comparison with other state-of-the-art methods on VAR and GARCH models.
# The competitor PCMCI methods are implemented in Python, using the Tigramite package; the code can be found in the corresponding python script.

# generate_data is the function that generated the data with sample size n, number of variables p (p=m from the manuscript) and structure='VAR' or 'GARCH' and heavy_tailed=TRUE or FALSE.
# For generating random graph we use the function erdos.renyi.game(p, 1/p, directed = TRUE,loops = FALSE) from igraph library
# We generate data with random graph + estimates the graph using our method + compute the distance between true graph and estimated graph + repeat 100 times and return the mean of the distances
# Results are saved for later processing.

source("./R/Main_functions.R")
source("./R/utils.R")
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
lag_future <- 1
lag_past <- 0


# Define the ranges for n and p
number_of_simulations <- 100 # 100
sequence_of_p <- c(2, 4, 7, 10, 20) # c(2, 4, 7, 10, 20)
sequence_of_n <- c(500, 5000) # c(500, 5000)

strategy <- "multisession" # "sequential", "multisession", "multicore"
n_workers <- 12 # future::availableCores() - 1

# Folder where results will be stored
results_folder <- "./Results/Simulations_comparison/tables/" # folder where results will be stored

check_directory(results_folder, recursive = TRUE, no_warning = TRUE)

set.seed(0)

start_time <- Sys.time()


generate_data <- function(n = 200, p = 4, structure = "VAR", heavy_tailed = TRUE) {
  generate_data_with_possible_NA <- function(n, p, structure, heavy_tailed) { # If the graph is too dense in high dimensions it can happen that the time series are not stationary but each X_{t+1}>X_t causing it to exponentially increase. How I dealt with it was that if this happens we just repeat while this is not the case. It is time-exhausting a bit but works

    x <- data.frame(matrix(rep(0, n * p), nrow = n))
    epsilon <- data.frame(matrix(rep(0, n * p), nrow = n))
    true_graph <- erdos.renyi.game(p, 1 / p, directed = TRUE, loops = FALSE)


    if (heavy_tailed) {
      if (structure == "VAR") {
        for (i in 1:p) {
          epsilon[, i] <- rpareto(n, 1, 1)
        }
      }
      if (structure == "GARCH") {
        for (i in 1:p) {
          epsilon[, i] <- rcauchy(n)
        }
      }
    } else {
      for (i in 1:p) {
        epsilon[, i] <- rnorm(n, 0, 1)
      }
    }


    if (structure == "VAR") {
      effect <- 0.3
      for (j in 3:n) {
        for (i in 1:p) {
          addition <- 0
          for (k in 1:p) {
            addition <- addition + (get.edge.ids(true_graph, c(k, i)) > 0) * effect * x[j - 1, k]
          }
          x[j, i] <- 0.3 * x[j - 1, i] + addition + epsilon[j, i]
        }
      }
    }


    if (structure == "GARCH") {
      effect <- 0.5
      if (heavy_tailed) {
        effect <- 05 
      }
      for (j in 3:n) {
        for (i in 1:p) {
          addition <- 0
          for (k in 1:p) {
            addition <- addition + (get.edge.ids(true_graph, c(k, i)) > 0) * effect * x[(j - sample(1:1, 1)), k]^2
          }
          x[j, i] <- ((0.1 + addition)^(0.5)) * epsilon[j, i]
        }
      }
    }

    return(list(data = x, true_graph = true_graph))
  }

  x <- generate_data_with_possible_NA(n = n, p = p, structure = structure, heavy_tailed = heavy_tailed)
  while (sum(is.na(x$data[n, ])) != 0) {
    x <- generate_data_with_possible_NA(n = n, p = p, structure = structure, heavy_tailed = heavy_tailed)
  }
  return(list(data = x$data, true_graph = x$true_graph))
}


distance_between_two_graphs <- function(graph1, graph2) {
  # graph1 <- graph_from_edgelist(graph1, directed = TRUE)
  # graph2 <- graph_from_edgelist(graph2, directed = TRUE)

  suppressWarnings({
    l <- length(E((graph1) %s% (graph2)))
  })

  return(length(E(graph1)) + length(E(graph2)) - 2 * l)
}


change_our_output_into_graph <- function(G, p, vtx_names = NULL) {
  if (is.list(G)) {
    G$G <- rbind(G$G, c(p, p))
    graph <- graph_from_edgelist(matrix(G$G, ncol = 2), directed = TRUE)
    # V(graph)$name <- names(x$data)
    if (!is.null(vtx_names)) {
      V(graph)$name <- vtx_names
    }
    graph <- graph - E(graph)[length(E(graph))]
    G$G <- G$G[-nrow(G$G), ]
  } else {
    graph <- make_empty_graph(p)
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

  x <- generate_data(n = n, p = p, structure = structure, heavy_tailed = heavy_tailed)
  
  G <- Extreme_causality_full_graph_estimate(x$data,
    lag_future = lag_future, lag_past = lag_past,
    nu_x = nu_x, q_y = q_y, q_z = q_z,
    both_tails = TRUE, p_value_based = FALSE
  )
  
  graph <- change_our_output_into_graph(G, p, vtx_names = names(x$data))
  one_result <- distance_between_two_graphs(graph, x$true_graph)

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


# for (n in sequence_of_n) {
#   final_average_result_for_n <- c()
#   for (p in sequence_of_p) {
#     result_for_p <- 0
#     for (one_simulation in seq(number_of_simulations)) {
#       x <- generate_data(n = n, p = p, structure = structure, heavy_tailed = heavy_tailed)
#       G <- Extreme_causality_full_graph_estimate(x$data)
#       graph <- change_our_output_into_graph(G, p, vtx_names = names(x$data))
#       one_result <- distance_between_two_graphs(graph, x$true_graph)

#       result_for_p <- c(result_for_p, one_result)
#     }
#     cat("p=", p, " and result =", mean(result_for_p), "\n")
#     final_average_result_for_n <- c(final_average_result_for_n, mean(result_for_p))
#   }
#   final_average_results[[as.character(n)]] <- final_average_result_for_n
# }


end_time <- Sys.time()
print(end_time - start_time)
# write execution time to a file
total_exec_time <- difftime(end_time, start_time, units = "hours")
write(paste("Total execution time:", as.numeric(total_exec_time), "hours"),
  file = paste0(results_folder, "execution_time_ExtremeCausal_", structure, "_", ht_str, ".txt")
)


# # Time measurement - how long does it take to estimate one graph?
# x <- generate_data(n = 500, p = 20, structure = structure, heavy_tailed = heavy_tailed)
# start_time <- Sys.time()

# G <- Extreme_causality_full_graph_estimate(x$data)

# end_time <- Sys.time()
# print(end_time - start_time)


# # Wait for 2 mins and then shutdown the machine
# Sys.sleep(120)
# system('shutdown -t 30 -s')# For Windows
# # system('sudo poweroff')# For Linux (with non-root sudo user)
