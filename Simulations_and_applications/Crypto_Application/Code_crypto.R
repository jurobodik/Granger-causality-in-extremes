# Crypto-Stock application of Granger causality in extremes

source("./R/Main_functions.R")
# source("./R/old_Main_functions_Juraj.R")
source("./R/utils.R")
library(EnvStats)
library(readxl)
library(knitr)
library(foreach)
library(doFuture)
library(igraph)


do_plots <- TRUE
force_recompute <- FALSE
dataset_choice <- "last_10days" # "last_day" or "last_10days"

use_parallel_doFuture <- TRUE
strategy <- "multisession" # "sequential", "multisession", "multicore"
n_workers <- 12 # future::availableCores() - 1

## Data import and preprocessing
results_folder <- "./Results/Crypto_Application/"
# save_folder <- paste0(results_folder, "graph_saves/")
save_folder <- paste0(results_folder, "graph_saves_parallel/")
plot_folder <- paste0(results_folder, "figures/")

check_directory(save_folder, recursive = TRUE, no_warning = TRUE)
check_directory(plot_folder, recursive = TRUE, no_warning = TRUE)

if(dataset_choice == "last_day") {
  dataset_path <- "./data/Crypto_Application/data_last_day.csv"
} else if(dataset_choice == "last_10days") {
  dataset_path <- "./data/Crypto_Application/data_last_10days_nona.csv"
} else {
  stop("Invalid dataset choice. Please choose either 'last_day' or 'last_10days'.")
}

set.seed(0)

data_origin <- read.csv(file = dataset_path)
head(data_origin)

log_returns <- function(prices) {
  return(diff(log(prices)))
}

data <- c()
for (i in 1:14) {
  # only works if the data is ordered by timestamp and there are no missing values for part of the assets.
  x <- data_origin[data_origin$Asset_ID == i, 4]
  x <- log_returns(x)
  data <- cbind(data, x)
}
data <- as.data.frame(data)
names(data) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14")
head(data)
data <- abs(data) # we take absolute values since we both tails are of interest

asset_names <- c(
  "Binance Coin", "Bitcoin", "BCH", "Cardano", "Dogecoin",
  "EOS.IO", "Ethereum", "Ethereum Classic", "IOTA", "Litecoin",
  "Maker", "Monero", "Stellar", "TRON"
)

# Just the simple plot
if(do_plots){
  plot(data$x9, type = "l", col = 2, lty = 1, lwd = 2, xlab = "Time [minutes]", ylab = "Log-return values", main = "IOTA and Binance")
  lines(data$x1, type = "l", lty = 1, lwd = 2)
}


## Graph estimation

# output graph postprocessing function
postprocess_estimated_graph <- function(Estimated_graph, dumping_factor = 0, V_names = asset_names) {
  G <- list(G = Estimated_graph$G, weights = Estimated_graph$weights)
  G$G <- Estimated_graph$G[which(Estimated_graph$weights >= dumping_factor), ]
  G$weights <- Estimated_graph$weights[which(Estimated_graph$weights >= dumping_factor)]
  graph <- igraph::graph_from_edgelist(G$G)
  igraph::E(graph)$edge.width <- 10 * (G$weights)
  igraph::V(graph)$name <- V_names
  return(graph)
}

plot_estimated_graph <- function(graph, graph_layout=NULL, ..., save_path = NULL) {
  edge_widths <- igraph::E(graph)$edge.width
  edge_widths <- edge_widths/max(edge_widths) * 3 + 0.5
  
  if(is.null(graph_layout)){
    # graph_layout <- igraph::layout_with_fr(graph)
    # graph_layout <- igraph::layout_with_fr(graph, weight = edge_widths)
    # graph_layout <- igraph::layout_with_fr(graph, coords = igraph::layout_as_tree(graph),
    #                                        niter=1000, grid="grid")
    # graph_layout <- igraph::layout_as_tree(graph)
    # graph_layout <- igraph::layout_nicely(graph)
    # graph_layout <- igraph::layout_with_sugiyama(graph, hgap=5, vgap=5, maxiter=1000)
    graph_layout <- igraph::layout_with_graphopt(graph, start = igraph::layout_as_tree(graph),
                                                 charge = 0.5, mass = 30, niter = 1000)
    # graph_layout <- igraph::layout_with_kk(graph)
    # graph_layout <- igraph::layout_on_grid(graph, width=5)
  }
  if(!is.null(save_path)) {
    pdf(file = save_path) #, width = 8, height = 8)
  }
  igraph::plot.igraph(graph,
    layout = graph_layout, 
    margin = c(0, 0, 0, 0),
    asp = 0.7,
    # vertex.label = igraph::V(graph)$name,
    vertex.label = gsub(' ', '\n ', igraph::V(graph)$name),
    # vertex.color = "lightgrey",
    vertex.color = "white",
    vertex.label.color = "black",
    vertex.label.cex = 1,
    # vertex.size = 22,
    vertex.size = 30,
    edge.width = edge_widths,
    edge.color = "black",
    # edge.arrow.mode = 3,
    # edge.curved = 0.2,
    edge.curved = 0,
    edge.arrow.size = 1.2,
    edge.arrow.width = 0.8,
    # main = title,
    ...
  )
  if(!is.null(save_path)){
    dev.off()
    knitr::plot_crop(save_path)
  } 
}

tkplot_estimated_graph <- function(graph, graph_layout=NULL, ...) {
  edge_widths <- igraph::E(graph)$edge.width
  edge_widths <- edge_widths/max(edge_widths) * 3 + 0.5
  
  if(is.null(graph_layout)){
    graph_layout <- igraph::layout_with_graphopt(graph, start = igraph::layout_as_tree(graph),
                                                 charge = 0.5, mass = 30, niter = 1000)
  }
  
  igraph::tkplot(graph,
    layout = graph_layout, 
    margin = c(0, 0, 0, 0),
    asp = 0.7,
    # vertex.label = igraph::V(graph)$name,
    vertex.label = gsub(' ', '\n ', igraph::V(graph)$name),
    # vertex.color = "lightgrey",
    vertex.color = "white",
    vertex.label.color = "black",
    vertex.label.cex = 1,
    # vertex.size = 22,
    vertex.size = 30,
    edge.width = edge_widths,
    edge.color = "black",
    # edge.arrow.mode = 3,
    # edge.curved = 0.2,
    edge.curved = 0,
    edge.arrow.size = 1.2,
    edge.arrow.width = 0.8,
    # main = title,
    ...
  )
}




# 1. Estimated graph using Algorithm 2 incorporating Algorithm 1
checkpoint_path_1 <- paste0(save_folder, "Estimated_graph_", dataset_choice, "_A2_l1_A1.rds")
if(force_recompute || !file.exists(checkpoint_path_1)){
  if(use_parallel_doFuture){
    Estimated_graph_1 <- Extreme_causality_graph_estimate_parallel(
      w = data, lag_future = 1, both_tails = TRUE,
      strategy = strategy, n_workers = n_workers
    )
  } else {
    Estimated_graph_1 <- Extreme_causality_full_graph_estimate(
      w = data, lag_future = 1, both_tails = TRUE
    )
  }
  safe_save_rds(Estimated_graph_1, file = checkpoint_path_1)
} else {
  Estimated_graph_1 <- readRDS(checkpoint_path_1)
}



# 2. Estimated graph using Algorithm 2 incorporating the testing procedure, where an edge is present only if its corresponding p-value is below 0.05
# (Note that this takes a few hours to compute)
checkpoint_path_2 <- paste0(save_folder, "Estimated_graph_", dataset_choice, "_A2_l1_p05.rds")
if(force_recompute || !file.exists(checkpoint_path_2)){
  if(use_parallel_doFuture){
    Estimated_graph_2 <- Extreme_causality_graph_estimate_parallel(
      w = data, lag_future = 1, both_tails = TRUE, 
      p_value_based = TRUE, p_value_cutoff = 0.05,
      # p_value_based = TRUE, p_value_cutoff = 0.06,
      strategy = strategy, n_workers = n_workers
    )
  } else {
    Estimated_graph_2 <- Extreme_causality_full_graph_estimate(
      w = data, lag_future = 1, both_tails = TRUE, 
      p_value_based = TRUE, p_value_cutoff = 0.05
      # p_value_based = TRUE, p_value_cutoff = 0.06
    )
  }
  safe_save_rds(Estimated_graph_2, file = checkpoint_path_2)
} else {
  Estimated_graph_2 <- readRDS(checkpoint_path_2)
}



# 3. The same (Algorithm 2 with test procedure), with lag = 30
checkpoint_path_3 <- paste0(save_folder, "Estimated_graph_", dataset_choice, "_A2_l30_p05.rds")
if(force_recompute || !file.exists(checkpoint_path_3)){
  if(use_parallel_doFuture){
    Estimated_graph_3 <- Extreme_causality_graph_estimate_parallel(
      w = data, lag_future = 30, both_tails = TRUE,
      p_value_based = TRUE, p_value_cutoff = 0.05,
      strategy = strategy, n_workers = n_workers
    )
  } else {
    Estimated_graph_3 <- Extreme_causality_full_graph_estimate(
      w = data, lag_future = 30, both_tails = TRUE,
      p_value_based = TRUE, p_value_cutoff = 0.05
    )
  }
  safe_save_rds(Estimated_graph_3, file = checkpoint_path_3)
} else {
  Estimated_graph_3 <- readRDS(checkpoint_path_3)
}



# 4. Finally, using Algorithm 2 incorporating Algorithm 1 with lag 30
set.seed(123) # It was run separately
checkpoint_path_4 <- paste0(save_folder, "Estimated_graph_", dataset_choice, "_A2_l30_A1.rds")
if(force_recompute || !file.exists(checkpoint_path_4)){
  if(use_parallel_doFuture){
    Estimated_graph_4 <- Extreme_causality_graph_estimate_parallel(
      w = data, lag_future = 30, both_tails = TRUE,
      strategy = strategy, n_workers = n_workers
    )
  } else {
    Estimated_graph_4 <- Extreme_causality_full_graph_estimate(
      w = data, lag_future = 30, both_tails = TRUE
    )
  }
  safe_save_rds(Estimated_graph_4, file = checkpoint_path_4)
} else {
  Estimated_graph_4 <- readRDS(checkpoint_path_4)
}



if(do_plots){
  # We use igraph package to visualize it:
  dumping_factor <- 0 # increase if interested only in the strongest edges
  graph_1 <- postprocess_estimated_graph(Estimated_graph_1, dumping_factor = dumping_factor, V_names = asset_names)
  graph_2 <- postprocess_estimated_graph(Estimated_graph_2, dumping_factor = dumping_factor, V_names = asset_names)
  graph_3 <- postprocess_estimated_graph(Estimated_graph_3, dumping_factor = dumping_factor, V_names = asset_names)
  graph_4 <- postprocess_estimated_graph(Estimated_graph_4, dumping_factor = dumping_factor, V_names = asset_names)
  
  # graphs_default_layout <- igraph::layout_with_graphopt(graph_2, start = igraph::layout_as_tree(graph_2),
  #                                                       charge = 0.5, mass = 30, niter = 1000)
  # graphs_default_layout <- graphs_default_layout[c(4,7,3,1,5,6,2,8,14,10,11,12,13,9),]
  # graphs_default_layout <- rbind(c(-248.939187,  814.0000),  # "Binance Coin"
  #                                c(  12.964700,  293.0000),  # "Bitcoin"
  #                                c(-278.574257, -223.7301),  # "BCH"
  #                                c(-732.500000,  673.0000),  # "Cardano"
  #                                c(-481.500000,  227.2756),  # "Dogecoin"
  #                                c( 390.000000, -945.5917),  # "EOS.IO"
  #                                c(  -9.695144, -658.0000),  # "Ethereum"
  #                                c(-771.500000, -181.2956),  # "Ethereum Classic"
  #                                c( 743.500000,  700.7325),  # "IOTA"
  #                                c( 250.295514,  734.0000),  # "Litecoin"
  #                                c( 711.500000, -234.8546),  # "Maker"
  #                                c( 524.532835,  249.0118),  # "Monero"
  #                                c(-554.930701, -640.0000),  # "Stellar"
  #                                c( 226.364516, -157.0000))  # "TRON"
  # graphs_default_layout <- rbind(c(-248.93,  814.00),  # "Binance Coin"
  #                                c(  12.96,  293.00),  # "Bitcoin"
  #                                c(-335.57, -233.73),  # "BCH"
  #                                c(-732.50,  673.00),  # "Cardano"
  #                                c(-481.50,  227.27),  # "Dogecoin"
  #                                c( 420.00, -650.00),  # "EOS.IO"
  #                                c(  -9.69, -658.00),  # "Ethereum"
  #                                c(-771.50, -181.29),  # "Ethereum Classic"
  #                                c( 743.50,  700.73),  # "IOTA"
  #                                c( 250.29,  734.00),  # "Litecoin"
  #                                c( 701.50, -234.85),  # "Maker"
  #                                c( 554.53,  249.01),  # "Monero"
  #                                c(-490.00, -640.00),  # "Stellar"
  #                                c( 290.00, -140.00))  # "TRON"
  
  graphs_default_layout <- rbind(c(  12.69, -598.00),  # "Binance Coin"
                                 c(-248.93,  814.00),  # "Bitcoin"
                                 c(-251.93, -220.73),  # "BCH"
                                 c( 420.00, -650.00),  # "Cardano"
                                 c(-481.50,  227.27),  # "Dogecoin"
                                 c(-732.50,  673.00),  # "EOS.IO"
                                 c(  12.96,  293.00),  # "Ethereum"
                                 c(-771.50, -210.29),  # "Ethereum Classic"
                                 c( 743.50,  700.73),  # "IOTA"
                                 c( 250.29,  734.00),  # "Litecoin"
                                 c( 290.00, -140.00),  # "Maker"
                                 c( 674.53,  249.01),  # "Monero"
                                 c(-490.00, -640.00),  # "Stellar"
                                 c( 701.50, -234.85))  # "TRON"
  
  
  plot_estimated_graph(graph_1, graph_layout = graphs_default_layout, main = "Using Algorithm 1 and 1 min lag")
  plot_estimated_graph(graph_1, graph_layout = graphs_default_layout, # main = "Using Algorithm 1 and 1 min lag", 
                       save_path = paste0(plot_folder, "Estimated_graph_", dataset_choice, "_A2_l1_A1.pdf"))
  # tkplot_estimated_graph(graph_1, graph_layout = graphs_default_layout)
  
  plot_estimated_graph(graph_2, graph_layout = graphs_default_layout, main = "Using p-value based edge selection and 1 min lag")
  plot_estimated_graph(graph_2, graph_layout = graphs_default_layout, # main = "Using p-value based edge selection and 1 min lag", 
                       save_path = paste0(plot_folder, "Estimated_graph_", dataset_choice, "_A2_l1_p05.pdf"))
  # tkplot_estimated_graph(graph_2, graph_layout = graphs_default_layout)
  
  plot_estimated_graph(graph_3, graph_layout = graphs_default_layout, main = "Using p-value based edge selection and 30 min lag")
  plot_estimated_graph(graph_3, graph_layout = graphs_default_layout, # main = "Using p-value based edge selection and 30 min lag", 
                       save_path = paste0(plot_folder, "Estimated_graph_", dataset_choice, "_A2_l30_p05.pdf"))
  # tkplot_estimated_graph(graph_3, graph_layout = graphs_default_layout)
  
  plot_estimated_graph(graph_4, graph_layout = graphs_default_layout, main = "Using Algorithm 2 and 30 min lag")
  plot_estimated_graph(graph_4, graph_layout = graphs_default_layout, # main = "Using Algorithm 2 and 30 min lag", 
                       save_path = paste0(plot_folder, "Estimated_graph_", dataset_choice, "_A2_l30_A1.pdf"))
  # tkplot_estimated_graph(graph_4, graph_layout = graphs_default_layout)
}



# Fun fact: We asked ChatGPT on May 2024 to draw causal graph between the variables. This is the result:

# Bitcoin --> Ethereum --> ERC-20 Tokens (e.g., Maker, IOTA)
#|           |
#             +--> Ethereum Classic
#|
#  +--> Litecoin --> Dogecoin
#|
#  +--> Monero
#|
#  +--> Binance Coin
#|
#  +--> Cardano --> Stellar
#|
#  +--> BCH
#|
#  +--> EOS.IO
