# Crypto-Stock application of Granger causality in extremes

source("./R/Main_functions.R")
# source("./R/old_Main_functions_Juraj.R")
source("./R/utils.R")
library(readxl)
library(igraph)


## Data import and preprocessing

dataset_path <- "./data/Crypto_Application/data_last_day.csv"
results_folder <- "./Results/Crypto_Application/"

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


# Just the simple plot
plot(data$x9, type = "l", col = 2, lty = 1, lwd = 2, xlab = "Time [minutes]", ylab = "Log-return values", main = "IOTA and Binance")
lines(data$x1, type = "l", lty = 1, lwd = 2)


## Graph estimation

# 1. Estimated graph using Algorithm 2 incorporating Algorithm 1
Estimated_graph_1 <- Extreme_causality_full_graph_estimate(w = data, lag_future = 1, both_tails = TRUE)

safe_save_rds(Estimated_graph_1, file = paste0(results_folder, "Estimated_graph_1.rds"))

# We use igraph package to visualize it:
dumping_factor <- 0 # increase if interested only in the strongest edges
G <- list(G = Estimated_graph_1$G, weights = Estimated_graph_1$weights)
G$G <- Estimated_graph_1$G[which(Estimated_graph_1$weights >= dumping_factor), ]
G$weights <- Estimated_graph_1$weights[which(Estimated_graph_1$weights >= dumping_factor)]
graph <- graph_from_edgelist(G$G)
E(graph)$edge.width <- 10 * (G$weights)
V(graph)$name <- c("Binance Coin", "Bitcoin", "BCH", "Cardano", "Dogecoin", "EOS.IO", "Ethereum", "Ethereum Classic", "IOTA", "Litecoin", "Maker", "Monero", "Stellar", "TRON")
plot.igraph(graph,
  layout = layout_with_fr(graph),
  vertex.label = V(graph)$name,
  edge.width = E(graph)$edge.width,
  vertex.color = "lightgrey",
  vertex.label.cex = 1,
  vertex.size = 22,
  main = "Applying Algorithm 1, lag = 1 min"
)



# 2. Estimated graph using Algorithm 2 incorporating the testing procedure, where an edge is present only if its corresponding p-value is below 0.05
# (Note that this takes a few hours to compute)
Estimated_graph_2 <- Extreme_causality_full_graph_estimate(w = data, lag_future = 1, both_tails = TRUE, p_value_based = TRUE, p_value_cutoff = 0.06)

safe_save_rds(Estimated_graph_2, file = paste0(results_folder, "Estimated_graph_2.rds"))

# We use igraph package to visualize it:
G <- list(G = Estimated_graph_2$G, weights = Estimated_graph_2$weights)
G$G <- Estimated_graph_2$G[which(Estimated_graph_2$weights >= dumping_factor), ]
G$weights <- Estimated_graph_2$weights[which(Estimated_graph_2$weights >= dumping_factor)]
graph <- graph_from_edgelist(G$G)
E(graph)$edge.width <- 10 * (G$weights)
V(graph)$name <- c(
  "Binance Coin", "Bitcoin", "BCH", "Cardano", "Dogecoin",
  "EOS.IO", "Ethereum", "Ethereum Classic", "IOTA", "Litecoin",
  "Maker", "Monero", "Stellar", "TRON"
)
plot.igraph(graph)


# 3. Finally, the same with lag = 30
Estimated_graph_3 <- Extreme_causality_full_graph_estimate(
  w = data, lag_future = 30, both_tails = TRUE,
  p_value_based = TRUE, p_value_cutoff = 0.05
)

safe_save_rds(Estimated_graph_3, file = paste0(results_folder, "Estimated_graph_3.rds"))

G <- list(G = Estimated_graph_3$G, weights = Estimated_graph_3$weights)
G$G <- Estimated_graph_3$G[which(Estimated_graph_3$weights >= dumping_factor), ]
G$weights <- Estimated_graph_3$weights[which(Estimated_graph_3$weights >= dumping_factor)]
graph <- graph_from_edgelist(G$G)
E(graph)$edge.width <- 10 * (G$weights)
V(graph)$name <- c(
  "Binance Coin", "Bitcoin", "BCH", "Cardano", "Dogecoin",
  "EOS.IO", "Ethereum", "Ethereum Classic", "IOTA", "Litecoin",
  "Maker", "Monero", "Stellar", "TRON"
)
plot.igraph(graph)



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
