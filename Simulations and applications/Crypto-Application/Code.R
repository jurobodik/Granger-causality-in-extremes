#Crypto-Stock application of Granger causality in extremes 

#main function needed to be uploaded for this application is 'Extreme_causality_full_graph_estimate' - this function can be found in the main file. 

library("readxl")
library(igraph)


data_origin = read.csv(file = 'data_last_day.csv')
head(data_origin)

log_returns <- function(prices) {return(diff(log(prices)))}

data = c()
for (i in 1:14) {
  x = data_origin[data_origin$Asset_ID==i,4]
  x = log_returns(x)
  data = cbind(data,x ) 
}
data = as.data.frame(data)
names(data) = c('x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8', 'x9', 'x10', 'x11', 'x12', 'x13', 'x14')
head(data)


#Just the simple plot 
plot(data$x9, type = 'l', col=2, lty = 1, lwd = 2, xlab = 'Time [minutes]', ylab = 'Log-return values', main = 'IOTA and Binance')
lines(data$x1, type = 'l', lty = 1, lwd = 2)



lag_future=1; 
data = abs(data) 

Estimated_graph=Extreme_causality_full_graph_estimate(w=data, lag_future = lag_future, both_tails = TRUE)



#We use igraph package to visualize it:
dumping_factor = 0; #First figure
dumping_factor = 0.2;  #Second figure

G=list(G= Estimated_graph$G, weights = Estimated_graph$weights)
G$G <- Estimated_graph$G[which(Estimated_graph$weights >= dumping_factor), ]; G$weights <- Estimated_graph$weights[which(Estimated_graph$weights >= dumping_factor)]
graph <- graph_from_edgelist(G$G)
E(graph)$edge.width <- 10*(G$weights)
V(graph)$name =  c("Binance Coin", "Bitcoin", "BCH", "Cardano", "Dogecoin", "EOS.IO", "Ethereum", "Ethereum Classic", "IOTA", "Litecoin", "Maker", "Monero","Stellar", "TRON" )


plot.igraph(graph, 
            layout = layout_with_fr(graph), 
            vertex.label = V(graph)$name, 
            edge.width = E(graph)$edge.width, 
            vertex.color = 'lightgrey', 
            vertex.label.cex = 1, 
            vertex.size =22, 
            main = 'Applying Algorithm 1, lag = 1 min')




#If you want to obtain the estimates based on p-values, you can just try copmuting it for some specific edges that you are interested in
#In general, the size of Gamma coefficients correspond to the p-values, so it is not necessary to compute all of them. 
#For example, is there a line from x1->x2?
x = data$x1; y = data$x2; z = data.frame(data$x3, data$x4, data$x5, data$x6, data$x7, data$x8, data$x9, data$x10, data$x11, data$x12, data$x13, data$x14)
Extreme_causality_test(x = x, y = y, z = z, lag_future = 30,  p_value_computation = TRUE)





#We do the same with lag = 30
lag_future=30
Estimated_graph=Extreme_causality_full_graph_estimate(w=data, lag_future = lag_future)



#We use igraph package to visualize it:
dumping_factor = 0.3
G=list(G= Estimated_graph$G, weights = Estimated_graph$weights)
G$G <- Estimated_graph$G[which(Estimated_graph$weights >= dumping_factor), ]; G$weights <- Estimated_graph$weights[which(Estimated_graph$weights >= dumping_factor)]
graph <- graph_from_edgelist(G$G)
E(graph)$edge.width <- 10*(G$weights)
V(graph)$name =  c("Binance Coin", "Bitcoin", "BCH", "Cardano", "Dogecoin", "EOS.IO", "Ethereum", "Ethereum Classic", "IOTA", "Litecoin", "Maker", "Monero","Stellar", "TRON" )


plot.igraph(graph, 
            layout = layout_with_fr(graph), 
            vertex.label = V(graph)$name, 
            edge.width = E(graph)$edge.width, 
            vertex.color = 'lightgrey', 
            vertex.label.cex = 1, 
            vertex.size =22, 
            main = 'Applying Algorithm 1, lag = 1 min')





#Fun fact: We asked ChatGPT on May 2024 to draw causal graph between the variables. This is the result: 

#Bitcoin --> Ethereum --> ERC-20 Tokens (e.g., Maker, IOTA, TRON, etc.)
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


