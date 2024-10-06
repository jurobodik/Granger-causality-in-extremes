#Simulations from section 6 where we compared our method with other state-of-the-art methods on VAR and GARCH models
#Python code can be found in the other file that deals with other methods from the Tigramite package

#generate_data is the function that generated the data with sample size n, number of variables p (p=m from the manuscript) and structure='VAR' or 'GARCH' and heavy_tailed=TRUE or FALSE. 
#For generating random graph we use the function erdos.renyi.game(p, 1/p, directed = TRUE,loops = FALSE) from igraph library
#Lines 100-150 are the main part. We generate data with random graph + estimates the graph using our method + compute the distance between true graph and estimated graph + repeat 100 times and return the mean of the distances
#We maually rewrote the values for each combination of 'structure' and 'heavy_tailed' into the table (easier since we dealt with python + R results)
#The graph with final results can be found here at lines 150+


library(EnvStats)
library(ExtremeRisks)
library(igraph)


generate_data = function( n=200, p=4, structure = 'VAR', heavy_tailed =TRUE){
  
  generate_data_with_possible_NA = function(n, p, structure, heavy_tailed){  #If the graph is too dense in high dimensions it can happen that the time series are not stationary but each X_{t+1}>X_t causing it to exponentially increase. How I dealt with it was that if this happens we just repeat while this is not the case. It is time-exhausting a bit but works
    
    x = data.frame(matrix(rep(0, n * p), nrow = n))
    epsilon = data.frame(matrix(rep(0, n * p), nrow = n))
    true_graph <- erdos.renyi.game(p, 1/p, directed = TRUE,loops = FALSE)
    
    
    
    if(heavy_tailed ==FALSE){
      for (i in 1:p) {epsilon[,i] = rnorm(n, 0,1)  }}
    if(heavy_tailed ==TRUE){
      if(structure == 'VAR') { for (i in 1:p) {epsilon[,i] = rpareto(n, 1,1)  }}
      if(structure == 'GARCH') { for (i in 1:p) {epsilon[,i] = rcauchy(n)  }}
    }
    
    
    
    
    if(structure == 'VAR'){
      effect = 0.3
      for (j in 3:n) { 
        for (i in 1:p) {
          addition = 0; for (k in 1:p) { addition = addition + (get.edge.ids(true_graph, c(k, i))>0)*effect*x[j-1,k] }
          x[j,i]  = 0.3*x[j-1,i] + addition + epsilon[j,i]
        }}
    }
    
    
    if(structure == 'GARCH'){
      effect = 0.5; if(heavy_tailed ==TRUE){effect = 05;}
      for (j in 3:n) { 
        for (i in 1:p) {
          addition = 0; for (k in 1:p) { addition = addition + (get.edge.ids(true_graph, c(k, i))>0)*effect*x[(j-sample(1:1, 1)),k]^2 }
          x[j,i]  =(    (0.1 + addition )^(0.5)     )*epsilon[j,i]
        }}
    }
    
    return( list(data = x, true_graph = true_graph)   ) 
  }
  
  x = generate_data_with_possible_NA(n=n, p=p, structure=structure, heavy_tailed=heavy_tailed)
  while (sum(is.na(x$data[n,]))!=0 ) {x = generate_data_with_possible_NA(n=n, p=p, structure=structure, heavy_tailed=heavy_tailed)}
  return( list(data = x$data, true_graph = x$true_graph)   ) 
}



distance_between_two_graphs = function(graph1, graph2){
  # graph1 <- graph_from_edgelist(graph1, directed = TRUE)
  #  graph2 <- graph_from_edgelist(graph2, directed = TRUE) 
  
  suppressWarnings({ l= length(E(   (graph1) %s% (graph2)  ))   }) 
  
  return(     length(E(graph1)) + length(E(graph2)) - 2*l      )
}



change_our_output_into_graph = function(G, p){
  if(is.list(G)){ G$G = rbind(G$G, c(p,p)); graph <- graph_from_edgelist(matrix(G$G, ncol = 2),  directed = TRUE);  V(graph)$name<- names(x$data); 
  graph = graph - E(graph)[length(E(graph))]; G$G = G$G[-nrow(G$G),];
  }else{ graph=   make_empty_graph(p)}
  return(graph)
}


















#Select structure and heavy_tailed and run the following code to obtain the results in one setting
structure = 'VAR'
heavy_tailed='TRUE'
number_of_simulations = 100



sequence_of_p = c(2, 4, 7, 10, 20)
sequence_of_n = c(500, 5000)
final_result =  list()

for (n in sequence_of_n) {
  
  final_result_for_n = c()
  for (p in sequence_of_p) {
    
    result_for_p = 0
    for (one_simulation in 1:number_of_simulations) {
      
      x=generate_data(n=n, p=p, structure = structure, heavy_tailed=heavy_tailed)
      G=Extreme_causality_full_graph_estimate(x$data)
      graph = change_our_output_into_graph(G, p)
      one_result =   distance_between_two_graphs(graph, x$true_graph)
      
      result_for_p = c(result_for_p, one_result)
      
    }  
    cat('p=', p,' and result =', mean(result_for_p), '\n')
    final_result_for_n = c(final_result_for_n, mean(result_for_p))
  }
  final_result[[as.character(n)]] = final_result_for_n
}

final_result
















#Here are the resulting numbers that you should get by running the previous codes 
#We divide everything by p(p-1) because there are p*(p-1) arrows to estimate
standardisation = sequence_of_p*(sequence_of_p-1) 


# Define the data vectors
data_vectors <- list(
  # OUR Method (VAR heavy), n=500 and n=5000 respectively
  c(0.0198, 0.1485, 0.732, 1.5346, 6.5754257) / standardisation,
  c(0.000, 0.0693, 0.1287, 0.4554, 1.485148) / standardisation,
  # OUR Method (GARCH heavy), n=500 and n=5000 respectively
  c(0.079, 0.316, 1.336, 2.940, 13.019) / standardisation,
  c(0.009, 0.069, 0.376, 1.128, 5.546) / standardisation,
  # OUR Method (VAR non-heavy), n=500 and n=5000 respectively
  c(0.69306, 1.93069, 4.5247, 7.5643, 23.831) / standardisation,
  c(0.38, 1.31, 2.95, 4.32, 9.81) / standardisation,
  # OUR Method (GARCH non-heavy), n=500 and n=5000 respectively
  c(0.603, 1.831, 4.287, 7.079, 20.336) / standardisation,
  c(0.366, 1.207, 2.742, 4.316, 8.266) / standardisation,
  
  # PCMCI COR (VAR heavy), n=500 and n=5000 respectively
  c(0.15, 1.61, 5.84, 12.26, 54.5) / standardisation,
  c(0.18, 2.63, 8.7, 17.71, 66.75) / standardisation,
  # PCMCI COR (GARCH heavy), n=500 and n=5000 respectively
  c(0.74, 2.9, 9.28, 18.7, 73.60) / standardisation,
  c(0.66, 3.21, 9.54, 18.62, 72) / standardisation,
  # PCMCI COR (VAR non-heavy), n=500 and n=5000 respectively
  c(0.12, 1.26, 5.42, 12.05, 52.87) / standardisation,
  c(0.12, 1.25, 5.05, 11.86, 48) / standardisation,
  # PCMCI COR (GARCH non-heavy), n=500 and n=5000 respectively
  c(0.84, 3.52, 9.93, 19.89, 74.83) / standardisation,
  c(0.8, 3.44, 9.92, 20.27, 76.2) / standardisation,
  
  # PCMCI GPDC (VAR heavy), n=500 and n=5000 respectively
  c(0.39, 2.91, 9.43, NA, NA) / (sequence_for_GPDC * (sequence_for_GPDC - 1)),
  c(NA, NA, NA, NA, NA),
  # PCMCI GPDC (GARCH heavy), n=500 and n=5000 respectively
  c(0.57, 2.74, 7.87, NA, NA) / (sequence_for_GPDC * (sequence_for_GPDC - 1)),
  c(NA, NA, NA, NA, NA),
  # PCMCI GPDC (VAR non-heavy), n=500 and n=5000 respectively
  c(0.11, 0.77, 3.09, NA, NA) / (sequence_for_GPDC * (sequence_for_GPDC - 1)),
  c(NA, NA, NA, NA, NA),
  # PCMCI GPDC (GARCH non-heavy), n=500 and n=5000 respectively
  c(0.58, 2.74, 7.87, NA, NA) / (sequence_for_GPDC * (sequence_for_GPDC - 1)),
  c(NA, NA, NA, NA, NA),
  # Random data
  rnorm(5, 0.5, 0.025), #0.025 is approximatelly the variance of the random graph error.
  rnorm(5, 0.5, 0.025),
  rnorm(5, 0.5, 0.025),
  rnorm(5, 0.5, 0.025),
  rnorm(5, 0.5, 0.025),
  rnorm(5, 0.5, 0.025),
  rnorm(5, 0.5, 0.025),
  rnorm(5, 0.5, 0.025)
)



library(ggplot2)
library(dplyr)
library(tidyr)

data <- data.frame()

# Fill the data frame
for (i in seq_along(data_vectors)) {
  # Determine method, setting_label, and n_label
  method_index <- (i - 1) %/% 8 + 1
  setting_index <- ((i-1)%/%2) %% 4 + 1
  n_label_index <- ifelse(i %% 2 == 1, 1, 2)
  
  # Create a temporary data frame for this vector
  temp_df <- expand.grid(
    p = sequence_of_p,
    method =c("Our method", "PCMCI cor", "PCMCI gpdc", "Random")[method_index],
    n_label = c("n = 500", "n = 5000")[n_label_index],
    setting_label = c("VAR heavy-tailed", "GARCH heavy-tailed", "VAR Gaussian", "GARCH Gaussian")[setting_index]
  )
  
  # Assign the values from data_vectors to the temporary data frame
  temp_df$structInterv_dist <- data_vectors[[i]]
  
  # Append to the main data frame
  data <- rbind(data, temp_df)
}

# Plotting
ggplot(data, aes(x = p, y = structInterv_dist, color = method)) +
  geom_line(size = 0.7, linetype = "solid") +
  geom_point(size = 2, shape = 16) + # Use shape = 16 for solid circles
  ylab("Average error") +
  xlab('Number of variables') +
  facet_grid(setting_label ~ n_label) + # Transposed facets
  scale_colour_manual(values = c('black', "#EE6677", "#228833", "#1f77b4")) +
  theme(strip.text = element_text(size = 10), # Smaller font size for facet labels
        axis.text = element_text(size = 12), # Smaller font size for axis text
        axis.title = element_text(size = 12), # Smaller font size for axis titles
        legend.text = element_text(size = 10), # Smaller font size for legend text
        legend.title = element_blank(), # Hide legend title
        legend.position = "right", # Move legend to the right
        panel.grid.major = element_line(size = 0.5, linetype = "dashed", color = "grey80"),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1) + # Ensure each facet is squared
  scale_x_continuous(breaks = seq(0, 20, by = 5)) + # Regular breaks
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))
#7x10 in export pdf












































#Time measurement - how long does it take to estimate one graph? 
x=generate_data(n=500, p=20, structure = structure, heavy_tailed=heavy_tailed)
start_time <- Sys.time()

G=Extreme_causality_full_graph_estimate(x$data)

end_time <- Sys.time()
print(end_time - start_time)










