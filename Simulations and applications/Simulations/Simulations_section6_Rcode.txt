#Simulations from section 6 where we compared our method with other state-of-the-art methods on VAR and GARCH models
#Python code can be found in the other file that deals with other methods from the Tigramite package

#generate_data is the function that generated the data with sample size n, number of variables p (p=m from the manuscript) and structure='VAR' or 'GARCH' and heavy_tailed=TRUE or FALSE. 
#For generating random graph we use the function erdos.renyi.game(p, 1/p, directed = TRUE,loops = FALSE) from igraph library
#Lines 100-150 are the main part. We generate data with random graph + estimates the graph using our method + compute the distance between true graph and estimated graph + repeat 100 times and return the mean of the distances
#We maually rewrote the values for each combination of 'structure' and 'heavy_tailed' into excel (easier since we dealt with python + R results)
#The table with final results can be found here at lines 150+


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
















standardisation = sequence_of_p*(sequence_of_p-1) #Because there are p*(p-1) arrows to estimate

OUR_n_500_VAR_heavy = c(0.0198, 0.1485, 0.732, 1.5346, 6.5754257)/standardisation
OUR_n_5000_VAR_heavy = c(0.000, 0.0693, 0.1287, 0.4554, 1.485148)/standardisation

OUR_n_500_GARCH_heavy = c(0.079, 0.316,1.336, 2.940, 13.019 )/standardisation
OUR_n_5000_GARCH_heavy = c(0.009, 0.069, 0.376, 1.128, 5.546)/standardisation

OUR_n_500_VAR_nonheavy = c(0.69306, 1.93069, 4.5247, 7.5643, 23.831)/standardisation
OUR_n_5000_VAR_nonheavy = c(0.48, 1.31, 2.95, 4.32, 9.81)/standardisation

OUR_n_500_GARCH_nonheavy = c(0.603, 1.831, 4.287, 7.079, 20.336)/standardisation
OUR_n_5000_GARCH_nonheavy = c(0.366, 1.207, 2.742, 4.316, 8.266)/standardisation




PCMCI_COR_n_500_VAR_heavy = c( 0.15,  1.61,  5.84, 12.26, 54.5 )/standardisation
PCMCI_COR_n_5000_VAR_heavy = c( 0.18,  2.63,  8.7 , 17.71, 66.75)/standardisation

PCMCI_COR_n_500_GARCH_heavy = c(0.74,  2.9 ,  9.28, 18.7, 73.60)/standardisation
PCMCI_COR_n_5000_GARCH_heavy = c(0.66, 3.21, 9.54, 18.62, ???)/standardisation

PCMCI_COR_n_500_VAR_nonheavy = c(0.12, 1.26,  5.42,  12.05, 52.87)/standardisation
PCMCI_COR_n_5000_VAR_nonheavy = c(0.12, 1.25, 5.05, 11.86, ???)/standardisation

PCMCI_COR_n_500_GARCH_nonheavy = c(0.94,  3.52,  9.93, 19.89, 74.83)/standardisation
PCMCI_COR_n_5000_GARCH_nonheavy = c( 0.8 ,  3.44,  9.92, 20.27, 76.2)/standardisation




PCMCI_GPDC_n_500_VAR_heavy = c(0.39, 3.23, 9.43, ???, ???)/standardisation
PCMCI_GPDC_n_5000_VAR_heavy = c(???, ???, ???, ???, ???)/standardisation

PCMCI_GPDC_n_500_GARCH_heavy = c(0.57, 2.94, 7.87, ???, ???)/standardisation
PCMCI_GPDC_n_5000_GARCH_heavy = c(???, ???, ???, ???, ???)/standardisation

PCMCI_GPDC_n_500_VAR_nonheavy = c(0.11, 0.77, 3.09, ???, ???)/standardisation
PCMCI_GPDC_n_5000_VAR_nonheavy =c(???, ???, ???, ???, ???)/standardisation

PCMCI_GPDC_n_500_GARCH_nonheavy = c(0.58, 3.79, 7.87, ???, ???)/standardisation
PCMCI_GPDC_n_5000_GARCH_nonheavy = c(???, ???, ???, ???, ???)/standardisation





















































