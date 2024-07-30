# Welcome to the Granger Causality in Extremes Package!
# This package contains the two main R functions from our framework.
# For any questions or an error, contact juraj.bodik@unil.ch.

######################### Extreme_causality_test ###############################################################
# This function estimates if tails/extremes of X cause Y given Z

# Function Inputs:
# 'x'                     A numeric vector representing the first time series (potential cause).
# 'y'                     A numeric vector representing the second time series (potential effect).
# 'z'                     A data.frame of potential confounders. Set to NULL if there are no confounders.
# 'lag_future'            The time delay for the effect from x to y. This is the coefficient 'p' in Appendix A of the manuscript.
# 'p_value_computation'   Set to TRUE to compute the p-value for the hypothesis H_0: X does not cause Y in extremes given Z. If p_value < 0.05, we conclude that X causes Y given Z.
# 'bootstrap_repetitions' The number of bootstrap repetitions for p-value computation. More repetitions yield more precise p-values but require longer computation time.
# 'both_tails'            Set to TRUE to consider both large and extremely negative values. For example, in GARCH models, both tails are of interest, while in VAR models, only large values might be relevant.
# 'nu_x'                  The coefficient tau_X or k_n in the manuscript, defined as "k=floor(n^{nu_x})". If strong hidden confounding is expected, set nu_x to 0.4 or 0.5.
# 'q_y'                   The coefficient tau_y = q_y * n, describing the conditioning on Y_t. For large auto-correlation in Y, set q_y to 0.1 or less. Note that in the manuscript, q_y is defined as 1 - q_y.
# 'q_z'                   The coefficient tau_z = q_z * n, describing the conditioning on Z_t. This is irrelevant if Z is NULL. For strong confounding effects, set q_z to 0.2 or 0.3. Note that in the manuscript, q_z is defined as 1 - q_z.
# 'lag_past'              The lag from Z to (X, Y). If the common cause has different lags to X and Y, it may cause spurious causality between X and Y. Ensure lag_past is larger than this lag.
# 'choice_of_F'           Choice of F in the coefficient. Leave default unless you want to reproduce the results from the manuscript

# Function Outputs:
# 'output'                Either 'Evidence of causality' or 'No causality' based on Algorithm 1 from the manuscript.
# 'p_value_tail'          This is not shown if p_value_computation==FALSE. Rejection indicates evidence of causality in tail.  It corresponds to the p-value for the hypothesis H_0: X does not cause Y in tail given Z, based on bootstrapping. Often p_value==1 which means that CTC<baseline
# 'p_value_extreme'       This is not shown if p_value_computation==FALSE. Rejection indicates evidence of causality in extremes. It corresponds to the p-value for the hypothesis H_0: hat{Gamma}_{X-->Y | Z} < (1 + 3 * hat{Gamma}^{baseline}_{X-->Y | Z}) / 4. 
# 'CTC'                   The coefficient hat{Gamma}_{X-->Y | Z}.
# 'baseline'              The baseline coefficient hat{Gamma}^{baseline}_{X-->Y | Z}.

######################### Extreme_causality_full_graph_estimate ###############################################################
# This function estimates the causal graph (path diagram) between a set of time series. 

# Function Inputs:
# 'w'                     A data.frame of all time series, which should be numeric and of the same length.
# 'lag_future'            Same as in Extreme_causality_test.
# 'both_tails'            Same as in Extreme_causality_test.
# 'nu_x'                  Same as in Extreme_causality_test.
# 'q_y'                   Same as in Extreme_causality_test.
# 'q_z'                   Same as in Extreme_causality_test.
# 'lag_past'              Same as in Extreme_causality_test.

# Function Outputs:
# 'G$G'                   A graph defined by its edges. Each row corresponds to an edge from the first column pointing to the second column. Use graph <- graph_from_edgelist(G$G) from the igraph library to obtain the graph environment.
# 'G$weights'             Weights corresponding to each edge, representing how close the coefficient hat{Gamma}_{X-->Y | Z} is to 1. If hat{Gamma}_{X-->Y | Z} = 1, the weight is 1. The weight is 0 if hat{Gamma}_{X-->Y | Z} = (1 + hat{Gamma}^{baseline}_{X-->Y | Z}) / 2.

# Load required libraries
library(igraph)   # For visualizing the final graph estimates
library(EnvStats) # Or any other package that can generate Pareto noise

## Example: Generating a 4-dimensional VAR time series with lag=2
# n = 5000
# epsilon_x = rpareto(n, 1, 1)
# epsilon_y = rpareto(n, 1, 1)
# epsilon_z1 = rpareto(n, 1, 1)
# epsilon_z2 = rpareto(n, 1, 1)
# x = rep(0, n); y = rep(0, n); z1 = rep(0, n); z2 = rep(0, n)
# for (i in 3:n) {
#   z1[i] = 0.5 * z1[i - 1] + epsilon_z1[i]
#   z2[i] = 0.5 * z2[i - 1] + epsilon_z2[i]
#   x[i] = 0.5 * x[i - 1] + 0.5 * z1[i - 2] + 0.5 * z2[i - 1] + epsilon_x[i]
#   y[i] = 0.5 * y[i - 1] + 0.5 * z1[i - 1] + 0.5 * z2[i - 2] + 0.2 * x[i - 1] + epsilon_y[i]
# }
# z = data.frame(z1, z2)
# w = data.frame(z1, z2, x, y)

## Running the extreme causality tests
# Extreme_causality_test(x, y, z, lag_future = 2, p_value_computation = FALSE) 
# Extreme_causality_test(y, x, z, lag_future = 2, p_value_computation = FALSE)

## Estimating the full causality graph
# G = Extreme_causality_full_graph_estimate(w, lag_future = 2) #Try it out also with lag = 1. You will see that the lagged edges dissapear

## Visualizing the graph using igraph
# graph <- graph_from_edgelist(G$G)
# V(graph)$name <- names(w)
# plot(graph, layout = layout_nicely(graph), vertex.label = V(graph)$name)





Extreme_causality_test = function(x, y, z=NULL, lag_future=1, lag_past=0, nu_x = 0.3, q_y = 0.2, q_z = 0.1, both_tails = TRUE, instant=FALSE, p_value_computation = FALSE, bootstrap_repetitions=50, choice_of_F = 0.5){
  
  n = length(x)
  z=data.frame(z)
  d = ncol(z) ;
  if (d >= 2) {q_z <- (2/d)*q_z};
  tau_y = q_y * n; 
  tau_z = q_z * n
   if(both_tails==TRUE){ 
    if(!all(x>=0))  x=abs(x-median(x)); 
    if(!all(y>=0))  y =abs(y-median(y)); 
    if(d>=1 & !all(z>=0))  z = abs(z - apply(z, 2, median))
    }
  
  

  F_u = function(imput, output){ifelse(output < quantile(imput, choice_of_F), 0, (ecdf(imput)(output)))} #This is the F^{truc}_Y(t) function that we opt for
  
  
  CTC_baseline = function(x, y, z=NULL, lag_future=lag_future, lag_past = lag_past, tau_y = tau_y, tau_z = tau_z, instant=instant){
    
    run_future_max <- function(x, k=3, instant=TRUE){ #instant= do we want to consider Y_0?
      if (instant==TRUE) {q=0}else{q=1}
      
      n=length(x);  y=c()  
      for (i in 1:n) {y=rbind(y, max(x[(min(n,i+q)):(min(n,i+k))]))}
      
      return(y)
    }
    
    x_to_y_masking= c()
    for (i in 1:(lag_past+1)) { yy = c( rep(0, i-1), y)
    x_to_y_masking= c(x_to_y_masking,  which(yy>sort(y)[(n-tau_y)])) 
    
    if(d>0){ for (j in 1:d){  
      zz = c( rep(0, i-1), z[,j])
      x_to_y_masking= c(x_to_y_masking,  which(zz>sort(z[,j])[(n-tau_z)]))}
    }
    x_to_y_masking = unique(x_to_y_masking)
    }
    
    future_y = run_future_max(y, lag_future, instant = instant)
    #baseline1 =mean(F_u(y[-x_to_y_masking], future_y[-x_to_y_masking]))
    baseline2=mean(F_u(y, future_y[-x_to_y_masking]))
    
    return( baseline2  )
  }
  
  
  CTC_masked = function(x, y, z=NULL, lag_future, lag_past, tau_y, tau_z, instant){
    
    run_future_max <- function(x, k=3, instant=TRUE){ #instant= do we want to consider Y_0?
      if (instant==TRUE) {q=0}else{q=1}
      
      n=length(x);  y=c()  
      for (i in 1:n) {y=rbind(y, max(x[(min(n,i+q)):(min(n,i+k))]))}
      
      return(y)
    }
    
    
    x_to_y_masking= c()
    for (i in 1:(lag_past+1)) { yy = c( rep(0, i-1), y)
    x_to_y_masking= c(x_to_y_masking,  which(yy>sort(y)[(n-tau_y)])) 
    
    if(d>0){ for (j in 1:d){  
      zz = c( rep(0, i-1), z[,j])
      x_to_y_masking= c(x_to_y_masking,  which(zz>sort(z[,j])[(n-tau_z)]))}
    }
    x_to_y_masking = unique(x_to_y_masking)
    }
    
    n=length(x); k = round((n-length(x_to_y_masking))^(nu_x) )
    
    new_x = x[-x_to_y_masking]
    
    future_y = run_future_max(y, lag_future, instant = instant)
    
    future_y = future_y[-x_to_y_masking]
    
    top_x=which(new_x>=sort(new_x)[length(new_x)-k+1])
    
    
    #return( mean(F_u(y[-x_to_y_masking], future_y[top_x]))  )
    return( mean(F_u(y, future_y[top_x]))  )
    
  }
  
  
  
  
  switcher_for_test <- function(x, number_of_blocks=15){ #resampling block-wise
    
    n=nrow(x)
    m=n%/%number_of_blocks #length of one block
    
    y=c()
    for (i in 1:number_of_blocks) {#choose one random block with the beginning uniformly chosen from 1:(n-m)
      kocka=sample(1:(n-m), 1)
      for (j in 1:m) {
        y=rbind(y,x[kocka+j,])  #Add this block to our resampled series
      }
      
    }
    #what to do with the ending if it is not divisible? we just add one random block with smaller length to obtain times eris with length $n$ again
    if (ncol(x)==1) {  #Code for one-dimensional time series
      if(n%%number_of_blocks!=0){y=c(y,x[((number_of_blocks*m+1):n),])}} 
    else{
      if(n%%number_of_blocks!=0){
        k=n-number_of_blocks*m
        kocka=sample(1:(n-k), 1)
        
        for (j in 1:k) {
          y=rbind(y,x[kocka+j,])  
        }}
    }
    return(data.frame(y))
  }
  
  
  
  baseline = CTC_baseline(x, y, z=z, lag_future=lag_future, lag_past = lag_past, tau_y = tau_y, tau_z = tau_z, instant=instant)
  CTC = CTC_masked(x,y,z,lag_future=lag_future,lag_past = lag_past, tau_y = tau_y, tau_z = tau_z, instant=instant)
  
  
  if(p_value_computation == FALSE){
    if(CTC<=(1+baseline)/2) return(data.frame(output = 'No causality', CTC=CTC , baseline = baseline))
    if(CTC>(1+baseline)/2) return(data.frame(output = 'Evidence of causality', CTC=CTC , baseline = baseline))
  }
  
  
  
  if(CTC<=baseline){result = 0; result2=1} 
  if(CTC>baseline){
    result = c()
    result2 = c()
    
    for (i in 1:bootstrap_repetitions) {
      if(d==0) original_data = data.frame(x,y) else original_data = data.frame(x,y, z)
      new_time_series=switcher_for_test(original_data, number_of_blocks=floor(sqrt(length(x))))
      tilde_x = new_time_series[,1] ;  tilde_y =new_time_series[,2] ; if(d==0){tilde_z = NULL};  if(d>0){tilde_z =  data.frame(new_time_series[,3:(d+2)])};  
      result  = c(result, CTC_masked(tilde_x, tilde_y, tilde_z, lag_future=lag_future,lag_past = lag_past, tau_y = tau_y, tau_z = tau_z, instant=instant ) )
      result2 = c(result2, CTC_baseline(tilde_x, tilde_y, tilde_z, lag_future=lag_future,lag_past = lag_past, tau_y = tau_y, tau_z = tau_z, instant=instant ) )
    }
  }
  
  
  output = 'No causality';  if(CTC>(1+baseline)/2)   output = 'Evidence of causality'; 
  return( data.frame( output=output, p_value_tail = mean(result <=  result2) , p_value_extreme = mean(result <=  (1+3*result2)/4), CTC=CTC , baseline = baseline))
}



Extreme_causality_full_graph_estimate = function(w, lag_future=1, lag_past=0, nu_x = 0.3,  q_y = 0.2, q_z = 0.1, instant=FALSE, both_tails = TRUE){
  
  m=ncol(w)
  
  find_parents = function(G, vertex){
    result = c()
    for (i in 1:nrow(G)) { 
      if(G[i,2]==vertex){result = c(result, G[i,1])}
    }
    return(result)
  }
  
  compute_edge_weights = function(CTC, baseline){ (CTC -((1+baseline)/2))/(1-((1+baseline)/2))  }
  
  #Step 1: Pairwise    
  G = c()  
  for (i in 1:m) {
    for (j in (1:m)[-i]) {
      x=w[,i];y=w[,j]
      CTC=Extreme_causality_test(x,y,z=NULL,  nu_x = nu_x,   q_y =q_y, q_z=q_z, lag_future = lag_future, lag_past = lag_past,both_tails=both_tails)
      if(CTC$output =='Evidence of causality') G=rbind(G, c(i,j))    
    }  
  }
  
  #Step 2: Multivariate  
  if( all(G == FALSE)  ){return('Result: Empty graph')}
  if( !all(G == FALSE)  ){ #if G is non-empty
    indexes_to_erase = .Machine$integer.max
    edges_weights = c()
    for (i in 1:nrow(G)) {
      x = w[,G[i,1]]
      y = w[,G[i,2]]
      
      z_indexes=intersect(find_parents(G, G[i,1]), find_parents(G, G[i,2]))  
      if(all(z_indexes == FALSE)) {z=NULL}  
      if(!all(z_indexes == FALSE)) {z=data.frame(w[,z_indexes])}  
      
      CTC=Extreme_causality_test(x,y,z=z,  nu_x = nu_x,   q_y =q_y, q_z=q_z, lag_future = lag_future, lag_past = lag_past, both_tails=both_tails)
      if(CTC$output == 'No causality'){ indexes_to_erase = c(indexes_to_erase, i)}else{ edges_weights = c(edges_weights, compute_edge_weights(CTC$CTC, CTC$baseline))}
    }}
  
  return(list(G = G[-indexes_to_erase,], weights = edges_weights))  
}






