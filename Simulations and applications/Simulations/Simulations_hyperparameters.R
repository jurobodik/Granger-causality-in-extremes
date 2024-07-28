#This code serves for reproducibility purposes concerning the simulations about the choice of hyper-parameters
#You need to first upload the function 'Extreme_causality_test' from the other file
#Lines 1-200 corresponds to the choice of q_F; that is, optimal choice of F
#Lines 200-400 corresponds to the choice of tau_X
#Lines 400-600 corresponds to the choice of tau_Y
#Lines 600-700 corresponds to the choice of tau_Z
#Lines 700-1000 correspond to the final plots using ggplot
library(EnvStats) #To generate Pareto noise
set.seed(2)

#Simulations for q_F
#Run lines 11-70 to obtain the first graph. Lines 70-130 to obtain the third graph. To obtain the second or the fourth graph, uncomment the lines that are marked with `` Un-comment this line ''

#VAR model
result = c()
sekv = seq(0.01, 0.61, by=0.1) #This is the sequence for alpha in the Gaussian model
#sekv = seq(0.01, 0.21, by=0.05) #sequence for Pareto noise. Un-comment this line for the heavy-tailed case result

for (q in c(0, 0.3,  0.5,  0.7)) {
  result2=c()
  cat('Time remaining: q = ',  q, '\n')  
  
  for (alpha in sekv) {
    result1=c()
    
    for (k in 100:1) {
      n=500
      
      epsilon_x=rnorm(n, 0,1)  
      epsilon_y=rnorm(n, 0,1)  
      epsilon_z=rnorm(n, 0,1)  
      
#      epsilon_x=rpareto(n, 1, 1)  # Un-comment this line for the heavy-tailed case result
#      epsilon_y=rpareto(n, 1, 1)  # Un-comment this line for the heavy-tailed case result
#      epsilon_z=rpareto(n, 1, 1)  # Un-comment this line for the heavy-tailed case result
      
      x=rep(0, n);y=rep(0, n);z=rep(0, n)
      
      for (i in 3:n) {
        z[i]=0.5*z[i-1]  + epsilon_z[i]
        x[i]=0.5*x[i-1]  + 0.5*z[i-1] + epsilon_x[i]
        y[i]=0.5*y[i-1]  + 0.5*z[i-1] + alpha*x[i-1] + epsilon_y[i]
      }
      
      
      CTC1 = Extreme_causality_test(x,y,z,    q_y =0.1, q_z=0.1, lag_future = 1, lag_past = 0, p_value_computation = FALSE, choice_of_F=q, both_tails = FALSE)
      CTC2 = Extreme_causality_test(y,x,z,    q_y =0.1, q_z=0.1, lag_future = 1, lag_past = 0, p_value_computation = FALSE, choice_of_F=q, both_tails = FALSE)
      result1 = c(result1, sum(CTC1$output=='Evidence of causality', CTC2$output=='No causality'))
    }
    result2 = c(result2, sum(result1)) 
  }
  result = rbind(result, result2)
}


result = result/2

# Plot the lines
plot(result[1,]~sekv, type = 'l', col=1, lwd=3, ylim=c(min(result), max(result)), xlab='alpha_x', ylab="Performance", main = 'VAR non-heavy-tailed case')
for (i in 2:nrow(result)) {
  lines(result[i,]~sekv, type = 'l', col=i, lwd=3, lty=i-1)
}


# Create legend
legend("bottomright", 
       legend=c("q_F=0", "q_F=0.3", "q_F=0.5",  "q_F=0.7"),
       col=c(1, 2, 3, 4),
       lty=c(1, 1, 2, 3), lwd=3)





#GARCH model
result = c()
sekv = seq(0.001, 0.011, by=0.002)  #This is the sequence for alpha
#sekv = seq(0.01, 5.01, by=1)  #sequence for normal noise.  Un-comment this line for the non-heavy-tailed case result


for (q in c(0, 0.3,  0.5,  0.7)) {
  result2=c()
  cat('Time remaining: ',  q, '\n')  
  
  for (alpha in sekv) {
    result1=c()
    
    
    for (k in 100:1) {
      n=500
      
      epsilon_x=rcauchy(n, 0,1)
      epsilon_y=rcauchy(n, 0,1)
      epsilon_z=rcauchy(n, 0,1)
      
#      epsilon_x=rnorm(n, 0,1)  #Un-comment this line for the non-heavy-tailed case result
#      epsilon_y=rnorm(n, 0,1)  #Un-comment this line for the non-heavy-tailed case result
#      epsilon_z=rnorm(n, 0,1)  #Un-comment this line for the non-heavy-tailed case result 
      
      x=rep(0, n);y=rep(0, n);z=rep(0, n)
      
      for (i in 3:n) {
        z[i]=((0.1+ 0.1*z[i-1]^2)^(0.5))*epsilon_z[i]
        x[i]=((0.1+ 0.1*x[i-1]^2 +  0.1*z[i-1]^2)^(0.5))*epsilon_x[i]
        y[i]=((0.1+ 0.1*y[i-1]^2 +  0.1*z[i-1]^2 + alpha*x[i-1]^2 )^(0.5))*epsilon_y[i]
      }
      
      
      CTC1 = Extreme_causality_test(x,y,z,    q_y =0.1, q_z=0.1, lag_future = 1, lag_past = 0, p_value_computation = FALSE, choice_of_F=q, both_tails = TRUE)
      CTC2 = Extreme_causality_test(y,x,z,    q_y =0.1, q_z=0.1, lag_future = 1, lag_past = 0, p_value_computation = FALSE, choice_of_F=q, both_tails = TRUE)
      result1 = c(result1, sum(CTC1$output=='Evidence of causality', CTC2$output=='No causality'))
    }
    result2 = c(result2, sum(result1)) 
  }
  result = rbind(result, result2)
}

result = result/2

# Plot the lines
plot(result[1,]~sekv, type = 'l', col=1, lwd=3, ylim=c(min(result), max(result)), xlab='alpha_x', ylab="Performance", main = 'GARCH heavy-tailed case')
for (i in 2:nrow(result)) {
  lines(result[i,]~sekv, type = 'l', col=i, lwd=3, lty=i-1)
}


# Create legend
legend("bottomright", 
       legend=c("q_F=0", "q_F=0.3", "q_F=0.5",  "q_F=0.7"),
       col=c(1, 2, 3, 4),
       lty=c(1,1, 2, 3), lwd=3)





































































#Simulations for k_n
#You can simply run all lines from here to line 400 to obtain the decired graphs
#This function generate the considered 4 models
generate_series = function(n=500, heavy_tailed = TRUE, VAR_or_GARCH = 'VAR'){
  
  x=rep(0, n);y=rep(0, n);z=rep(0, n)
  
  if(heavy_tailed == FALSE & VAR_or_GARCH == 'VAR'){
    epsilon_x=rnorm(n, 0,1) 
    epsilon_y=rnorm(n, 0,1) 
    epsilon_z=rnorm(n, 0,1) 
    
    
    for (i in 3:n) {
      z[i]=0.5*z[i-1]  + epsilon_z[i]
      x[i]=0.5*x[i-1]  + 0.5*z[i-1] + epsilon_x[i]
      y[i]=0.5*y[i-1]  + 0.5*z[i-1] + 0.5*x[i-1] + epsilon_y[i]
    }
  }
  
  if(heavy_tailed == TRUE & VAR_or_GARCH == 'VAR'){
    epsilon_x=rpareto(n, 1,1) 
    epsilon_y=rpareto(n, 1,1) 
    epsilon_z=rpareto(n, 1,1) 
    
    
    for (i in 3:n) {
      z[i]=0.5*z[i-1]  + epsilon_z[i]
      x[i]=0.5*x[i-1]  + 0.5*z[i-1] + epsilon_x[i]
      y[i]=0.5*y[i-1]  + 0.5*z[i-1] + 0.1*x[i-1] + epsilon_y[i]
    }
  }
  
  
  
  if(heavy_tailed == FALSE & VAR_or_GARCH == 'GARCH'){
    epsilon_x=5*rnorm(n, 0,1) 
    epsilon_y=rnorm(n, 0,1) 
    epsilon_z=rnorm(n, 0,1) 
    
    
    for (i in 3:n) {
      z[i]=((0.1+ 0.1*z[i-1]^2)^(0.5))*epsilon_z[i]
      x[i]=((0.1+ 0.1*x[i-1]^2 +  0.1*z[i-1]^2)^(0.5))*epsilon_x[i]
      y[i]=((0.1+ 0.1*y[i-1]^2 +  0.1*z[i-1]^2 + 5*x[i-1]^2 )^(0.5))*epsilon_y[i]
    }
    x=abs(x); y=abs(y); z=abs(z)
  }
  
  if(heavy_tailed == TRUE & VAR_or_GARCH == 'GARCH'){
    epsilon_x=rcauchy(n, 0,1) 
    epsilon_y=rcauchy(n, 0,1) 
    epsilon_z=rcauchy(n, 0,1) 
    
    for (i in 3:n) {
      z[i]=((0.1+ 0.1*z[i-1]^2)^(0.5))*epsilon_z[i]
      x[i]=((0.1+ 0.1*x[i-1]^2 +  0.1*z[i-1]^2)^(0.5))*epsilon_x[i]
      y[i]=((0.1+ 0.1*y[i-1]^2 +  0.1*z[i-1]^2 + 0.5*x[i-1]^2 )^(0.5))*epsilon_y[i]
    }
    x=abs(x); y=abs(y); z=abs(z)
  }
  
  return(data.frame(x,y,z))
}



#No hidden confounder case
result = c()
sekv_n = c(200, 400, 600)
sekv_tau = c(0.2, 0.3, 0.4, 0.5, 0.6)



for (nu_x in sekv_tau) {
  result2=c()
  cat('Time remaining: ',  nu_x, '\n')  
  
  for (n in sekv_n) {
    result1=c()
    alpha = 0.5
    for (heavy_tailed in c(FALSE, TRUE)) {
      for (VAR_or_GARCH in c('VAR', 'GARCH')) {  
        for (k in 100:1) {
          
          data = generate_series(n=n,  heavy_tailed = heavy_tailed, VAR_or_GARCH = VAR_or_GARCH)
          x=data$x; y=data$y; z=data$z  
          
          CTC1 = Extreme_causality_test(x,y,z,  q_y =0.2, q_z=0.1,  nu_x = nu_x, both_tails = FALSE)
          CTC2 = Extreme_causality_test(y,x,z,  q_y =0.2, q_z=0.1,  nu_x = nu_x, both_tails = FALSE)
          
          result1 = c(result1, sum(CTC1$output=='Evidence of causality', CTC2$output=='No causality'))
        }}}
    result2 = c(result2, sum(result1)) 
  }
  result = rbind(result, result2)
}

result = result / 8
par(mfrow = c(1, 1))

# Plot the lines
plot(result[,1]~sekv_tau, type = 'l', col=1, lwd=3, ylim=c(min(result), max(result)), xlab='k_n', ylab="Performance", main = 'Hidden confounder case')
for (i in 2:ncol(result)) {
  lines(result[,i]~sekv_tau, type = 'l', col=i, lwd=3, lty=i-1)
}

# Create legend
legend("bottomright", 
       legend=c("n=200", "n=400","n=600"),
       col=c(1, 2, 3),
       lty=c(1, 1, 2), lwd=3)


result_no_hidden = result




#Hidden confounder case

result = c()
sekv_n = c(200, 400, 600)
sekv_tau = c(0.2, 0.3, 0.4, 0.5, 0.6)

for (nu_x in sekv_tau) {
  result2=c()
  cat('Time remaining: ',  nu_x, '\n')  
  
  for (n in sekv_n) {
    result1=c()
    for (heavy_tailed in c(FALSE, TRUE)) {
      for (VAR_or_GARCH in c('VAR', 'GARCH')) {  
        for (k in 200:1) {
          
          data = generate_series(n=n,  heavy_tailed = heavy_tailed, VAR_or_GARCH = VAR_or_GARCH)
          x=data$x; y=data$y; 
          
          CTC1 = Extreme_causality_test(x,y,z=NULL,    q_y =0.1, nu_x = nu_x, both_tails = FALSE)
          CTC2 = Extreme_causality_test(y,x,z=NULL,    q_y =0.1, nu_x = nu_x, both_tails = FALSE)
          
          result1 = c(result1, sum(CTC1$output=='Evidence of causality', CTC2$output=='No causality'))
        }}}
    result2 = c(result2, sum(result1)) 
  }
  result = rbind(result, result2)
}

result = result / 16
result_hidden = result

par(mfrow = c(1, 1))

# Plot the lines
plot(result[,1]~sekv_tau, type = 'l', col=1, lwd=3, ylim=c(min(result), max(result)), xlab='k_n', ylab="Performance", main = 'Hidden confounder case')
for (i in 2:ncol(result)) {
  lines(result[,i]~sekv_tau, type = 'l', col=i, lwd=3, lty=i-1)
}

# Create legend
legend("bottomright", 
       legend=c("n=200", "n=400","n=600"),
       col=c(1, 2, 3),
       lty=c(1, 1, 2), lwd=3)




par(mfrow = c(1, 2))

# Plot the lines
plot(result_no_hidden[,1]~sekv_tau, type = 'l', col=1, lwd=3, ylim=c(min(result_no_hidden), max(result_no_hidden)), xlab='k_n', ylab="Performance", main = 'No Hidden confounder')
for (i in 2:ncol(result_no_hidden)) {
  lines(result_no_hidden[,i]~sekv_tau, type = 'l', col=i, lwd=3, lty=i-1)
}

# Create legend
legend("bottomleft", 
       legend=c("n=200", "n=400","n=600"),
       col=c(1, 2, 3),
       lty=c(1, 1, 2), lwd=3)

# Plot the lines
plot(result_hidden[,1]~sekv_tau, type = 'l', col=1, lwd=3, ylim=c(min(result_hidden), max(result_hidden)), xlab='k_n', ylab="Performance", main = 'Hidden confounder')
for (i in 2:ncol(result_hidden)) {
  lines(result_hidden[,i]~sekv_tau, type = 'l', col=i, lwd=3, lty=i-1)
}

# Create legend
legend("topleft", 
       legend=c("n=200", "n=400","n=600"),
       col=c(1, 2, 3),
       lty=c(1, 1, 2), lwd=3)








#Graphs for optimal q_Y
#Note that q_Y = 1-q_Y from the manuscript. Slightly confusing
generate_series = function(n=500, alpha_y=0.5,  heavy_tailed = TRUE, VAR_or_GARCH = 'VAR'){
  
  x=rep(0, n);y=rep(0, n);z=rep(0, n)
  
  if(heavy_tailed == FALSE & VAR_or_GARCH == 'VAR'){
    epsilon_x=rnorm(n, 0,1) 
    epsilon_y=rnorm(n, 0,1) 
    epsilon_z=rnorm(n, 0,1) 
    
    
    for (i in 3:n) {
      z[i]=0.5*z[i-1]  + epsilon_z[i]
      x[i]=0.5*x[i-1]  + 0.5*z[i-1] + epsilon_x[i]
      y[i]=alpha_y*y[i-1]  + 0.5*z[i-1] + 0.5*x[i-1] + epsilon_y[i]
    }
  }
  
  if(heavy_tailed == TRUE & VAR_or_GARCH == 'VAR'){
    epsilon_x=rpareto(n, 1,1) 
    epsilon_y=rpareto(n, 1,1) 
    epsilon_z=rpareto(n, 1,1) 
    
    
    for (i in 3:n) {
      z[i]=0.5*z[i-1]  + epsilon_z[i]
      x[i]=0.5*x[i-1]  + 0.5*z[i-1] + epsilon_x[i]
      y[i]=alpha_y*y[i-1]  + 0.5*z[i-1] + 0.1*x[i-1] + epsilon_y[i]
    }
  }
  
  
  
  if(heavy_tailed == FALSE & VAR_or_GARCH == 'GARCH'){
    epsilon_x=rnorm(n, 0,1) 
    epsilon_y=rnorm(n, 0,1) 
    epsilon_z=rnorm(n, 0,1) 
    
    
    for (i in 3:n) {
      z[i]=((0.1+ 0.1*z[i-1]^2)^(0.5))*epsilon_z[i]
      x[i]=((0.1+ 0.1*x[i-1]^2 +  0.1*z[i-1]^2)^(0.5))*epsilon_x[i]
      y[i]=((alpha_y/5+ 0.1*y[i-1]^2 +  0.1*z[i-1]^2 + 10*x[i-1]^2 )^(0.5))*epsilon_y[i]
    }
    x=abs(x); y=abs(y); z=abs(z)
  }
  
  if(heavy_tailed == TRUE & VAR_or_GARCH == 'GARCH'){
    epsilon_x=rcauchy(n, 0,1) 
    epsilon_y=rcauchy(n, 0,1) 
    epsilon_z=rcauchy(n, 0,1) 
    
    for (i in 3:n) {
      z[i]=((0.1+ 0.1*z[i-1]^2)^(0.5))*epsilon_z[i]
      x[i]=((0.1+ 0.1*x[i-1]^2 +  0.1*z[i-1]^2)^(0.5))*epsilon_x[i]
      y[i]=((alpha_y/5+ 0.1*y[i-1]^2 +  0.1*z[i-1]^2 + 0.5*x[i-1]^2 )^(0.5))*epsilon_y[i]
    }
    x=abs(x); y=abs(y); z=abs(z)
  }
  
  return(data.frame(x,y,z))
}



result = c()
n=500
sekv_alpha_y = c(0.1, 0.3, 0.5, 0.7)
sekv_q_y = c(0.01, 0.1, 0.2, 0.3, 0.4)



for (q_y in sekv_q_y) {
  result2=c()
  cat('Time remaining: ',  q_y, '\n')  
  
  for (alpha_y in sekv_alpha_y) {
    result1=c()
    for (heavy_tailed in c(FALSE, TRUE)) {
      for (VAR_or_GARCH in c('VAR', 'GARCH')) {  
        for (k in 100:1) {
          
          data = generate_series(n=n,  alpha_y =  alpha_y, heavy_tailed = heavy_tailed, VAR_or_GARCH = VAR_or_GARCH)
          x=data$x; y=data$y; z=data$z
          
          CTC1 = Extreme_causality_test(x,y,z=z,     q_y =q_y,  nu_x = 0.3, both_tails = FALSE)
          CTC2 = Extreme_causality_test(y,x,z=z,     q_y =q_y,  nu_x = 0.3, both_tails = FALSE)
          
          result1 = c(result1, sum(CTC1$output=='Evidence of causality', CTC2$output=='No causality'))
        }}}
    result2 = c(result2, sum(result1)) 
  }
  result = rbind(result, result2)
}
result = result / 8
sekv_q_y = 1-sekv_q_y  #Note that q_Y = 1-q_Y from the manuscript. Slightly confusing





par(mfrow = c(1, 2))
plot(result[,1]~sekv_q_y, type = 'l', col=1, lwd=3, ylim=c(min(result), max(result)), xlab='q_Y', ylab="Performance", main = 'n=500')
for (i in 2:ncol(result)) {
  lines(result[,i]~sekv_q_y, type = 'l', col=i, lwd=3, lty=i-1)
}
legend("bottomleft", 
       legend=c("alpha_y=0.1", "alpha_y=0.3","alpha_y=0.5", 'alpha_y=0.7'),
       col=c(1, 2, 3, 4),
       lty=c(1, 1, 2, 3), lwd=3)


result500 = result









result = c()
n=10000
sekv_alpha_y = c(0.1, 0.3, 0.5, 0.7)
sekv_q_y = c(0.01, 0.1, 0.2, 0.3, 0.4)



for (q_y in sekv_q_y) {
  result2=c()
  cat('Time remaining: ',  q_y, '\n')  
  
  for (alpha_y in sekv_alpha_y) {
    result1=c()
    for (heavy_tailed in c(FALSE, TRUE)) {
      for (VAR_or_GARCH in c('VAR', 'GARCH')) {  
        for (k in 100:1) {
          
          data = generate_series(n=n,  alpha_y =  alpha_y, heavy_tailed = heavy_tailed, VAR_or_GARCH = VAR_or_GARCH)
          x=data$x; y=data$y; z=data$z
          
          CTC1 = Extreme_causality_test(x,y,z=z,     q_y =q_y,  nu_x = 0.3, both_tails = FALSE)
          CTC2 = Extreme_causality_test(y,x,z=z,     q_y =q_y,  nu_x = 0.3, both_tails = FALSE)
          
          result1 = c(result1, sum(CTC1$output=='Evidence of causality', CTC2$output=='No causality'))
        }}}
    result2 = c(result2, sum(result1)) 
  }
  result = rbind(result, result2)
}
result = result / 8
sekv_q_y = 1-sekv_q_y




# Plot the lines
plot(result[,1]~sekv_q_y, type = 'l', col=1, lwd=3, ylim=c(min(result), max(result)), xlab='q_Y', ylab="Performance", main = 'n=10000')
for (i in 2:ncol(result)) {
  lines(result[,i]~sekv_q_y, type = 'l', col=i, lwd=3, lty=i-1)
}

# Create legend
legend("bottomleft", 
       legend=c("alpha_y=0.1", "alpha_y=0.3","alpha_y=0.5", 'alpha_y=0.7'),
       col=c(1, 2, 3, 4),
       lty=c(1, 1, 2, 3), lwd=3)






























#Graphs for optimal q_Z
#Here there is only heavy-tailed VAR model. The other models had roughly the same results; you can just rewrite the data-generating process by the function 'generate_series'
result = c()
n=1000
sekv_alpha_z = c(0.1, 0.5, 1, 2)
sekv_q_z = c(0, 0.01, 0.05, 0.1,  0.2, 0.5)



for (alpha_z in sekv_alpha_z) {
  cat('Time remaining: ', alpha_z, '\n')
  result1=rep(0, length(sekv_q_z))
  for (k in 200:1) {
    
    
    epsilon_x=rpareto(n, 1,1)
    epsilon_y=rpareto(n, 1,1)
    epsilon_z=rpareto(n, 1,1)
    
    x=rep(0, n);y=rep(0, n);z=rep(0, n)
    
    for (i in 6:n) {
      z[i]=0.5*z[i-1]  + epsilon_z[i]
      x[i]=0.5*x[i-1]  + alpha_z*z[i-2] + epsilon_x[i]
      y[i]=0.5*y[i-1]  + alpha_z*z[i-1] + 0.1*x[i-1] + epsilon_y[i]
    }
    
    result2=c()
    for (q_z in sekv_q_z) {
      
      CTC1 = Extreme_causality_test(x,y,z=z,     q_z =q_z,  nu_x = 0.3, both_tails = FALSE)
      CTC2 = Extreme_causality_test(y,x,z=z,     q_z =q_z,  nu_x = 0.3, both_tails = FALSE)
      
      
      result2 = c(result2, sum(CTC1$output=='Evidence of causality', CTC2$output=='No causality'))
    }
    result1 = result1+result2
  }
  result = rbind(result, result1)
}

result = result / 4
sekv_q_z = 1-sekv_q_z#Note that q_z = 1-q_z from the manuscript. Slightly confusing
par(mfrow = c(1, 1))

# Plot the lines
plot(result[1,]~sekv_q_z, type = 'l', col=1, lwd=3, ylim=c(min(result), max(result)), xlab='q_Z', ylab="Performance", main = 'Choice of q_Z')
for (i in 2:nrow(result)) {
  lines(result[i,]~sekv_q_z, type = 'l', col=i, lwd=3, lty=i-1)
}

# Create legend
legend("bottomleft", 
       legend=c('alpha_z=0.1', 'alpha_z=0.5', 'alpha_z=1', 'alpha_z=2'),
       col=c(1, 2, 3, 4),
       lty=c(1, 1, 2, 3), lwd=3)












































#We recieved a complaint about the quality of the pictures, so we rewrote here all the plots in ggplot style
#It is super ugly how we rewrote the obtained final values by hand (using chatGPT), but it was the fastest way
library(ggplot2)
library(dplyr)
library(tidyr)

# PLOT 1
result_plot_1 <- rbind(c(100, 117, 121, 144, 160, 172, 180),   
                       c(101, 127, 133, 160, 172, 175, 182),   
                       c(99, 122, 145, 165, 172, 182, 186),
                       c(105, 109, 143, 157, 165, 173, 178))/2

result1 <- data.frame(
  alpha_X = seq(0.01, 0.61, by = 0.1),
  `q_F0` = result_plot_1[1,],
  `q_F0.3` = result_plot_1[2,],
  `q_F0.5` = result_plot_1[3,],
  `q_F0.7` = result_plot_1[4,]
)

# PLOT 2
result_plot_2 <- rbind(c(118, 169, 187, 195, 195),   
                       c(124, 178, 189, 197, 194),   
                       c(114, 185, 190, 198, 196),
                       c(105, 148, 174, 182, 190))/2

# Create data frames for each result set
result2 <- data.frame(
  alpha_X = seq(0.01, 0.21, by = 0.05),
  `q_F0` = result_plot_2[1,],
  `q_F0.3` = result_plot_2[2,],
  `q_F0.5` = result_plot_2[3,],
  `q_F0.7` = result_plot_2[4,]
)


result_plot_3 <- rbind(c(103, 141, 159, 160, 168, 163),   
                       c(98, 144, 153, 160, 169, 163),   
                       c(99, 150, 157, 163, 169, 169),
                       c(102, 148, 158, 160, 164, 168))/2

# Create data frames for each result set
result3 <- data.frame(
  alpha_X =  seq(0.01, 5.01, by = 1),
  `q_F0` = result_plot_3[1,],
  `q_F0.3` = result_plot_3[2,],
  `q_F0.5` = result_plot_3[3,],
  `q_F0.7` = result_plot_3[4,]
)

result_plot_4 <- rbind(c(189, 195, 196, 196, 192, 191),   
                       c(191, 190, 192, 192, 196, 195),   
                       c(188, 197, 197, 198, 196, 195),
                       c(183, 192, 195, 195, 197, 195))/2


# Create data frames for each result set
result4 <- data.frame(
  alpha_X = seq(0.01, 0.41, by = 0.075),
  `q_F0` = result_plot_4[1,],
  `q_F0.3` = result_plot_4[2,],
  `q_F0.5` = result_plot_4[3,],
  `q_F0.7` = result_plot_4[4,]
)

result_plot_5 <- cbind(c(81.50, 83.25, 84.75),   
                       c(83.25, 87.75, 90.75),   
                       c(81.25, 85.00, 89.00),
                       c(76.00, 80.25, 87.00),
                       c(70.70, 74.25, 81.75))

# Create data frames for each result set
result5 <- data.frame(
  k_n = seq(0.2, 0.6, by = 0.1),
  `n200` = result_plot_5[1,],
  `n400` = result_plot_5[2,],
  `n600` = result_plot_5[3,]
)

result_plot_6 <- cbind(c(78.00, 81.50, 81.75),   
                       c(81.25, 81.25, 84.50),   
                       c(82.00, 85.50, 86.75),
                       c(82.75, 86.75, 88.75),
                       c(79.00, 83.50, 87.00))

result6 <- data.frame(
  k_n = seq(0.2, 0.6, by = 0.1),
  `n200` = result_plot_6[1,],
  `n400` = result_plot_6[2,],
  `n600` = result_plot_6[3,]
)


result_plot_7 <- cbind(c(85.250, 85.625, 84.875, 81.625),   
                       c(93.250, 91.500, 91.625, 88.875),   
                       c(93.875, 92.275, 90.750, 87.250),
                       c(92.375, 90.000, 88.250, 81.375),
                       c(90.250, 88.000, 86.250, 76.500))


# Create data frames for each result set
result7 <- data.frame(
  q_Y = 1-c(0.01, 0.1, 0.2, 0.3, 0.4),
  `alpha_Y0.1` = result_plot_7[1,],
  `alpha_Y0.3` = result_plot_7[2,],
  `alpha_Y0.5` = result_plot_7[3,],
  `alpha_Y0.7` = result_plot_7[4,]
)

result_plot_8 <- cbind(c(91.875, 91.500, 91.625, 90.625),   
                       c(97.875, 97.500, 97.625, 96.875),   
                       c(99.750, 99.625, 99.250, 99.250),
                       c(99.500, 99.500, 97.750, 96.500),
                       c(99.250, 99.000, 97.250, 95.500))


result8 <- data.frame(
  q_Y = 1-c(0.01, 0.1, 0.2, 0.3, 0.4),
  `alpha_Y0.1` = result_plot_8[1,],
  `alpha_Y0.3` = result_plot_8[2,],
  `alpha_Y0.5` = result_plot_8[3,],
  `alpha_Y0.7` = result_plot_8[4,]
)

result_plot_9 <- rbind(c(98.75, 99.75, 100.00, 100.00, 99.75, 96.75),   
                       c(89.75, 97.00, 98.25, 98.50, 97.75, 93.75),   
                       c(79.75, 88.00, 96.25, 96.00, 95.25, 87.25),
                       c(59.25, 70.25, 76.75, 82.25, 87.50, 78.00))

result9 <- data.frame(
  q_z = 1- c(0, 0.01, 0.05, 0.1, 0.2, 0.5),
  `alpha_Z0.1` = result_plot_9[1,],
  `alpha_Z0.5` = result_plot_9[2,],
  `alpha_Z1` = result_plot_9[3,],
  `alpha_Z2` = result_plot_9[4,]
)















# FIRST PLOT (first 4 results)
data_long <- bind_rows(
  result1 %>% pivot_longer(cols = -alpha_X, names_to = "q_F", values_to = "value") %>% mutate(plot = "VAR heavy-tailed case"),
  result2 %>% pivot_longer(cols = -alpha_X, names_to = "q_F", values_to = "value") %>% mutate(plot = "VAR Gaussian case"),
  result3 %>% pivot_longer(cols = -alpha_X, names_to = "q_F", values_to = "value") %>% mutate(plot = "GARCH heavy-tailed case"),
  result4 %>% pivot_longer(cols = -alpha_X, names_to = "q_F", values_to = "value") %>% mutate(plot = "GARCH Gaussian case")
)
data_long$plot <- factor(data_long$plot, levels = c("VAR heavy-tailed case","VAR Gaussian case", "GARCH heavy-tailed case","GARCH Gaussian case"))
# Plot the data
ggplot(data_long, aes(x = alpha_X, y = value, color = q_F)) +
  geom_line(size = 0.7, linetype = "solid") +
  geom_point(size = 2, shape = 16) + # Use shape = 16 for solid circles
  ylab("Average error") +
  xlab(expression(alpha[X])) +
  facet_wrap(~plot, ncol = 2, scales = "free_x") + # Facet into a 2x2 grid
  scale_colour_manual(values = c('black', "#EE6677", "#228833", "#1f77b4"),
                      labels = c(expression(q[F] == 0), expression(q[F] == 0.3),
                                 expression(q[F] == 0.5), expression(q[F] == 0.7))) +
  theme(
    strip.text = element_text(size = 10), # Smaller font size for facet labels
    axis.text = element_text(size = 12), # Smaller font size for axis text
    axis.title = element_text(size = 15), # Smaller font size for axis titles
    legend.text = element_text(size = 10), # Smaller font size for legend text
    legend.title = element_blank(), # Hide legend title
    legend.position = "right", # Move legend to the right
    panel.grid.major = element_line(size = 0.5, linetype = "dashed", color = "grey80"),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(2, "lines"),
    aspect.ratio = 1 # Ensure each facet is squared
  ) 

#5.5x8.3 in export pdf















# SECOND PLOT (5 and 6 results)
data_long <- bind_rows(
  result5 %>% pivot_longer(cols = -k_n, names_to = "n", values_to = "value") %>% mutate(plot = "No hidden confounder"),
  result6 %>% pivot_longer(cols = -k_n, names_to = "n", values_to = "value") %>% mutate(plot = "Hidden confounder")
)
data_long$plot <- factor(data_long$plot, levels = c("No hidden confounder", "Hidden confounder"))

# Plot the data
ggplot(data_long, aes(x = k_n, y = value, color = n)) +
  geom_line(size = 0.7, linetype = "solid") +
  geom_point(size = 2, shape = 16) + # Use shape = 16 for solid circles
  ylab("Performance") +
  xlab(expression(k[n])) +
  facet_wrap(~plot, ncol = 2) + # Facet into a 1x2 grid
  scale_colour_manual(values = c('black', "#EE6677", "#228833"),
                      labels = c(expression(n == 200), expression(n == 400),
                                 expression(n == 600))) +
  theme(
    strip.text = element_text(size = 10), # Smaller font size for facet labels
    axis.text = element_text(size = 12), # Smaller font size for axis text
    axis.title = element_text(size = 12), # Smaller font size for axis titles
    legend.text = element_text(size = 10), # Smaller font size for legend text
    legend.title = element_blank(), # Hide legend title
    legend.position = "right", # Move legend to the right
    panel.grid.major = element_line(size = 0.5, linetype = "dashed", color = "grey80"),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(2, "lines"),
    aspect.ratio = 1 # Ensure each facet is squared
  ) +
  scale_x_continuous(breaks = seq(0.2, 0.6, by = 0.1)) + # Regular breaks for x-axis
  scale_y_continuous(breaks = seq(45, 95, by = 10)) # Regular breaks for y-axis

#3x5 export pdf











# THIRD PLOT (7 and 8 results)
data_long <- bind_rows(
  result7 %>% pivot_longer(cols = -q_Y, names_to = "alpha_Y", values_to = "value") %>% mutate(plot = "n=500"),
  result8 %>% pivot_longer(cols = -q_Y, names_to = "alpha_Y", values_to = "value") %>% mutate(plot = "n=10000")
)
data_long$plot <- factor(data_long$plot, levels = c("n=500", "n=10000"))

# Plot the data
ggplot(data_long, aes(x = q_Y, y = value, color = alpha_Y)) +
  geom_line(size = 0.7, linetype = "solid") +
  geom_point(size = 2, shape = 16) + # Use shape = 16 for solid circles
  ylab("Performance") +
  xlab(expression(q[Y])) +
  facet_wrap(~plot, ncol = 2) + # Facet into a 1x2 grid
  scale_colour_manual(values = c('black', "#EE6677", "#228833", "#1f77b4"),
                      labels = c(expression(alpha[Y] == 0.1), expression(alpha[Y] == 0.3),
                                 expression(alpha[Y] == 0.5), expression(alpha[Y] == 0.7))) +
  theme(
    strip.text = element_text(size = 10), # Smaller font size for facet labels
    axis.text = element_text(size = 12), # Smaller font size for axis text
    axis.title = element_text(size = 12), # Smaller font size for axis titles
    legend.text = element_text(size = 10), # Smaller font size for legend text
    legend.title = element_blank(), # Hide legend title
    legend.position = "right", # Move legend to the right
    panel.grid.major = element_line(size = 0.5, linetype = "dashed", color = "grey80"),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(1, "lines"),
    aspect.ratio = 1 # Ensure each facet is squared
  )

#3x5 export pdf















# LAST PLOT (9th result)
data_long <- result9 %>% 
  pivot_longer(cols = -q_z, names_to = "alpha_Z", values_to = "value")

ggplot(data_long, aes(x = q_z, y = value, color = alpha_Z)) +
  geom_line(size = 0.7, linetype = "solid") +
  geom_point(size = 2, shape = 16) + # Use shape = 16 for solid circles
  ylab("Performance") +
  xlab(expression(q[Z])) +
  ggtitle(expression("Choice of " * q[Z])) + # Add title with math expression
  scale_colour_manual(values = c('black', "#EE6677", "#228833", "#1f77b4"),
                      labels = c(expression(alpha[Z] == 0.1), expression(alpha[Z] == 0.5),
                                 expression(alpha[Z] == 1), expression(alpha[Z] == 2))) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5), # Center title and adjust size/style
    strip.text = element_text(size = 10), # Smaller font size for facet labels
    axis.text = element_text(size = 12), # Smaller font size for axis text
    axis.title = element_text(size = 12), # Smaller font size for axis titles
    legend.text = element_text(size = 10), # Smaller font size for legend text
    legend.title = element_blank(), # Hide legend title
    legend.position = "right", # Move legend to the right
    panel.grid.major = element_line(size = 0.5, linetype = "dashed", color = "grey80"),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1 # Ensure each facet is squared
  )
#3x5 export pdf


