#This code serves for reproducibility purposes concerning the simulations about the choice of hyper-parameters
#You need to first upload the function 'Extreme_causality_test' from the other file
#Lines 1-200 corresponds to the choice of q_F; that is, optimal choice of F
#Lines 200-400 corresponds to the choice of tau_X
#Lines 400-600 corresponds to the choice of tau_Y
#Lines 600-800 corresponds to the choice of tau_Z
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





