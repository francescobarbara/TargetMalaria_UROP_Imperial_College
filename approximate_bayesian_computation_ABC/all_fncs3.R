#FUNCTIONS THAT YOU NEED

library("Rcpp")
Rcpp::sourceCpp("C:/Users/angus/OneDrive/Desktop/urop/C++/simulation_3.cpp")

############################################################################

quadratic <- function(a,b,c){
  if(delta(a,b,c) > 0){ # first case D>0
    x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
    x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
    result = c(x_1,x_2)
    return(result)
  }
  else if(delta(a,b,c) == 0){ # second case D=0
    x = -b/(2*a)
    return(x)
  }
  else {"There are no real roots."} # third case D<0
}

# Constructing delta
delta<-function(a,b,c){
  b^2-4*a*c
}

###########################################################################

## NEEDED TO CREATE THE L^2 

likelihood_test = function(p, q, size1, size2, dof = 2){
  # Given the 2 vector of alleles frequencies and the final size2 
  # of the two population, and the degrees of freedom of the chi_sq distribution
  # It returns the p-value for the likelihood ratio test I developed (see pdfs),
  # which is indeed our joint r^2 measure
  
  # creating the frequencies of the 2x2x2 cube from the single-population frequencies
  p = p*size1/(size1+size2)
  q = q*size2/(size1+size2)
  #print(sum(p) + sum(q))
  
  #creating a 3x2 matrix of the marginals for the 2x2x2 contingency table
  marginals = marginals(p,q)
  #print(marginals)
  
  #initializing the vectors that contain the predicted frequencies under H0
  # i.e. they are the product of the 3 marginals for each cell
  exp_p = numeric(4)
  exp_q = numeric(4)
  
  exp_p[1] = marginals[1,1]*marginals[2,1]*marginals[3,1]
  exp_p[2] = marginals[1,2]*marginals[2,1]*marginals[3,1]
  exp_p[3] = marginals[1,1]*marginals[2,2]*marginals[3,1]
  exp_p[4] = marginals[1,2]*marginals[2,2]*marginals[3,1]
  
  exp_q[1] = marginals[1,1]*marginals[2,1]*marginals[3,2]
  exp_q[2] = marginals[1,2]*marginals[2,1]*marginals[3,2]
  exp_q[3] = marginals[1,1]*marginals[2,2]*marginals[3,2]
  exp_q[4] = marginals[1,2]*marginals[2,2]*marginals[3,2]
  
  #print(c('marginals', marginals))
  #print(c('exp_p', exp_p))
  #print(c('exp_q',exp_q))
  
  #calculating the likelihood ratio test statistic (as illustrated in the pdf)
  test_statistic = 0
  for (i in 1:4){
    if (p[i] != 0){
      test_statistic = test_statistic + 2*p[i]*log(p[i]/exp_p[i])
    }
    if (q[i] != 0){
      test_statistic = test_statistic + 2*q[i]*log(q[i]/exp_q[i])
    }
  }
  #print(test_statistic)
  return(pchisq(test_statistic, df = dof))
}




marginals = function(p,q){
  # Given two vectors of size 4 (the allele frequencies for the 2 sub-populations,
  # which we already normalize by size in the fnc above),
  # it creates a 2x2x2 contingency table and return the marginals as a 3x2 matrix
  arr = array(c(p, q), dim = c(2,2,2))
  
  marginals = matrix(nrow = 3,ncol = 2)
  
  marginals[1,1] = arr[1,1,1] + arr[1,2,1] + arr[1,1,2] + arr[1,2,2]
  marginals[1,2] = arr[2,1,1] + arr[2,2,1] + arr[2,1,2] + arr[2,2,2]
  marginals[2,1] = arr[1,1,1] + arr[2,1,1] + arr[1,1,2] + arr[2,1,2]
  marginals[2,2] = arr[1,2,1] + arr[2,2,1] + arr[1,2,2] + arr[2,2,2]
  marginals[3,1] = arr[1,1,1] + arr[2,1,1] + arr[1,2,1] + arr[2,2,1]
  marginals[3,2] = arr[1,1,2] + arr[2,1,2] + arr[1,2,2] + arr[2,2,2]
  
  return(marginals)
}


################################################################################

simulator_2_pop = function(p0, q0, N, M, c, migration){
  
  
  return (simulation(p0, q0, N, M, c, migration))
  
}  


################################################################################

create_synthetic_observed_data = function(N, M, migration, recombination_rates,
                                          generations_to_sim, p0, q0){
  n = length(recombination_rates)
  real_values = matrix(0, n, 3)
  
  for (j in 1:n){
    c = recombination_rates[j]
    out = simulator_2_pop(p0, q0, rep(N, generations_to_sim), 
                          rep(M, generations_to_sim), c, migration)
    p_values = matrix(unlist(out[1]), nrow = 4)[1:4,generations_to_sim]
    q_values = matrix(unlist(out[2]), nrow = 4)[1:4,generations_to_sim]
    r_sq_pop1 = unlist(out[5])[generations_to_sim]
    r_sq_pop2 = unlist(out[6])[generations_to_sim]
    l_sq = likelihood_test(p_values, q_values, N, M, 2)
    real_values[j,] = c(l_sq, r_sq_pop1, r_sq_pop2)
  }
  return(real_values)
}

####################################################################################################################

######################################################################################

simulate_exp_less_05 = function(number_samples, rate = 5){
  out = numeric(number_samples)
  for (i in 1:number_samples){
    sim_val = 1
    while (sim_val > 0.5){
      sim_val = rexp(1, rate = rate)
    }
    out[i] = sim_val
  }
  return( out )
} 

######################################################################################
######################################################################################
######################################################################################
## NEW APPROACH WITH ECDF

compare = function(estimates, real_values, recombination_rates, simulations, p0_estimate, 
                   q0_estimate, generations_to_sim = 100){
  
  N = estimates[1]
  M = estimates[2]
  migr = estimates[3]
  migration = matrix(c(1-migr, migr, migr, 1-migr), nrow = 2)
  
  ecdfs1 = matrix(0,simulations,10)
  ecdfs2 = matrix(0,simulations,10)
  ecdfs3 = matrix(0,simulations,10)
  
  for(i in 1:10){
    c = 0.05*i
    simulated_r_sq = matrix(0, simulations, 3)
    for (j in 1:simulations){
      
      out = simulator_2_pop(p0_estimate, q0_estimate, rep(N, generations_to_sim), 
                            rep(M, generations_to_sim), c, migration)
      p_values = matrix(unlist(out[1]), nrow = 4)[1:4,generations_to_sim]
      q_values = matrix(unlist(out[2]), nrow = 4)[1:4,generations_to_sim]
      r_sq_pop1 = unlist(out[5])[generations_to_sim]
      r_sq_pop2 = unlist(out[6])[generations_to_sim]
      l_sq = likelihood_test(p_values, q_values, N, M, 2)
      simulated_r_sq[j,] = c(l_sq, r_sq_pop1, r_sq_pop2)
    }
    
    if (any(is.na(simulated_r_sq[,2]))){
      simulated_r_sq[is.na(simulated_r_sq[,2]) , 2] = 1
    }
    if (any(is.na(simulated_r_sq[,3]))){
      simulated_r_sq[is.na(simulated_r_sq[,3]) , 3] = 1
    }
    ecdfs1[,i] = simulated_r_sq[,1]
    ecdfs2[,i] = simulated_r_sq[,2]
    ecdfs3[,i] = simulated_r_sq[,3]
  }
  
  
  
  
  n = length(recombination_rates)
  residuals = matrix(0, n, 3)
  for (i in 1:n){
    index = recombination_rates[i]/0.05
    
    ecdf1 = ecdf(ecdfs1[,index])
    ecdf2 = ecdf(ecdfs2[,index])
    ecdf3 = ecdf(ecdfs3[,index])
    
    
    res1 = 8*(ecdf1(real_values[i,1]) - 0.5) 
    res2 = 1*(ecdf2(real_values[i,2]) - 0.5)
    res3 = 1*(ecdf3(real_values[i,3]) - 0.5) 
    
    residuals[i,] = c(res1, res2, res3)
  }
  sum_of_res = abs(sum(residuals))
  test = 0
  mse = sum_of_res / (10*n)
  
  return( c(test, mse))
}


####################################################################################################################
df_and_plots = function(full_mat, values_per_layer, simulated_values_per_layer, test_number = 666, save_df = FALSE){
  par(pty="s")
  df = data.frame(full_mat)
  colnames(df) = c('population_1', 'population_2', 'migration', 'similarity', 'mse')
  
  rbPal <- colorRampPalette(c('red','blue'))
  
  #This adds a column of color values
  # based on the y values
  df$col <- rbPal(1000)[as.numeric(cut(df$mse,breaks = 1000))]
  
  plot(df$population_1, df$population_2, pch = 16, col = df$col, log='xy', 
       xlab = 'Population 1 estimate', ylab = 'Population 2 estimate', main = 'Dissimilarity of simulated triplets')
  
  plot(df$population_1, df$migration, pch = 16, col = df$col, log='x', 
       xlab = 'Population 1 estimate', ylab = 'Migration estimate', main = 'Dissimilarity of simulated triplets')
  
  plot(df$population_2, df$migration, pch = 16, col = df$col, log='x', 
       xlab = 'Population 2 estimate', ylab = 'Migration estimate', main = 'Dissimilarity of simulated triplets')
  #pch = 20 è più piccolo
  
  # Now plotting the layers as colours
  
  n = dim(full_mat)[1]
  layer_index = numeric(n)
  
  for (i in 1:values_per_layer){
    layer_index[i] = 0
  }
  number_layers = (n - values_per_layer) / simulated_values_per_layer
  
  for (j in 1:number_layers){
    for( i in (4 + simulated_values_per_layer*(j-1) + 1) : (4 + simulated_values_per_layer*j)){
      layer_index[i] = j
    }
  }
  
  rbPal <- colorRampPalette(c('blue', 'red'))
  
  df$layer = layer_index
  
  df$col2 = rbPal(number_layers)[as.numeric(cut((df$layer),breaks = number_layers))]
  
  plot(df$population_1, df$population_2, pch = 16, col = df$col2, log='xy', 
       xlab = 'Population 1 estimate', ylab = 'Population 2 estimate', main = 'Layers convergence')
  
  plot(df$population_1, df$migration, pch = 16, col = df$col2, log='x', 
       xlab = 'Population 1 estimate', ylab = 'Migration estimate', main = 'Layers convergence')
  
  plot(df$population_2, df$migration, pch = 16, col = df$col2, log='x', 
       xlab = 'Population 2 estimate', ylab = 'Migration estimate', main = 'Layers convergence')
  df = df[,c(1,2,3,4,5,7)]
  
  if (save_df){
    
    
    path = paste("C:/Users/angus/OneDrive/Desktop/urop/ABC/tests_csv/test_", as.character(test_number), '.csv')
    write.csv(x=df, file= path)
  }
  return (df)
}

###################################################################################################################
####################################################################
#DA USARE PER ZERO
sequential_mc_exp_migration = function(real_values, recombination_rates, values_per_layer, 
                                       simulated_values_per_layer, simulations_per_est, max_depth, improvement_threshold,
                                       p0_estimate, q0_estimate, generations_to_sim = 100){
  
  n = length(recombination_rates)
  
  upperbound_pop1 = numeric(n)
  upperbound_pop2 = numeric(n)
  
  for (i in 1:n){
    a = 2*(real_values[i,2]) - 2*(real_values[i,2]) * ((1-recombination_rates[i])^2)
    b = (real_values[i,2]) * ((1-recombination_rates[i])^2) * 
      (2*recombination_rates[i] + 1) - (1+recombination_rates[i]^2)
    d = real_values[i,2] * (1-recombination_rates[i])^2 * recombination_rates[i]
    
    upperbound_pop1[i] = quadratic(a,b,d)[1]
    
    a = 2*(real_values[i,3]) - 2*(real_values[i,3]) * ((1-recombination_rates[i])^2)
    b = (real_values[i,3]) * ((1-recombination_rates[i])^2) * 
      (2*recombination_rates[i] + 1) - (1+recombination_rates[i]^2)
    d = real_values[i,3] * (1-recombination_rates[i])^2 * recombination_rates[i]
    
    upperbound_pop2[i] = quadratic(a,b,d)[1]
  }
  
  #upperbound_pop1 = 2*max(upperbound_pop1)
  #upperbound_pop2 = 2*max(upperbound_pop2)
  
  upperbound_pop1 = 1*median(upperbound_pop1)  #to change in 3
  upperbound_pop2 = 1*median(upperbound_pop2)
  
  print(c(upperbound_pop1, upperbound_pop2))
  estimates_previous_layer = matrix(0, values_per_layer, 5)
  #N, M, migr, p-value, mse
  estimates_previous_layer[,1] = runif(values_per_layer, 0, upperbound_pop1)
  estimates_previous_layer[,2] = runif(values_per_layer, 0, upperbound_pop2)
  estimates_previous_layer[,3] = simulate_exp_less_05(values_per_layer, rate = 5)
  for (j in 1:values_per_layer){
    print(c(j, 'out of', values_per_layer, 'layer 1')) 
    estimates_previous_layer[j,c(4,5)] = compare(estimates_previous_layer[j,c(1,2,3)],
                                                 real_values, recombination_rates, simulations_per_est, p0_estimate, 
                                                 q0_estimate, generations_to_sim)
    
  }
  
  all_estimates = estimates_previous_layer
  tree_depth_count = 1
  improvement = Inf
  best_similarity = max(estimates_previous_layer[,4])
  lowest_mse = min(estimates_previous_layer[,5])
  standard_deviations = c(sd(estimates_previous_layer[,1]), sd(estimates_previous_layer[,2]), sd(estimates_previous_layer[,3]))
  
  while((tree_depth_count <= max_depth) & (improvement > improvement_threshold)){
    
    tree_depth_count = tree_depth_count +1
    print(c('implementing layer number', tree_depth_count))
    layer_best_similarity = 0
    layer_lowest_mse = Inf
    new_estimates = matrix(0, simulated_values_per_layer, 5)
    
    for (i in 1:simulated_values_per_layer){
      print(c(i, 'out of', values_per_layer, 'layer', tree_depth_count))
      index = sample(1:values_per_layer, 1, prob = 1/estimates_previous_layer[,5])
      vals = estimates_previous_layer[index,c(1,2,3)] + 1*standard_deviations * rnorm(3)    #### adjust?
      
      while ( !( (vals[1] > 0) & (vals[2] > 0)  & (vals[3] <= 0.5)  ) ){   #& (vals[3] > 0)
        
        vals = estimates_previous_layer[index,c(1,2,3)] + c(1,1,1)*standard_deviations * rnorm(3)  
      }
      
      if (vals[3] < 0){
        vals[3] = 0
      }
      
      new_estimates[i, c(1,2,3)] = vals
      
      new_estimates[i, c(4,5)] = compare(vals, real_values, recombination_rates, simulations_per_est, p0_estimate, 
                                         q0_estimate, generations_to_sim)
    }
    
    layer_best_similarity = max(layer_best_similarity, new_estimates[i, 4])
    layer_lowest_mse = min(layer_lowest_mse, new_estimates[i,5])
    improvement = (lowest_mse - layer_lowest_mse)/lowest_mse
    
    
    new_estimates = new_estimates[order(new_estimates[,5], decreasing = FALSE),]
    all_estimates = rbind(all_estimates, new_estimates)
    
    estimates_previous_layer = new_estimates[1:values_per_layer,]
    standard_deviations = c(sd(estimates_previous_layer[,1]), sd(estimates_previous_layer[,2]), sd(estimates_previous_layer[,3]))
    
  }
  
  #all_estimates = all_estimates[order(all_estimates[,5], decreasing = FALSE),]
  
  return(all_estimates)
}

######################################################################################


simulate_exp_less_05 = function(number_samples, rate = 5){
  out = numeric(number_samples)
  for (i in 1:number_samples){
    sim_val = 1
    while (sim_val > 0.5){
      sim_val = rexp(1, rate = rate)
    }
    out[i] = sim_val
  }
  return( out )
} 

######################################################################################

sequential_mc_exp_migration = function(real_values, recombination_rates, values_per_layer, 
                                       simulated_values_per_layer, simulations_per_est, max_depth, improvement_threshold,
                                       p0_estimate, q0_estimate, generations_to_sim = 100){
  
  n = length(recombination_rates)
  
  upperbound_pop1 = numeric(n)
  upperbound_pop2 = numeric(n)
  
  for (i in 1:n){
    a = 2*(real_values[i,2]) - 2*(real_values[i,2]) * ((1-recombination_rates[i])^2)
    b = (real_values[i,2]) * ((1-recombination_rates[i])^2) * 
      (2*recombination_rates[i] + 1) - (1+recombination_rates[i]^2)
    d = real_values[i,2] * (1-recombination_rates[i])^2 * recombination_rates[i]
    
    upperbound_pop1[i] = quadratic(a,b,d)[1]
    
    a = 2*(real_values[i,3]) - 2*(real_values[i,3]) * ((1-recombination_rates[i])^2)
    b = (real_values[i,3]) * ((1-recombination_rates[i])^2) * 
      (2*recombination_rates[i] + 1) - (1+recombination_rates[i]^2)
    d = real_values[i,3] * (1-recombination_rates[i])^2 * recombination_rates[i]
    
    upperbound_pop2[i] = quadratic(a,b,d)[1]
  }
  
  #upperbound_pop1 = 2*max(upperbound_pop1)
  #upperbound_pop2 = 2*max(upperbound_pop2)
  
  upperbound_pop1 = 1*median(upperbound_pop1)  #1
  upperbound_pop2 = 1*median(upperbound_pop2)  #1
  
  print(c(upperbound_pop1, upperbound_pop2))
  estimates_previous_layer = matrix(0, values_per_layer, 5)
  #N, M, migr, p-value, mse
  estimates_previous_layer[,1] = runif(values_per_layer, 0, upperbound_pop1)
  estimates_previous_layer[,2] = runif(values_per_layer, 0, upperbound_pop2)
  estimates_previous_layer[,3] = simulate_exp_less_05(values_per_layer, rate = 5)
  for (j in 1:values_per_layer){
    print(c(j, 'out of', values_per_layer, 'layer 1')) 
    estimates_previous_layer[j,c(4,5)] = compare(estimates_previous_layer[j,c(1,2,3)],
                                                 real_values, recombination_rates, simulations_per_est, p0_estimate, 
                                                 q0_estimate, generations_to_sim)
    
  }
  
  all_estimates = estimates_previous_layer
  tree_depth_count = 1
  improvement = Inf
  best_similarity = max(estimates_previous_layer[,4])
  lowest_mse = min(estimates_previous_layer[,5])
  standard_deviations = c(sd(estimates_previous_layer[,1]), sd(estimates_previous_layer[,2]), sd(estimates_previous_layer[,3]))
  
  while((tree_depth_count <= max_depth) & (improvement > improvement_threshold)){
    
    tree_depth_count = tree_depth_count +1
    print(c('implementing layer number', tree_depth_count))
    layer_best_similarity = 0
    layer_lowest_mse = Inf
    new_estimates = matrix(0, simulated_values_per_layer, 5)
    
    for (i in 1:simulated_values_per_layer){
      print(c(i, 'out of', values_per_layer, 'layer', tree_depth_count))
      index = sample(1:values_per_layer, 1, prob = 1/estimates_previous_layer[,5])
      vals = estimates_previous_layer[index,c(1,2,3)] + 1*standard_deviations * rnorm(3)
      
      while ( !( (vals[1] > 0) & (vals[2] > 0) & (vals[3] > 0) & (vals[3] <= 0.5)  ) ){
        
        vals = estimates_previous_layer[index,c(1,2,3)] + c(1,1,1)*standard_deviations * rnorm(3)  
      }
      
      
      new_estimates[i, c(1,2,3)] = vals
      
      new_estimates[i, c(4,5)] = compare(vals, real_values, recombination_rates, simulations_per_est, p0_estimate, 
                                         q0_estimate, generations_to_sim)
      print(new_estimates[i,])
    }
    
    layer_best_similarity = max(layer_best_similarity, new_estimates[i, 4])
    layer_lowest_mse = min(layer_lowest_mse, new_estimates[i,5])
    improvement = (lowest_mse - layer_lowest_mse)/lowest_mse
    
    
    new_estimates = new_estimates[order(new_estimates[,5], decreasing = FALSE),]
    all_estimates = rbind(all_estimates, new_estimates)
    
    estimates_previous_layer = new_estimates[1:values_per_layer,]
    standard_deviations = c(sd(estimates_previous_layer[,1]), sd(estimates_previous_layer[,2]), sd(estimates_previous_layer[,3]))
    
  }
  
  #all_estimates = all_estimates[order(all_estimates[,5], decreasing = FALSE),]
  
  return(all_estimates)
}

######################################################################################
