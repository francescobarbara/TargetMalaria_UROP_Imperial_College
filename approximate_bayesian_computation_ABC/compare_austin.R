compare = function(estimates, real_values, recombination_rates, simulations, p0_estimate, 
                   q0_estimate, generations_to_sim = 100){
  
  #INPUTS:
  # estimates is a 3x1 vector, represents the (N,M,migratio) estimate you wanna 
  #                             compare to the real data
  
  # real_values is a length(recombination_rates)x3 matrix (usually generated with
  # " create_synthetic_observed_data "); represents the oserved 3 measures for every
  # pair of loci
  
  # recombination_rates are the known rec. rates.
  # REMARK: Here and in the final analyses I only used recombination rates 
  # 0.05, 0.10, ..., 0.45, 0.50 becuase this way I only have to simulate 10 
  # approximated cdfs (speed otherwise would have made it impossible)
  # However, if you look for example at 'old_compare', there I have a better but
  # extremely more expensive approach that works with all the recombination rates
  
  # simulations: are the number of points you want to use when creating each approximate
  # cdf, in the presentation I showed that 1000 is very good, 100 is okay-ish, but 
  # the variance in 3-4 times bigger than expected in this case
  
  # p0_estimate, just use c(0.25, ..., 0.25) it's the initial 4 freq to use in the cpp code
  # same for population 2 is q0_estimate
  
  # generations_to_sim in the cpp code, I generally use 50, to make sure c = 0.05 
  # has converged
  
  # OUTPUT:
  # A tuple of the form c(0, x) where x is the value of the dissimilarity 
  # measure developed
  
  N = estimates[1]
  M = estimates[2]
  migr = estimates[3]
  migration = matrix(c(1-migr, migr, migr, 1-migr), nrow = 2)
  
  # this vectors will store the simulated values for the 3 measures
  ecdfs1 = matrix(0,simulations,10)
  ecdfs2 = matrix(0,simulations,10)
  ecdfs3 = matrix(0,simulations,10)
  
  for(i in 1:10){
    # I am only using 'multiples of 5' recombination rates
    # so I only need to store 10 ecdfs for each of the 3 measures
    c = 0.05*i
    simulated_r_sq = matrix(0, simulations, 3)
    for (j in 1:simulations){
      
      # Here c is fixed, and am creating the 100/1000 datapoints that will make
      # up the cdf
      
      out = simulator_2_pop(p0_estimate, q0_estimate, rep(N, generations_to_sim), 
                            rep(M, generations_to_sim), c, migration)
      p_values = matrix(unlist(out[1]), nrow = 4)[1:4,generations_to_sim]
      q_values = matrix(unlist(out[2]), nrow = 4)[1:4,generations_to_sim]
      # calculating the three measures, r^2 is calculated in C++, so just need
      # to take it out of the output 'likelihood_test'
      # For l^2 I need to call 
      r_sq_pop1 = unlist(out[5])[generations_to_sim]
      r_sq_pop2 = unlist(out[6])[generations_to_sim]
      l_sq = likelihood_test(p_values, q_values, N, M, 2)
      
      #storing these measures
      simulated_r_sq[j,] = c(l_sq, r_sq_pop1, r_sq_pop2)
    }
    
    # small subtlety here, in the r^2 you might be dividing by zero, turning these
    # NaN values into 1 (essentially this happens if you have like 50 individuals
    # and hence p1 could be equal to zero)
    
    if (any(is.na(simulated_r_sq[,2]))){
      simulated_r_sq[is.na(simulated_r_sq[,2]) , 2] = 1
    }
    if (any(is.na(simulated_r_sq[,3]))){
      simulated_r_sq[is.na(simulated_r_sq[,3]) , 3] = 1
    }
    
    # appending the datapoints for a specific recombination rat
    ecdfs1[,i] = simulated_r_sq[,1]
    ecdfs2[,i] = simulated_r_sq[,2]
    ecdfs3[,i] = simulated_r_sq[,3]
  }
  
  
  
  # Now I calculate the dissimilarity measure using the approximate cdfs I am
  # storing in ecdfs1, 2, 3
  
  n = length(recombination_rates)
  residuals = matrix(0, n, 3)
  for (i in 1:n){
    #tell me which column in the ecdf matrix I need to use
    # 0.05 -> 1, 0.10 -> 2 etc
    
    index = recombination_rates[i]/0.05
    
    # We have the vector of datapoints, we now use built-in function ecdf
    # to create an approximate cdf
    ecdf1 = ecdf(ecdfs1[,index])
    ecdf2 = ecdf(ecdfs2[,index])
    ecdf3 = ecdf(ecdfs3[,index])
    
    #calculating the distance from the median in a weighted manner
    res1 = 8*(ecdf1(real_values[i,1]) - 0.5) 
    res2 = 1*(ecdf2(real_values[i,2]) - 0.5)
    res3 = 1*(ecdf3(real_values[i,3]) - 0.5) 
    
    residuals[i,] = c(res1, res2, res3)
  }
  
  
  sum_of_res = abs(sum(residuals[,1])) + abs(sum(residuals[,2])) + abs(sum(residuals[,3]))
  test = 0
  mse = sum_of_res / (10*n)
  
  return( c(test, mse))
}
