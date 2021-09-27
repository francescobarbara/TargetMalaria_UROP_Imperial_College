# This one was an earlier approach for a dissimilarity function, it uses squared distances
# and returns the mean squared error


compare = function(estimates, real_values, recombination_rates, simulations, p0_estimate, 
                   q0_estimate, generations_to_sim = 100){
  
  #INPUTS:
  # estimates is a 3x1 vector, represents the (N,M,migration) estimate you wanna 
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
  # A tuple of the form c(y, x) where x is the value of the dissimilarity 
  # measure developed
  
  N = estimates[1]
  M = estimates[2]
  migr = estimates[3]
  migration = matrix(c(1-migr, migr, migr, 1-migr), nrow = 2)
  
  n = length(recombination_rates)
  residuals = matrix(0, n, 3)
  for (i in 1:n){
    c = recombination_rates[i]
    # initializing the matrix that will store the simulated 3 LD measures
    simulated_r_sq = matrix(0, simulations, 3)
    for (j in 1:simulations){
      
      out = simulator_2_pop(p0_estimate, q0_estimate, rep(N, generations_to_sim), 
                            rep(M, generations_to_sim), c, migration)
      p_values = matrix(unlist(out[1]), nrow = 4)[1:4,generations_to_sim]
      q_values = matrix(unlist(out[2]), nrow = 4)[1:4,generations_to_sim]
      r_sq_pop1 = matrix(unlist(out[5]), nrow = 4)[generations_to_sim]
      r_sq_pop2 = matrix(unlist(out[6]), nrow = 4)[generations_to_sim]
      l_sq = likelihood_test(p_values, q_values, N, M, 2)
      # adding row in the matrix
      simulated_r_sq[j,] = c(l_sq, r_sq_pop1, r_sq_pop2)
    }
    
    #calculating means and std devs for the 3 measures
    mu1 = mean(simulated_r_sq[,1])
    mu2 = mean(simulated_r_sq[,2])
    mu3 = mean(simulated_r_sq[,3])
    
    sigma1 = sd(simulated_r_sq[,1]) 
    sigma2 = sd(simulated_r_sq[,2])
    sigma3 = sd(simulated_r_sq[,3])
    
    res1 = ((real_values[i][1] - mu1)/ sigma1)^2
    res2 = ((real_values[i][1] - mu2)/ sigma2)^2
    res3 = ((real_values[i][1] - mu3)/ sigma3)^2
    
    residuals[i,] = c(res1, res2, res3)
  }
  #returning the mse and the test value (you have a chisq distribution with 3n dof)
  sum_of_res = sum(residuals)
  test = 1 - pchisq(sum_of_res, 3*n)
  mse = sum_of_res / (3*n)
  
  return( c(test, mse))
}

