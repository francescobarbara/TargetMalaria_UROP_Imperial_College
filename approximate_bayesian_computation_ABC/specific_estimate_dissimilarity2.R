#ANOTHER COMPARISON, just simulationg synthetic data,
# then compare with a chosen estimate

N = 1000
M = 1000
migration = matrix(c(0.9, 0.1, 0.1,0.9), nrow = 2)
recombination_rates = c(0.1, 0.05, 0.4)
generations_to_sim = 100

p0 = c(0.25,0.25,0.25,0.25)
q0 = c(0.25,0.25,0.25,0.25)

real_values = create_synthetic_observed_data(N, M, migration, recombination_rates,
                                                        generations_to_sim, p0, q0)


estimates = c(1000, 1000, 0.4)  #try 10,000 and m = 0.5
simulations = 20
p0_estimate = c(0.25,0.25,0.25,0.25)
q0_estimate = c(0.25,0.25,0.25,0.25)

compare(estimates, real_values, recombination_rates, simulations, p0_estimate, 
                   q0_estimate, generations_to_sim = 100)
