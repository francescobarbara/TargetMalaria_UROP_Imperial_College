
# Here I wanna see, what the algorithm outputs for the dissimilarity of
# a specifically given estimate.
# I pick the true parameters, create the synthetic real data, then choose my estimate
# and see what the value of the dissimilarity function is

##########################
N = 1000
M = 1000
migration = matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)
recombination_rates = c(0.1, 0.05, 0.4, 0.15, 0.5, 0.35, 0.2, 0.11, 0.175, 0.45)
generations_to_sim = 50
p0 = c(0.25,0.25,0.25,0.25)
q0 = c(0.25,0.25,0.25,0.25)
p0_estimate = c(0.25,0.25,0.25,0.25)
q0_estimate = c(0.25,0.25,0.25,0.25)
simulations = 100


recombination_rates = sample(c(0.05, 0.1, 0.15,0.2, 0.25,0.3, 0.35, 0.4, 0.45, 0.5),
                             size = 1000, replace = TRUE)
real_values = create_synthetic_observed_data(N, M, migration, recombination_rates,
                                             generations_to_sim, p0, q0)
estimates = c(1000, 1000, 0.1)

MSE = compare(estimates, real_values, recombination_rates, simulations, p0_estimate, 
                   q0_estimate, generations_to_sim)


