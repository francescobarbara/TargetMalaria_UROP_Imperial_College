# Here I am trying to see how much the uniform approximation holds, and the consequent
# Normal test that I introduced to estimate the probability that 'real data can come
# from an estimated triplet'
# Essentially the question is, how many points do I need to simulate to get an
# empirical cdf that behaves well?
# To this aim, for each i in [1,...,no_sim] , I store the value of the dissimilarity easure
# After all the simulatons, I check to see whether the variance is the same
# as the assumed one under the Uniform assumpution

simulations = 1000
pairs_loci = 100
###############################################################

N = 1000
M = 1000
migration = matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)
#recombination_rates = c(0.1, 0.05, 0.4, 0.15, 0.5, 0.35, 0.2, 0.11, 0.175, 0.45)
recombination_rates = sample(c(0.05, 0.1, 0.15,0.2, 0.25,0.3, 0.35, 0.4, 0.45, 0.5),
                             size = pairs_loci, replace = TRUE)

generations_to_sim = 50

p0 = c(0.25,0.25,0.25,0.25)
q0 = c(0.25,0.25,0.25,0.25)
p0_estimate = c(0.25,0.25,0.25,0.25)
q0_estimate = c(0.25,0.25,0.25,0.25)
###########################################################################

real_values = create_synthetic_observed_data(N, M, migration, recombination_rates,
                                             generations_to_sim, p0, q0)

estimates = c(1000, 1000, 0.1)
compare(estimates, real_values, recombination_rates, simulations, p0_estimate, 
        q0_estimate, generations_to_sim)

2*(1 - pnorm(sqrt(12*100 / 0.66) * 0.07))

no_sim = 100
resulted_diss = numeric(no_sim)
for (i in 1:no_sim){
  print(i)
  real_values = create_synthetic_observed_data(N, M, migration, recombination_rates,
                                               generations_to_sim, p0, q0)
  resulted_diss[i] = compare(estimates, real_values, recombination_rates, simulations, p0_estimate, 
                             q0_estimate, generations_to_sim)[2]
}

mean(resulted_diss) #0.0014  -0.005    -0.007
sd(resulted_diss)  #0.081   0,012      0.0085
sqrt(0.66/(12*1000)) #0.02   0.0074    0.0074 expected

hist(resulted_diss, breaks = 30, col = 'lightblue', xlab = 'simulated dissimilarity values', 
     main = 'Dissimilarity fnc, 1000 simulations for 1000 pairs ol loci')
