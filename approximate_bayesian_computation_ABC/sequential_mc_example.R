
# TEST 


# initializing the true parameters, we want to predict
N = 1000
M = 1000
migration = matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)

#do not use this random rec rates if you are using the 'all_fncs4.R' compare function
#recombination_rates = c(0.1, 0.05, 0.4, 0.15, 0.5, 0.35, 0.2, 0.11, 0.175, 0.45)

# use this instead
recombination_rates = sample(c(0.05, 0.1, 0.15,0.2, 0.25,0.3, 0.35, 0.4, 0.45, 0.5),
                             size = 1000, replace = TRUE)


# values_per_layer is the number of 'good points' you keep at each layer
# simulated_values_per_layer is the number of points you simulate at each layer

# simulations_per_est is the number of points you want to use to create
# the approximate cdf

# max_depth is the maximum number of layers you wanna simulate

# improvement threshold, tells you to stop if you are not improvin by at leeast
# a factor of the in-puted value, I always put -100, so that the algo doesn't stop early

# p0_estimate, q0_estimate the initial frequencies you want your C++ code to start simulating from
# generations_to_sim = 100 the number of genrations you want your C++ codeto simulate

generations_to_sim = 50

p0 = c(0.25,0.25,0.25,0.25)
q0 = c(0.25,0.25,0.25,0.25)

values_per_layer = 5
simulated_values_per_layer = 30
simulations_per_est = 100

max_depth = 4
improvement_threshold = -10000

p0_estimate = c(0.25,0.25,0.25,0.25)
q0_estimate = c(0.25,0.25,0.25,0.25)
generations_to_sim = 50

##########################
#creating a synthetic dataset of the 3 measures
#nx3 matrix created using 'create_synthetic_observed_data'
real_values = create_synthetic_observed_data(N, M, migration, recombination_rates,
                                             generations_to_sim, p0, q0)
# doing the sequential MC
full_mat = sequential_mc_exp_migration(real_values, recombination_rates, values_per_layer, 
                         simulated_values_per_layer, simulations_per_est, max_depth, improvement_threshold,
                         p0_estimate, q0_estimate, generations_to_sim )
#plotting and creating a dataframe with the results
final_df = df_and_plots(full_mat, values_per_layer, simulated_values_per_layer, 1, FALSE)

#if you want to save the dataframe
path = paste("C:/Users/angus/OneDrive/Desktop/urop/ABC/tests_csv/test_", as.character(7), '.csv')
write.csv(x=final_df, file= path)

###############################################################################