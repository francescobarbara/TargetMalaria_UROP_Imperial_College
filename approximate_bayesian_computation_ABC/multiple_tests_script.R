# Here I want to do multiple (14 in my case) sequential Montecarlos
# for different true parameters
# The script was run overnight and it is set up in a way such that it
# saves each dataframe once that specific run of the MC algo is finished

generations_to_sim = 50

p0 = c(0.25,0.25,0.25,0.25)
q0 = c(0.25,0.25,0.25,0.25)

values_per_layer = 5
simulated_values_per_layer = 30
simulations_per_est = 100

max_depth = 6
improvement_threshold = -10000

p0_estimate = c(0.25,0.25,0.25,0.25)
q0_estimate = c(0.25,0.25,0.25,0.25)
generations_to_sim = 50

##################################

triplets1 = c(1000, 1500, 4000, 500, 5000, 1000, 1250, 1300, 600, 2000, 
              1500, 1100, 3400, 1500, 800, 2600, 6000, 1000, 600, 2400)
triplets2 = c(1000, 1500, 4000, 500, 1000, 500, 250, 1700, 2000, 600, 
              300,  3000, 1000, 1500, 7000,1000, 2000, 1000, 1500,1400)
triplets3 = c(0.1,  0,    0.22, 0.5, 0.32, 0.11, 0.05, 0.17, 0.22, 0.02, 
              0, 0.2, 0.5, 0.05, 0.1, 0.2, 0.5, 0, 0.05, 0.11, 0.22)

for (i in 14:length(triplets1)){
  N = triplets1[i]
  M = triplets2[i]
  m = triplets3[i]
  migration = matrix(c(1-m, m, m, 1-m), nrow = 2)
  recombination_rates = sample(c(0.05, 0.1, 0.15,0.2, 0.25,0.3, 0.35, 0.4, 0.45, 0.5),
                               size = 1000, replace = TRUE)
  
  # creating synthetic data for the parameters
  real_values = create_synthetic_observed_data(N, M, migration, recombination_rates,
                                               generations_to_sim, p0, q0)
  # doing the sequential MC algo
  full_mat = sequential_mc_exp_migration(real_values, recombination_rates, values_per_layer, 
                                         simulated_values_per_layer, simulations_per_est, max_depth, improvement_threshold,
                                         p0_estimate, q0_estimate, generations_to_sim )
  # saving the dataframe
  # REMARK: CHANGE THE PATH IN df_and_no_plots so that it saves the df in the right place
  final_df = df_and_no_plots(full_mat, values_per_layer, simulated_values_per_layer, i, TRUE)
}

#Now I want to create a new dataframe that gives us the confidence intervals
# and the point estimates for the 3 params
# the results are 'test_table.csv'
estimates = matrix(0, 11, 12)

for (i in 1:11){
  # importing the saves dataframes, created with the loop above
  df_full = read.csv(paste("C:/Users/angus/OneDrive/Desktop/urop/ABC/100_tests_csv/test_", as.character(i), '.csv'))
  
  # keeping the final layer, if you did more then 6 layers, change the number below accordingly
  df = df_full[df_full$layer == 6,]
  
  estimates[i,] = c(triplets1[i], triplets2[i], triplets3[i],
                    mean(df$population_1), mean(df$population_2), mean(df$migration),
                    min(df$population_1), max(df$population_1), min(df$population_2),
                    max(df$population_2), min(df$migration), max(df$migration) )
}

new_df = data.frame(estimates)
colnames(new_df) = c('N_real', 'M_real', 'migr_real', 'N_est', 'M_est', 'migr_est',
                     'N_min', 'N_max', 'M_min', 'M_max', 'migr_min', 'migr_max')

#saving this new dataframe
path = paste("C:/Users/angus/OneDrive/Desktop/urop/ABC/100_tests_csv/test_table", '.csv')
write.csv(x=new_df, file= path)
