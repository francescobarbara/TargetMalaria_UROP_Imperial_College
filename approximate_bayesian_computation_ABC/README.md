
**'all_fncs4.R'** has all the functions you need, just save **'simulation_3.cpp'** in your computer and modify the path in the first line of 'all_fncs4.R'.
I suggest running the whole script once, becuase all the functions are interconnected.

If you want to use the l^2 measure, then you need to use **'likelihood_test'**

If you want to use the dissimilarity function, then you need to run **'compare'**

If you want to try the sequential MC, just modify the example script **'sequential_mc_example.R'**

**'all_fncs2.R'** has a slightly difference than 'all_fncs4.R.
The difference is in the way points are sampled in the MC algo.
If the migration rate in not in [0, 0.5], then you resample,
this will give almost never values with migr = 0, which isn't great

**'old_compare.R'** was an earlier approach for a dissimilarity function, it uses squared distances
and returns the mean squared error

**'old_austin.R'** is the fixed dissimilarity function, essentially is compare function with the final line modified

**specific_estimate_dissimilarity.R**
outputs for the dissimilarity of
a specifically given estimate.
You pick the true parameters, create the synthetic real data, then choose the 3 estimates
and see what the value of the dissimilarity function is

**multiple_tests_script.R** Here I want to do multiple (14 in my case) sequential Montecarlos
for different true parameters
 The script was run overnight and it is set up in a way such that it
 saves each dataframe once that specific run of the MC algo is finished
 The reults are stored in **test_table .csv**
 
**uniform_assumption_test.R**  Here I am trying to see how much the uniform approximation holds, and the consequent
 Normal test that I introduced to estimate the probability that 'real data can come
from an estimated triplet'
Essentially the question is, how many points do I need to simulate to get an empirical cdf that behaves well?
 To this aim, for each i in [1,...,no_sim] , I store the value of the dissimilarity easure
 After all the simulatons, I check to see whether the variance is the same
 as the assumed one under the Uniform assumpution

