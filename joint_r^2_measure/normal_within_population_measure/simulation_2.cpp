#include <R.h>
#include <Rcpp.h>
#include <cmath>    /* needed for floor() */
using namespace Rcpp;
using namespace std;



// [[Rcpp::export]]
IntegerMatrix constrained_matrix(NumericVector p, int n){

    /* inputs: p = c(p1, p2, p3, p4) in [0,1], represent the (0,0), (0,1), (1,0) frequencies you want the simulated matrix to have
               n = number of diploid individuals
        output: 2x2n matrix with the required haplotype frequencies
    */

    // q0 is the number of alleles that will have form {0, 0}
    int q0 = 2*n*p[0];
    int q1 = 2*n*p[1];
    int q2 = 2*n*p[2];
    int q3 = 2*n - q0 - q1 - q2;

    // possible warning if n is too low
    if(abs(q3 - floor(2*n*(1 - p[0] - p[1] -p[2])))/ (2*n) > 0.1){
        Rcout << "imprecise approximation due to discretization error \n";
    }

    // sampling without replacement alleles {0, ..., 2n -1} 
    IntegerVector v = seq(0, 2*n -1);
    IntegerVector parents = sample(v, 2*n, false);
    
    // initializing empty matrix, which will contain the output 
    IntegerMatrix m(2 , 2*n);

    // unordered_map<int, IntegerVector> umap; //maps column index of the future matrix to the 2 alleles it will contain

    // giving the first q0 columns alleles (0, 0), then up to q0 + q1 alleles (0, 1), etc.
    int count = 0;
    for (IntegerVector::iterator i = parents.begin(); i != parents.end(); ++i){
        count += 1;
        if (count <= q0){
            m(_, *i) = IntegerVector::create(0, 0);
        }
        else if (count <= q0 + q1){
            m(_, *i) = IntegerVector::create(0, 1);    
        }
        else if (count <= q0 + q1 + q2){
            m(_, *i) = IntegerVector::create(1, 0);    
        }
        else{
            m(_, *i) = IntegerVector::create(1, 1);    
        }
    }
    
    return m;   
    
}





// [[Rcpp::export]]
List calculate_r_sq(IntegerMatrix m){

    /* inputs: a 2xn binary matrix
       output: list containing: [haplotype freq., raw Ld, r^2] for the matrix 
    */

    int n = m.ncol();
    NumericVector p (4); // vector of frequencies

    //calculating counts
    IntegerMatrix::iterator i = m.begin();
    while(i < m.end()) {
        if ( (*i == 0) & (*(i+1) == 0)) p[0] += 1;
        else if ((*i == 0) & (*(i+1) == 1)) p[1] += 1;
        else if ((*i == 1) & (*(i+1) == 0)) p[2] += 1;
        else p[3] += 1;
        ++i;
        ++i;
    }

    // transforming counts into frequencies
    p =  p/n;


    double raw_LD = p[0]*p[3] - p[1]*p[2];
    double r_sq = (pow(raw_LD, 2))/((p[0] + p[1])*(p[2] + p[3])*(p[0] + p[2])*(p[1] + p[3]));
    List out = List::create(p, raw_LD, r_sq);
    return out;
}







// [[Rcpp::export]]
IntegerMatrix generate_child_matrix(IntegerMatrix m, double c, int size){

    /* inputs: m = 2x2*n parent matrix
               c = recombination rate in [0, 0.5]
               size = number of individuals you want in child generation
        output: 2x2*size child matrix
    */

    // 2*size vector of uniform r.v.
    NumericVector unif = Rcpp::runif(2*size, 0.0, 1.0 );  //unif[i] will tell us if allele is coming from parents[i], or needs to recombine
    
    // number of parents
    int n = m.ncol()/2;  

    // parents[i] is the parent of the i-th allele in the output matrix; we are obvs sampling with replacement
    IntegerVector v = seq(0, n -1);
    IntegerVector parents = sample(v, 2*size, true);

    // this is our output matrix    
    IntegerMatrix children (2, 2*size);

    // pointer at the start of the children matrix
    IntegerVector::iterator i = children.begin();

    // pointer at the start of the sampled parents vector
    IntegerVector::iterator j = parents.begin();

    // main loop, using pointers for better performance
    for(NumericVector::iterator k = unif.begin(); k != unif.end(); ++k){

        // recombination occurs Ab  (wlog allele 0 of parent is AB, allele 1 of parent is ab)
        if(*k < c/2.0){
            *i = m(0, 2*(*j));
            ++i;
            *i = m(1, 2*(*j) + 1);
            ++i;
            ++j;
        }

        // recombination occurs aB
        else if(*k < c){
            *i = m(0, 2*(*j) + 1);
            ++i;
            *i = m(1, 2*(*j));
            ++i;
            ++j;
        }

        // no recombination AB
        else if(*k < (1 + c)/2.0){
            *i = m(0, 2*(*j));
            ++i;
            *i = m(1, 2*(*j));
            ++i;
            ++j;
        }

        // no recombination ab
        else{
            *i = m(0, 2*(*j) + 1);
            ++i;
            *i = m(1, 2*(*j) + 1);
            ++i;
            ++j;
        }
        
    }

    return children;
}







// [[Rcpp::export]]
IntegerMatrix merge_on_indices (IntegerMatrix m1, IntegerMatrix m2, IntegerVector index1, IntegerVector index2){

    /* inputs: m1, m2 matrices to merge
       inputs: index1 indicates which individuals from m1 to use
               for ex. if 2 is in index1, then the two gametes corresponding to individual 2 (i.e  m[:, 2*2], m[:, 2*2 +1])
                                          will be put in the output matrix
        output: a matrix which takes the individuals from m1, m2 according to index1, index2
    */     

    int len1 = index1.length();
    int len2 = index2.length();
    int n = len1 + len2;

    // we have n individuals, so matrix must have ncol = 2*n
    IntegerMatrix out (2, 2*n);

    // using pointers for better performance, first 2*len1 columns will come from m1
    // last 2*len2 columns will come from m2 
    IntegerMatrix::iterator i = out.begin();
    for (IntegerVector::iterator j = index1.begin(); j != index1.end(); ++j){
        *i = m1(0, 2*(*j));
        ++i;
        *i = m1(1, 2*(*j));
        ++i;
        *i = m1(0, 2*(*j) + 1);
        ++i;
        *i = m1(1, 2*(*j) + 1);
        ++i;
    }
    for(IntegerVector::iterator j = index2.begin(); j != index2.end(); ++j){
        *i = m2(0, 2*(*j));
        ++i;
        *i = m2(1, 2*(*j));
        ++i;
        *i = m2(0, 2*(*j) + 1);
        ++i;
        *i = m2(1, 2*(*j) + 1);
        ++i;
    }
    return out;    
}










// [[Rcpp::export]]
IntegerVector migration_vec(IntegerVector N, NumericMatrix m, int index1, int index2){

    /* inputs:  N = vector of population sizes for generations 1, 2, 3,...
                m = migration matrix
                index1, index2 indicate which entry of m you want. i.e. migration rate from population "index1" (which has size N) to population "index2"
        output: a vector like N, that indicates at each generation the number of individuals going 
                from population "index1" to population "index2"
    */

    IntegerVector out (N.length());
    double val = m(index1, index2);

    IntegerVector::iterator i = N.begin();
    for (IntegerVector::iterator j = out.begin(); j != out.end(); ++j){
        *j = floor(val * (*i));
        ++i;
    }


    return out;
}







// [[Rcpp::export]]
List simulation(NumericVector p0, NumericVector q0, IntegerVector N, IntegerVector M, double c, NumericMatrix migration){

    /*

    inputs: p0: 4x1 vector; initial haplotypes frequencies at generation 0 for population 1
    inputs: q0: 4x1 vector; initial haplotypes frequencies at generation 0 for population 2
    inputs: N: vector representing the size of population 1 at generation i, BEFORE MIGRATION OCCURS
    inputs: M: vector representing the size of population 2 at generation i, BEFORE MIGRATION OCCURS
            M and N should have the same size, which is the number of generations to simulate
    inputs: c is the recombination rate for the two loci
    inputs: migration: FIXED migration matrix; migration[0,1] represents the proportion of population 1 that migrates
                       to population 2 at each generation
  
    output: list of size 6 containing: simulated haplotypes frequencies for the two populations for all generations
                                     simulated D and r^2 values for the two populations for all generations

    */

    // number of generations
    int t = N.length();

    // initializing matrices of haplotype frequencies (to be filled during the loop)
    NumericMatrix p (4, t);
    NumericMatrix q (4, t);

    // initializing vectors of raw LD measures (length = number of generations)
    NumericVector D1 (t);
    NumericVector D2 (t);

    // initializing vectors of r^2 measures (length = number of generations)
    NumericVector r_sq1 (t);
    NumericVector r_sq2 (t);

    // calculating the migration flows at each generation
    // ex. size2[i] indicates how many individuals will migrate from pop1 to pop2 between generation i and generation i+1
    IntegerVector size1 = migration_vec(N, migration, 0, 0);   //stay in pop1
    IntegerVector size2 = N - size1;
    IntegerVector size3 = migration_vec(M, migration, 1, 1);   //stay in pop2
    IntegerVector size4 = M - size3;
    /*
    Rcout << size3 << "\n";
    Rcout << size1 << "\n";
    */

    // initialization: old_matrix1 will contain the parent generation 1 through the loop, new_matrix will containg the child generation
    IntegerMatrix old_matrix1 = constrained_matrix(p0, N[0]);
    IntegerMatrix old_matrix2 = constrained_matrix(q0, M[0]);
    IntegerMatrix new_matrix1( 2 , 3 );  // this is just a random initialization, to avoid having multiple objects 
    IntegerMatrix new_matrix2( 2 , 2);   // with same name but different pointers (which is what happens if we define new_matrix within the loop,
                                         // at every iteration, we would define a new object; in this way at each iteration we just modifiy this object)
    
    /*
    Rcout << "before loop" << "\n";
    Rcout << old_matrix1 << "\n";
    Rcout << old_matrix2 << "\n";
    */
    
    // calculating the 3 statistics for generation 0
    List stats_pop1 = calculate_r_sq(old_matrix1);
    List stats_pop2 = calculate_r_sq(old_matrix2);
     
    // stats_pop1[0];

    // updating the 3 statistics for generation 0
    p(_, 0) = static_cast<NumericVector>(stats_pop1[0]);
    q(_, 0) = static_cast<NumericVector>(stats_pop2[0]);  

    D1[0] = stats_pop1[1];
    D2[0] = stats_pop2[1];

    r_sq1[0] = stats_pop1[2];
    r_sq2[0] = stats_pop2[2];
    

    for (int i = 1; i < t; ++i){
        
        // sampling which individuals will stay in pop1, which will migrate to pop2
        IntegerVector v = seq(0, N[i - 1] - 1 );
        IntegerVector temp1 = sample(v, v.length(), false);
        IntegerVector pop_1_stays = temp1[Range(0, size1[i- 1] - 1)];
        IntegerVector pop_1_migrates = temp1[Range(size1[i - 1], v.length() - 1)];

        // same thing for pop2
        v = seq(0, M[i - 1] - 1 );                                             
        IntegerVector temp2 = sample(v, v.length(), false);
        IntegerVector pop_2_stays = temp2[Range(0, size3[i - 1] - 1)];
        IntegerVector pop_2_migrates = temp2[Range(size3[i - 1], v.length() - 1)];

        /*
        Rcout << "migration vecs" << "\n";
        Rcout << pop_1_stays << "\n";
        Rcout << pop_1_migrates << "\n";
        Rcout << pop_2_stays << "\n";
        Rcout << pop_2_migrates << "\n";

        Rcout << "old matrix" << "\n";
        Rcout << old_matrix1 << "\n";
        Rcout << old_matrix2 << "\n";
        */

        // creating the new populations after migration occurs, i.e. putting pop_1_stays and pop_2_migrates together etc.
        new_matrix1 = merge_on_indices(old_matrix1, old_matrix2, pop_1_stays, pop_2_migrates);
        new_matrix2 = merge_on_indices(old_matrix1, old_matrix2, pop_1_migrates, pop_2_stays);

        /*
        Rcout << "After migration, before reproduction" << "\n";
        Rcout << new_matrix1 << "\n";
        Rcout << new_matrix2 << "\n";
        */

        // making the new populations reproduce, size of the populations are specified by N[i], M[i]
        new_matrix1 = generate_child_matrix(new_matrix1, c, N[i]);    
        new_matrix2 = generate_child_matrix(new_matrix2, c, M[i]);

        /*
        Rcout << "After migration, after reproduction" << "\n";
        Rcout << new_matrix1 << "\n";
        Rcout << new_matrix2 << "\n";
        */

        // updating old matrix to generation i
        old_matrix1 = new_matrix1;
        old_matrix2 = new_matrix2;

        // calculating 3 statistics for generation i
        List stats_pop1 = calculate_r_sq(old_matrix1);
        List stats_pop2 = calculate_r_sq(old_matrix2);

        // updating the 3 statistics in the output vectors
        p(_,i) = static_cast<NumericVector>(stats_pop1[0]);
        q(_,i) = static_cast<NumericVector>(stats_pop2[0]);

        D1[i] = stats_pop1[1];
        D2[i] = stats_pop2[1];

        r_sq1[i] = stats_pop1[2];
        r_sq2[i] = stats_pop2[2];


    }
    
    List out = List::create(p, q, D1, D2, r_sq1, r_sq2);
    return out;

}












// [[Rcpp::export]]
List average_r_sq_individual_based(NumericVector p0, NumericVector q0, IntegerVector N, IntegerVector M, double c, NumericMatrix migration, int simulations){
    
    // This function simply takes the average r^2 over a given number of simulations
    // Returns a list of the two average r^2 for pop 1 and pop 2
    
    
    NumericVector r_sq_1 (N.length());
    NumericVector r_sq_2 (N.length());


    for (int n = 0; n < simulations; ++n){

        List l_temp = simulation(p0, q0, N, M, c, migration);
        NumericVector out1 = l_temp[4];
        NumericVector out2 = l_temp[5];    
           
        r_sq_1 += out1;
        r_sq_2 += out2;
 
        
    }

    r_sq_1 = r_sq_1 / simulations;
    r_sq_2 = r_sq_2 / simulations;

    List out = List::create(r_sq_1, r_sq_2);

    return out;
}


// [[Rcpp::export]]
double average_r_sq_at_equilibrium(NumericVector p0, NumericVector q0, IntegerVector N, IntegerVector M, double c, NumericMatrix migration, int simulations, int start_gen){
    
    /* Given the usual starting parameters plus "start_gen", it computes r^2 for all generations for "simulations" times
       Then it takes the avearage of values from "start_gen" to the final generation and returns the average over all simulations
       Essentially, we assume we are already in equilibrium at "start_gen", say 50, so if we simulate 100 generations, and take the average
       over generations 51,...,100, at each iteration we get 50 data points instead of 1! 
       The output is a single value, i.e. the equilibrium value of r^2 for the initial parameters
    */
    
    int n = N.length();
    int m = n-start_gen;
    
    IntegerVector v = seq(start_gen - 1, n-1);
    double r_sq = 0;    

    for (int n = 0; n < simulations; ++n){
        
        List l_temp = simulation(p0, q0, N, M, c, migration);
        NumericVector out1 = l_temp[4];
           
        for (int i = 0; i < m; ++i){
           
            r_sq += out1[start_gen + i - 1];
        }            
    }
    
    r_sq = r_sq / (simulations*(n - start_gen));
    return r_sq;
}