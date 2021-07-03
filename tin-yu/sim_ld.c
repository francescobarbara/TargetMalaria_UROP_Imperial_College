// INDIVIDUAL BASED LD SIMULATOR
// TO BE COMPILED. R CMD SHLIB sim_ld.c
// 002/07/2021
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

// SAMPLE PARENT FROM 0, 1, 2, ..., (N-1). UNIFORMLY. 
int sample_parent(int N)
{
	double temp=runif(0, (double) N);
	return floor(temp);
}

// A FUNCTION TO CALCULATE r2
SEXP cal_r2(SEXP dat)
{
	double p[4]={0, 0, 0, 0};
	int nrow=INTEGER(GET_DIM(dat))[0];
	Rbyte *h_dat=RAW(dat); 	Rbyte temp;
	for (int i=0; i<nrow; i++)
	{
		temp=2*h_dat[i]+h_dat[nrow+i];
		p[temp]=p[temp]+1;
	}
	double D=p[0]*p[3]-p[1]*p[2];
	return ScalarReal(D*D/((p[0]+p[1])*(p[2]+p[3])*(p[0]+p[2])*(p[1]+p[3])));
}

// THE MAIN FUNCTION CALLABLE FROM R .CALL()
// N IS A VECTOR OF POPULATION SIZES (AND DETERMINES HOW MANY GENERATIONS TO RUN FORWARD IN TIME)
// c IS THE RECOMBINATION RATE. initial IS A MATRIX SPECIFIED BY THE INITIAL HAPLOTYPE FREQ
SEXP sim_ld(SEXP initial, SEXP N, SEXP c)
{
	// GET THE RANDOM NUMBER GENERATOR STATE FROM R. 
	GetRNGstate();
	// NUMBER OF GENERATIONS
	int t=length(N);
	// result IS A BIG LIST CONTAINING ALL THE ALLELIC CONFIGURATIONS OVER TIME. 
	// CODED IN 8-BIT RAW (char) INTEGER
	SEXP result=PROTECT(allocVector(VECSXP, t+2));
	// INITIALISE. N_parent IS THE POPULATION SIZE OF THE PARENTAL GENERATION. WILL CHANGE OVER TIME. 
	int N_parent=INTEGER(GET_DIM(initial))[0]/2;
	SET_VECTOR_ELT(result, 0, initial);
	// INITIALSE r2 VECTOR. THIS WILL BE APPENDED TO result AT THE END OF THE PROGRAM. 
	SEXP r2=PROTECT(allocVector(REALSXP, t+1));
	double *h_r2=REAL(r2);
	h_r2[0]=asReal(cal_r2(VECTOR_ELT(result, 0)));
	// POINTER TO THE N (POPULATION SIZE) VECTOR
	int *h_N; h_N=INTEGER(N);
	// CREATE SOME VARIABLES TO DENOTE THE TWO PARENTS AND WHICH CHRONOSOMES THEY CONTRIBUTE
	// p1_l2 MEANS "PARENT 1, LOCUS 2"
	int parent1; int parent2;
	int p1_l1=0; int p1_l2=0; int p2_l1=0; int p2_l2=0;
	Rbyte *new_gen; Rbyte *parent_gen;
	// PROPAGATION
	for (int i=0; i<t; i++)
	{
		// SET THE MATRIX FOR THE NEXT GEN
		SET_VECTOR_ELT(result, i+1, PROTECT(allocMatrix(RAWSXP, 2*h_N[i], 2)));
		// POINTERS TO THE OFFSPRING (NEXT) AND PARENTAL GENERATIONS. PURE PERFORMANCE THING. 
		new_gen=RAW(VECTOR_ELT(result, i+1));
		parent_gen=RAW(VECTOR_ELT(result, i));
		// FOR EACH NEW INDIVIDUAL
		for (int k=0; k<h_N[i]; k++)
		{
			// SAMPLE THE TWO PARENTS
			parent1=sample_parent(N_parent);
			parent2=sample_parent(N_parent);
			// CHOOSE ONE CHROMOSOME PER PARENT
			p1_l1=runif(0, 1)<0.5; p1_l2=p1_l1;
			p2_l1=runif(0, 1)<0.5; p2_l2=p2_l1;
			// IF THERE IS RECOMBINATION THEN SWAP LOCUS 2
			if (runif(0, 1)<asReal(c)) {p1_l2=1-p1_l2;}
			if (runif(0, 1)<asReal(c)) {p2_l2=1-p2_l2;}
			// FILL IN NEW HAPLOTYPES. 2*k AND 2*k+1 ROWS BELONG TO THE SAME INDIVIDUAL
			// NOTE THE TWO POPULATION SIZES. 
			new_gen[2*k]=parent_gen[2*parent1+p1_l1];
			new_gen[2*h_N[i]+2*k]=parent_gen[2*N_parent+2*parent1+p1_l2];
			new_gen[2*k+1]=parent_gen[2*parent2+p2_l1];
			new_gen[2*h_N[i]+2*k+1]=parent_gen[2*N_parent+2*parent2+p2_l2];						
		}
		// CALCULATE r2
		h_r2[i+1]=asReal(cal_r2(VECTOR_ELT(result, i+1)));
		// FINALLY UPDATE N_parent;
		N_parent=INTEGER(N)[i];
	}
	// APPEND THE r2 MATRIX TO THE result BIG LIST, AND THEN RETURN EVERYTHING TO R. 
	SET_VECTOR_ELT(result, t+1, r2);
	PutRNGstate();
	UNPROTECT(t+2);
	return result;
}

