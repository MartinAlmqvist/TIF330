#include <math.h>
#include <gsl/gsl_rng.h>
#include "functions.h"
#define PI 3.141592653589


double function(
                double **neighbors, int *type, double *N_XY, double *E_ref, double *N_A, int n_atoms)
{
    double E = 0.;
    double N_CC = 0., N_ZZ = 0., N_CZ = 0.;
    int b;
    
    double n_A = n_atoms;
    
    for(int a = 0; a < n_atoms; ++a){
    	for(int i = 0; i < 8; ++i){
    		b = (int) neighbors[a][i] + n_atoms;
    		if(type[a] == 0 && type[b] == 0){
    			N_CC += 1.;
    		}
    		else if(type[a] == 1 && type[b] == 1){
    			N_ZZ += 1.;
    		}
    		else{
    			N_CZ += 1.;
    		}
    	
    	}
    	n_A -= type[a]; // Calculate amout of A atoms in the A lattice by subtracting B
    }
    N_XY[0] = N_CC;
    N_XY[1] = N_ZZ;
    N_XY[2] = N_CZ;
    N_XY[3] = n_A;
    // &N_A = n_A;
    
    for(int i = 0; i<3; ++i){
    	E += N_XY[i] * E_ref[i];
    }
    
    return E;
}


double MCMC_step_displace_all(double **N, int *type, int *type_change, double *N_XY, double *E_ref, double *N_A, double k_B, double T, int n_atoms, gsl_rng *rng)
{
	double E_change, E;
	int temp;
	for(int i = 0; i<2*n_atoms; ++i){
		type_change[i] = type[i];
	}
    
	int n_a, n_b;
	int change = 0;
	
	E = function(N, type, N_XY, E_ref, N_A, n_atoms);
	
	do{
    	//randomize atoms until we change two of them (Displace the walker)
		n_a = (int) (gsl_rng_uniform(rng) * 2 * n_atoms);
	    	n_b = (int) (gsl_rng_uniform(rng) * 2 * n_atoms);
	    	if(type_change[n_a] != type_change[n_b]){
	    		temp = type_change[n_a];
	    		type_change[n_a] = type_change[n_b];
	    		type_change[n_b] = temp;
	    		change = 1;
	    	}
	}while(change == 0);
    
	E_change = function(N, type_change, N_XY, E_ref, N_A, n_atoms);

	double delta = E_change - E; 
	double random_Number = gsl_rng_uniform(rng);
	double ratio = exp(-delta / (k_B * T));
    
	//Acceptance
	if(delta <= 0 || ratio >= random_Number){
		E = function(N, type_change, N_XY, E_ref, N_A, n_atoms);
        
		// Update the walker's position
		type[n_a] = type_change[n_a];
		type[n_b] = type_change[n_b];   
    }
    return E;
}


