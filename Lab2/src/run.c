#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "tools.h"
#include "lattice.h"
#include "functions.h"

int
run(
    int argc, double T,
    char *argv[]
   )
{
	// TODO add read xyz (removed due to fscanf bug)
	/*
		--- Part 2 ---
	*/
	
	int N = 10; 						// Cells in each direction
	int n_atoms = N * N * N; 				// Number of atoms in each sublattice, 2 per cell for a sc-lattice
	
	
	// Allocate arrays
	double **neighbors_A = create_2D_array(n_atoms, 8);		// Neighbor array for atoms in sublattice A
	int *types = (int*) malloc(2 * n_atoms * sizeof(int));	// Atom material in the lattice, atoms in both sublattices
	int *type_change = (int*) malloc(2 * n_atoms * sizeof(int));	// Array to store change
	double * E_ref = create_1D_array(3);				// Binding energies
	double *N_XY = create_1D_array(4);
	double *N_A = create_1D_array(1);
	
	// Constants
	double k_B = 8.617e-5;
	E_ref[0] = -0.436, E_ref[1] = -0.113, E_ref[2] = -0.294; 
	
	// To denote atomic type
	int Cu = 0;
    	int Zn = 1;
    	
	// Run variables
    	// double T = 400.;
    	unsigned int N_eq = 1e6;		// # of equilibratiion steps
	unsigned int N_steps = 4e6;		// # of steps after equilibration
	double E;
	
	// Define a filename with temperature
	char filename[50] = "MC_task2_short_";
	char filetype[] = ".csv";	
	sprintf(filename + strlen(filename), "%d", (int) T);
	strcat(filename, filetype);
    	
    	// Initiate neighbors and types
    	idx_nearest(neighbors_A, Cu, n_atoms, N);
    	
    	// Initiate random number for MCMC acceptance
    	const gsl_rng_type * rng_T;
    	gsl_rng * rng;
    	gsl_rng_env_setup();
	rng_T = gsl_rng_default;
    	rng = gsl_rng_alloc(rng_T);
    	time_t seed = time(NULL);
    	gsl_rng_set(rng, seed);
    	
    	
    	// Initate lattice with total order
	for (int i = 0; i <n_atoms; ++i){
		types[i] = Cu;
		types[i + n_atoms] = Zn;
	}
 	
 	printf("\n --- Equilibration --- \n");
	FILE *fp = fopen(filename, "w+");
	fprintf(fp, "Temp \t Eq steps \t Production steps\n");
	fprintf(fp, "%f \t %i \t %i \n", T, N_eq, N_steps);
	fprintf(fp, "E \t N_CC \t N_ZZ \t N_CZ \t N_A \n");	

	// Equilibration
	for(int i = 0; i < N_eq; i++){
		E = MCMC_step_displace_all(neighbors_A, types, type_change, N_XY, E_ref, N_A, k_B, T, n_atoms, rng);
		fprintf(fp, "%f \t %f \t %f \t %f \t %f \t %f \n", E, N_XY[0], N_XY[1], N_XY[2], N_XY[3], N_A[0]);
		if(i % 100000 == 0){
			printf("%i\n", i);
		}
	}
	
	printf("\n --- Equilibration done, starting production run --- \n");
	
	// Production run
	for(int i = 0; i < N_steps; i++){
		E = MCMC_step_displace_all(neighbors_A, types, type_change, N_XY, E_ref, N_A, k_B, T, n_atoms, rng);
		fprintf(fp, "%f \t %f \t %f \t %f \t %f \t %f \n", E, N_XY[0], N_XY[1], N_XY[2], N_XY[3], N_A[0]);
		if(i % 100000 == 0){
			printf("%i\n", i);
		}
	}
	
	fclose(fp);
	destroy_2D_array(neighbors_A, n_atoms);
	free(types);
	free(type_change);
	free(E_ref);
	free(N_XY);
	
	return 0;
}
