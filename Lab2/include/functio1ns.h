#pragma once
#include <stdio.h>
#include <gsl/gsl_rng.h>

// Define the result_t struct
typedef struct {
    double integral;
    double error;
} result_t;

// Define the result_t5 struct
typedef struct{
    double weight;
    double function_value;
    int accepted;
} result_t5;

// Function declarations
double MCMC_step_displace_all(double **N, int *type, int* type_change, double *N_XY, double *E_ref, double* N_A, double k_B, double T, int n_atoms, gsl_rng *rng);

double average( double *v,unsigned int len);
double variance(double *v,unsigned int len);

double autocorrelation(double *data, int data_len, int time_lag_ind);
double block_average(double *data,int data_len,int block_size);
/* ***************************************
*
* Perform the calculation of the
* unormalized probability function, i.e.,
* the weight used in the MCMC routine
*
* Parameters
* ----------
*  x - current walker position, size = 3
*
* Returns
* -------
* The weight of the current position
*
* ***************************************/
double weight(int *x);


/* ***************************************
*
* Perform the calculation of the
* function to be sampled
*
* Parameters
* ----------
*  x - current walker position, size = 3
*
* Returns
* -------
* The function value
*
* ***************************************/
double function(double **neighbors, int *type, double *N_XY, double *E_ref, double *N_A, int n_atoms);


/* ***************************************
*
* Perform the calculation of the
* function to be sampled
*
* - Use gsl_rng_uniform to displace the
*   walker.
* - The walker should be displaced
*   in the order x,y,z.
* - The draw for the acceptance condition
*   should be done after the draws
*   for the displacing of the walker
*
*
* Parameters
* ----------
*  x - current walker position, size = 3
*  delta - step size
*  k - GSL random number generator object
*
* Returns
* -------
* - The function should update the x parameter
*   to reflect if the move was accepted or
*   rejected.
* - result should contain the probability
*   of the "exiting" x parameter.
* - result should should contain the function
*   value of the "exiting" x parameter.
* - result should be 1 if the move was accepted
*   and 0 if it was rejected.
*
* ***************************************/
