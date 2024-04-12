#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <complex.h>

#include "tools.h"

void
elementwise_addition(
                     double *res,
                     double *v1,
                     double *v2,
                     unsigned int len
                    )
{
	for(int i=0; i<len; i++)
    	{
        	res[i] = v1[i] + v2[i];
    	}
}

void
elementwise_addition_matrix(
                     double **res,
                     double **A,
                     double **B,
                     unsigned int n,
                     unsigned int m
                    )
{
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<m; ++j)
		{
			res[i][j] = A[i][j] + B[i][j];
		}
	}
}

void
elementwise_multiplication(
                           double *res,
                           double *v1,
                           double *v2,
                           unsigned int len
                          )
{
	for(int i=0; i<len; i++)
    	{
        	res[i] = v1[i] * v2[i];
    	}
}

void
elementwise_multiplication_matrix(
                     double **res,
                     double **A,
                     double **B,
                     unsigned int n,
                     unsigned int m
                    )
{
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<m; ++j)
		{
			res[i][j] = A[i][j] * B[i][j];
		}
	}
}

void
addition_with_constant(
                       double *res,
                       double *v,
                       double constant,
                       unsigned int len
                       )
{
	for(int i=0; i<len; i++)
    	{
        	res[i] = v[i] + constant;
    	}
}

void
addition_with_constant_matrix(
                       double **res,
                       double **A,
                       double constant,
                       unsigned int n,
                       unsigned int m
                       )
{
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<m; ++j)
		{
			res[i][j] = A[i][j] + constant;
		}
	}
}

void
multiplication_with_constant(
                             double *res,
                             double *v,
                             double constant,
                             unsigned int len
                             )
{
	for(int i=0; i<len; i++)
    	{
        	res[i] = v[i] * constant;
    	}
}

void
multiplication_with_constant_matrix(
                             double **res,
                             double **A,
                             double constant,
                             unsigned int n,
                             unsigned int m
                             )
{
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<m; ++j)
		{
			res[i][j] = A[i][j] * constant;
		}
	}
}

double
dot_product(
            double *v1,
            double *v2,
            unsigned int len
           )
{
	double sum = 0;
	for(int i=0; i<len; i++)
    	{
        	sum += v1[i] * v2[i];
    	}
    	return sum;
}

double *
create_1D_array(
                unsigned int n
               )
{
	double* array = (double*)malloc(n * sizeof(double));
	if(array == NULL){
		fprintf(stderr, "Malloc failed\n");
		free(array);
		exit(EXIT_FAILURE);
	}
	return array;
}

double **
create_2D_array(
                unsigned int n,
                unsigned int m
               )
{
	double** array = (double**)malloc(n * sizeof(double*));
	if(array == NULL){
		fprintf(stderr, "Malloc failed\n");
		exit(EXIT_FAILURE);
	}
	
	//Continuous allocation of memory
	double* data = (double*)malloc(n * m * sizeof(double));
	if(data == NULL){
		fprintf(stderr, "Malloc failed\n");
		free(array);
		exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n; ++i){
    		array[i] = data + i*m;
	}
	return array;
}

void
destroy_2D_array(
                 double **array,
                 unsigned int n
                )

{
        free(array[0]);  // Free continous allocation
	free(array);
}

void
matrix_vector_multiplication(
                             double *result,
                             double **A,
                             double *b,
                             unsigned int n,
                             unsigned int m
                            )
{
	for(int i=0; i<n; ++i)
	{
		result[i] = 0.;
		for(int j=0; j<m; ++j)
		{
			result[i] += A[i][j] * b[j];
		}
	}
}

void
matrix_matrix_multiplication(
                             double **result,
                             double **A,
                             double **B,
                             unsigned int n,
                             unsigned int m,
                             unsigned int k
                            )
{
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<k; ++j)
		{
			result[i][j] = 0.;
			for(int l=0; l<m; ++l)
			{
				result[i][j] += A[i][l] * B[l][j];
			}
		}
	}
	
}

double
vector_norm(
            double *v1,
            unsigned int len
           )
{
	double norm = dot_product(v1, v1, len);
	norm = sqrt(fabs(norm));
	return norm;
}


void
normalize_vector(
                 double *v1,
                 unsigned int len
                )
{
	double norm = vector_norm(v1, len);
	multiplication_with_constant(v1, v1, 1./norm, len);
}

double
average(
        double *v1,
        unsigned int len
       )
{
	double sum = 0;
	for(int i = 0; i < len; ++i)
	{
		sum += v1[i];
	}
	return sum / (double) len;
}


double
standard_deviation(
                       double *v1,
                       unsigned int len
                  )
{
	double mean = average(v1, len);
	double *diff =  malloc(sizeof(double[len]));
	addition_with_constant(diff, v1, -mean, len);
	double sum = dot_product(diff, diff, len);
	
	free(diff);
	return sqrt(sum / len);
}

double
distance_between_vectors(
                         double *v1,
                         double *v2,
                         unsigned int len
                        )
{
	double * u = malloc(sizeof(double[len]));
	double dist;
	for(int i = 0; i<len; ++i)
	{
		u[i] = v1[i] - v2[i];
	}
	dist = vector_norm(u, len);
	free(u);
	return dist;
}

void
cumulative_integration(
                       double *res,
                       double *v,
                       double dx,
                       unsigned int v_len
                      )
{
	res[0] = 0;
	for(int i = 1; i<v_len; ++i)
	{
		res[i] = dx * (v[i] + v[i-1]) / 2. + res[i-1]; // cumulative??
	}
}


void
write_xyz(
          FILE *fp,
          char *symbol,
          double **positions,
          double **velocities,
          double alat,
          int natoms)
{
	fprintf(fp, "%i\nLattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\" ", natoms, alat, alat, alat);
    	fprintf(fp, "Properties=species:S:1:pos:R:3:vel:R:3 pbc=\"T T T\"\n");
	for(int i = 0; i<natoms; ++i)
	{
		fprintf(fp, "%s %lf %lf %lf ", symbol, positions[i][0], positions[i][1], positions[i][2]);
        	fprintf(fp, "%lf %lf %lf\n", velocities[i][0], velocities[i][1], velocities[i][2]);
	}
	
}

void fft_freq(
          double *res,
          int n,
          double timestep)
{
	int i;
	if(n % 2 == 0)
	{
		for(i = 0; i<n/2; ++i)
		{
			res[i] = (double) i / (timestep * (double) n);
		}
		
		for(; i<n; ++i)
		{
			res[i] = (double)(i-n) /  (timestep * (double) n);
		}
		//res[41] = 1.006291;
		//res[215] = -1.006291;
	}
	else
	{
		for(i = 0; i<=(n-1)/2; ++i)
		{
			res[i] = (double) i /  (timestep * (double) n);
		}
		
		for(; i<n; ++i)
		{
			res[i] = (double) (i-n) /  (timestep * (double) n);
		}
	}
	
	for(i = 0; i<n; ++i){
		res[i] = res[i] * 2 * M_PI;
	}

}

/* Freely given functions */
void
skip_line(FILE *fp)
{
    int c;
    while (c = fgetc(fp), c != '\n' && c != EOF);
}
/* Commented to be able to use Werror
void
read_xyz(
         FILE *fp,
         char *symbol,
         double **positions,
         double **velocities,
         double *alat)
{
    int natoms;
    if(fscanf(fp, "%i\nLattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 %lf\" ", &natoms, alat, alat, alat) == 0){
        perror("Error");
    }
    skip_line(fp);
    for(int i = 0; i < natoms; ++i){
        fscanf(fp, "%s %lf %lf %lf ",
                symbol, &positions[i][0], &positions[i][1], &positions[i][2]);
        fscanf(fp, "%lf %lf %lf\n",
                &velocities[i][0], &velocities[i][1], &velocities[i][2]);
    }
}

*/
void powerspectrum(
           double *res,
           double *signal,
           int n,
                   double timestep)
{
    /* Declaration of variables */
    int i;
    double *complex_coefficient = malloc(sizeof(double) * 2*n); // array for the complex fft data
    double *data_cp = malloc(sizeof(double) * n);

    /*make copy of data to avoid messing with data in the transform*/
    for (i = 0; i < n; i++) {
    data_cp[i] = signal[i];
    }

    /* Declare wavetable and workspace for fft */
    gsl_fft_real_wavetable *real;
    gsl_fft_real_workspace *work;

    /* Allocate space for wavetable and workspace for fft */
    work = gsl_fft_real_workspace_alloc(n);
    real = gsl_fft_real_wavetable_alloc(n);

    /* Do the fft*/
    gsl_fft_real_transform(data_cp, 1, n, real, work);

    /* Unpack the output into array with alternating real and imaginary part */
    gsl_fft_halfcomplex_unpack(data_cp, complex_coefficient,1,n);

    /*fill the output powspec_data with the powerspectrum */
    for (i = 0; i < n; i++) {
    res[i] = (complex_coefficient[2*i]*complex_coefficient[2*i]+complex_coefficient[2*i+1]*complex_coefficient[2*i+1]);
    res[i] *= timestep / n;
    }

    /* Free memory of wavetable and workspace */
    gsl_fft_real_wavetable_free(real);
    gsl_fft_real_workspace_free(work);
    free(complex_coefficient);
    free(data_cp);
}
