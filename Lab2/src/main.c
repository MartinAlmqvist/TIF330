#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "tools.h"



int
run(
    int argc,
    char *argv[]
   )
{
    // use finitie integral method


    //grid - allocate memorey should create a poitner that points at the "layers"
    int firstLayer = 10; // define the number of points at bottom
    int x = 1;
    int l = x;
    int numberOfLayers = 10;

    double** vectorHolder = create_1D_array(firstLayer);
    for (int idx = 0; idx<numberOfLayers; idx+2){
        vectorHolder[idx] = create_1D_array(numberOfLayers-idx);
        vectorHolder[idx+1] = create_1D_array(numberOfLayers-idx);
    }

    


    //save value log of ratio is freq * diffrence in time.




    // int N = 1000; //Change this 

    // //Calulate the best P for diffrent t. 
    // int iterations = 1000;
    // double p_optimal_low=-1.;
    // double p_optimal_high=0.;
    // double pdt = 1000000.0; //step size will be 1/pdt.
    // double Energy_low = 0.0;
    // double Energy_high = 0.0;
  


    //  FILE *H2a= fopen("H2a.csv", "w+");
    // double step_size = 1.0 / pdt; // Calculate the step size

    // for (int t = 1; t < iterations; t++) // t in kelvin.
    // {
    //     double temperture = t ;
    //     double initial_low = meanfieldFreeEnergy(N, -1.0+step_size, temperture);
    //     double initial_high = meanfieldFreeEnergy(N, 0., temperture);
    //     for (int p = 0; p <= pdt; p++) // Runs p from 0-1 with 1/pdt steps
    //     {
    //         double p_step_low = -1.0+step_size + p * step_size; 
    //         double p_step_high = p * step_size; 
    //         double temp_low = meanfieldFreeEnergy(N, p_step_low, temperture);
    //         double temp_high = meanfieldFreeEnergy(N, p_step_high, temperture);
    //         if (temp_low < initial_low)
    //         { // wants to test different p to minimize the energy
    //             p_optimal_low = p_step_low;
    //             initial_low = temp_low;
    //             Energy_low = meanfieldEnergy(N, p_optimal_low); 
    //         }
    //         if (temp_high < initial_high)
    //         { // wants to test different p to minimize the energy
    //             p_optimal_high = p_step_high;
    //             initial_high = temp_high;
    //             Energy_high = meanfieldEnergy(N, p_optimal_high); 
    //         }
    //     }
    //     fprintf(H2a, "%e \t %e \t %e \t %e\n", p_optimal_low, Energy_low, p_optimal_high,Energy_high);
    //     p_optimal_low=-1.;
    //     p_optimal_high=0.;
    // }
    // fclose(H2a); 

    // for(int t = 0; t < iterations; t++) //t in kelvin.
    // {
    //     double initial = meanfieldFreeEnergy(N,-1.,t);
    //     for(int p =-100000; p < pdt; p++) //Runs p from 0-1 with 1/pdt steps
    //     {
    //        double p_step = (double)p/pdt;
    //        double temp = meanfieldFreeEnergy(N,p_step,t);
    //        if(temp < initial)
    //        {// wants to test diffrent p to minimize the energy
    //             p_optimal = p_step;
    //             initial = temp;
    //             temperature = meanfieldEnergy(N, p_optimal); //wrong here, should be energy
               
    //        } 
           
    //     }
    //     fprintf(H2a,"%f \t %f\n" ,p_optimal,temperature);
    //     p_optimal = 0.0;
    // }



    return 0;
}
