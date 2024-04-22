#include <iostream>
#include <cmath>
#include <complex>
#include <vector>



// Potential function
double potential(double x, double y) {
    return -5 * pow(1 + pow(x / 5.0, 2) + pow(y / 4.0, 2), -4);
}

// Wave function at t=0
double initialWave(double x, double y) {
    return exp(-(pow(x - 1, 2) + pow(y - 1, 2)) / M_PI);
}

// Function to construct the kinetic energy matrix A
double kineticEnergy(double kx, double ky) {
    return kx * kx + ky * ky;
}

// Split-step method for time evolution
double splitStep(double A, double B, int t) {
    
    // Calculate the split-step propagator P_AB using the matrices A and B
    double Pab = exp(t/2*(A(t)))*exp(t/2*(B(t)))*exp(t/2*(B(0)))*exp(t/2*(A(0)));

    return Pab;
}

// Function to perform the Fourier transform (FFT)
double performFFT(double x) {
    return x;
}

int main() {
    // Example usage
    double kx = 1.0; 
    double ky = 1.0;
    const double xy_low = -10;
    const double xy_high = 10;
    const int t_start = 0;
    const int t_end = 100;
    double dx = 1;
    double dt = 1;
    double size = 20/dx;
        //Allocate A memorey
    double** A = new double*[size];
    for (int i = 0; i < size; ++i) {
        A[i] = new double[size];
    }
    //Allocate B memorey
    double** B = new double*[size];
    for (int i = 0; i < size; ++i) {
        B[i] = new double[size];
    }
    for (int i = 0; i < size; ++i) {
        double x = xy_low + i * dx;
        for (int j = 0; j < size; ++j) {
            double y = xy_low + j * dx;
            
            // Calculate potential energy at (x, y)
            B[i][j] = potential(x, y);
            
            // Calculate kinetic energy at (kx, ky)
            double kx = 1; //Need to find this
            double ky = 1;
            A[i][j] = kineticEnergy(kx, ky);
        }
    }
  
    double psi = initialWave(x,y);
    for (double t = t_start; t <= t_end; t += dt) {
        //Uppdate psi according to Split-step
        psi = splitStep(A, B, dt) * psi;
}


    return 0;
}
