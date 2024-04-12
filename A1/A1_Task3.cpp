//g++ -o Task3 A1_Task3.cpp

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>

// Constants
const double c = 3*pow(10,8);
const double dx = 0.2; // Spatial step size
const double dt = dx/(2*c); // Time step size
const int nx = 2/dx + 1;     // Number of space points
const int nt = 5/dt + 1;     // Number of time steps


// Set initial conditions
void initializeFields(double Ey[], double Bz[]) {
    for (int i = 0; i < nx; ++i) {
        // Range x from -0.1 to 0.1 
        double x = -1 + i * dx;

        if (x > -0.1 && x<0.1) {
            Ey[i] = sin(20 * M_PI * x);
            Bz[i] = sin(20 * M_PI * x);
        } else {
            Ey[i] = 0.0;
            Bz[i] = 0.0;
        }
    }
}

// Function to perform FDTD update
void updateFDTD(double Ey[], double Bz[]) {
    // Update Ey
    for (int k = 1; k < nx ; ++k) {
        // use previous value of E and B to define E at new time
        Ey[k] = Ey[k] - dt * (Bz[k] - Bz[k-1]) / dx;
    }

    // Update Bz
    for (int k = 1; k < nx ; ++k) {
        Bz[k] = Bz[k] - dt * (Ey[k] - Ey[k-1]) / dx;
    }
}

int main() {
    // File for storing data
    std::ofstream outfile_Ey("results_Ey.txt");
    std::ofstream outfile_Bz("results_Bz.txt");

    // Arrays to store electric and magnetic fields
    double Ey[nx], Bz[nx];

    // Initialize fields
    initializeFields(Ey, Bz);

    // Main loop for time evolution
    for (int t = 0; t < nt; ++t) {

        for (int i = 0; i < nx; ++i) {
            outfile_Ey << Ey[i] << " ";
        }
        outfile_Ey << std::endl;

        // Write current values of Bz to file
        for (int i = 0; i < nx; ++i) {
            outfile_Bz << Bz[i] << " ";
        }
        outfile_Bz << std::endl;

        // Update fields using FDTD method
        updateFDTD(Ey, Bz);
    }

    outfile_Ey.close();
    outfile_Bz.close();
    return 0;
}
