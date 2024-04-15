//g++ -o Task3 A1_Task3.cpp

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>

// Constants
const double dx = 0.001; // Spatial step size
const double dt = 0.001; // Time step size
const long nx = 2/dx + 1; // Number of space points
const long nt = 5/dt + 1; // Number of time steps

bool inf_space = false; // boolean for determining BCs

// Set initial conditions
void initializeFields(double Ey[], double Bz[]) {
    for (int i = 0; i < nx; ++i) {
        // Range x from -1 to 1 
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

// Perform the FDTD of E and B
void updateFDTD(double Ey[], double Bz[]) {
    // Update Ey
    for (int k = 1; k < nx ; ++k) {
        // Update E
        Ey[k] = Ey[k] + dt*(Bz[k-1] - Bz[k])/dx;
    }

    // Update B
    for (int k = 0; k < nx ; ++k) {
        Bz[k] = Bz[k] + dt*(Ey[k] - Ey[k+1]) / dx;
    }
}

// Update boundary conditions (of ideal conductor or infinite space)
void BoundaryConditions(double Ey[], double Bz[]) {
    // infinite space
    if (inf_space == true){
        Ey[0] = Ey[1];
        Ey[nx] = Ey[nx-1];
        Bz[0] = Bz[1];
        Bz[nx] = Bz[nx-1];
    }
    // ideal condutctor
    else{
        Ey[0] = 0.0;
        Ey[nx] = 0.0;
        Bz[0] = 0.0;
        Bz[nx] = 0.0;
    }
 }

int main() {
    // Files for storing data
    std::ofstream outfile_Ey("results_Ey.txt");
    std::ofstream outfile_Bz("results_Bz.txt");

    // Arrays to store E and B
    double Ey[nx], Bz[nx];
    // Initialize fields  at t = 0
    initializeFields(Ey, Bz);

    // Iterate over time
    for (int t = 0; t < nt; ++t) {

        // Save E and B data
        for (int i = 0; i < nx; ++i) {
            outfile_Ey << Ey[i] << " ";
            outfile_Bz << Bz[i] << " ";
        }

        outfile_Ey << std::endl;
        outfile_Bz << std::endl;

        // Update the fields
        updateFDTD(Ey, Bz);
        // Update boundary conditions
        BoundaryConditions(Ey, Bz);
    }

    outfile_Ey.close();
    outfile_Bz.close();
    return 0;
}
