#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>

// Constants

const double dx = 0.02; // Spatial step size
const double dt = 0.01; // Time step size
const int nx = 111;     // Number of space points
const int nt = 501;     // Number of time steps

// Set initial conditions
void initializeFields(double Ey[], double Bz[]) {
    for (int i = 0; i < nx; ++i) {
        // Range x from -0.1 to 0.1 
        double x = -0.1 + i * dx;
        if (x < 0.1) {
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
    for (int i = 1; i < nx - 1; ++i) {
        Ey[i] = Ey[i] - dt * (Bz[i + 1] - Bz[i]) / dx;
    }

    // Update Bz
    for (int i = 1; i < nx - 1; ++i) {
        Bz[i] = Bz[i] - dt * (Ey[i] - Ey[i - 1]) / dx;
    }
}

int main() {
    // File for storing data
    std::ofstream outfile("results.txt");

    // Arrays to store electric and magnetic fields
    double Ey[nx], Bz[nx];

    // Initialize fields
    initializeFields(Ey, Bz);

    // Main loop for time evolution
    for (int t = 0; t < nt; ++t) {

        // Write current values of Ey and Bz to file
        for (int i = 0; i < nx; ++i) {
            outfile << Ey[i] << " " << Bz[i] << std::endl;
        }

        // Update fields using FDTD method
        updateFDTD(Ey, Bz);
    }

    outfile.close();
    return 0;
}
