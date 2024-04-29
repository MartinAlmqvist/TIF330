#include <iostream>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <fstream>

// Define the potential function
double potential(double x, double y) {
    double term2 = pow((x / 5.0), 2);
    double term3 = pow((y / 4.0), 2);
    return -5.0 * pow(1 + term2 + term3, -4);
}

// Function to initialize the wave function
void initialize_wavefunction(std::complex<double>** psi, int N, double dx) {
    double x, y;
    for (int i = 0; i < N; ++i) {
        x = -10.0 + i * dx; // Ensure x stays within [-10, 10] given in the pdf
        for (int j = 0; j < N; ++j) {
            y = -10.0 + j * dx;
            psi[i][j] = std::exp(-(x - 1) * (x - 1) - (y - 1) * (y - 1));
        }
    }
}

// Function to apply the potential operator in real space
void apply_potential_operator(std::complex<double>** psi, int N, double dx) {
    double x, y;
    for (int i = 0; i < N; ++i) {
        x = -10.0 + i * dx;
        for (int j = 0; j < N; ++j) {
            y = -10.0 + j * dx;
            psi[i][j] *= std::exp(-std::complex<double>(0, potential(x, y)));
        }
    }
}

// Function to perform forward Fourier transform to momentum space
void transform_to_momentum_space(std::complex<double> **psi, int N) {
    fftw_complex *in;
    fftw_complex *out;
    fftw_plan p;
    
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    
    // Copy data from psi to in array
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            in[i * N + j][0] = psi[i][j].real();  // Real part
            in[i * N + j][1] = psi[i][j].imag();  // Imaginary part
        }
    }
    
    p = fftw_plan_dft_2d(N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    
    // Copy the transformed data back to psi
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            psi[i][j] = std::complex<double>(out[i * N + j][0], out[i * N + j][1]);
        }
    }
    
    fftw_free(in);
    fftw_free(out);
}

// Function to apply the kinetic energy operator in momentum space
void apply_kinetic_operator(std::complex<double>** psi, int N, double dk) {
    double kx, ky;
    for (int i = 0; i < N; ++i) {
        if (i <= N / 2) {
            kx = i * dk;
        } else {
            kx = (i - N) * dk;
        }
        for (int j = 0; j < N; ++j) {
            if (j <= N / 2) {
                ky = j * dk;
            } else {
                ky = (j - N) * dk;
            }
            double k2 = kx * kx + ky * ky;
            psi[i][j] *= std::exp(-std::complex<double>(0, k2 / 2));
        }
    }
}

// Function to perform inverse Fourier transform to position space
void transform_to_position_space(std::complex<double>** psi, int N) {
    fftw_complex *in;
    fftw_complex *out;
    fftw_plan p;
    
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    
    // Copy data from psi to in array
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            in[i * N + j][0] = psi[i][j].real();  // Real part
            in[i * N + j][1] = psi[i][j].imag();  // Imaginary part
        }
    }
    
    p = fftw_plan_dft_2d(N, N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    
    // Copy the transformed data back to psi
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            psi[i][j] = std::complex<double>(out[i * N + j][0] / (N * N), out[i * N + j][1] / (N * N));
        }
    }
    
    fftw_free(in);
    fftw_free(out);
}

int main() {
    const int N = 64; // Number of grid points
    const double L = 20.0; // Size of the spatial grid
    const double dx = L / N; // Spatial step
    const double dk = 2 * M_PI / L; // Fourier step
    const int tmax = 100;

    // Allocate memory for the wave function
    std::complex<double>** psi = new std::complex<double>*[N];
    for (int i = 0; i < N; ++i) {
        psi[i] = new std::complex<double>[N];
    }

    // Initialize the wave function
    initialize_wavefunction(psi, N, dx);
    
    // Open file for writing
    std::ofstream outfile("wave_function_values.csv");
    if (!outfile.is_open()) {
        std::cerr << "Error opening file for writing." << std::endl;
        return 1;
    }
    
    for (int t = 0; t < tmax; ++t) {
        // Apply the potential operator in real space
        apply_potential_operator(psi, N, dx);

        // Transform the wave function to momentum space
        transform_to_momentum_space(psi, N);

        // Apply the kinetic energy operator in momentum space
        apply_kinetic_operator(psi, N, dk);

        // Transform the wave function back to position space
        transform_to_position_space(psi, N);
        
        // Save the time and wave function value at point (0.1, 0) to CSV
        outfile << t << "," << std::real(psi[N/2][N/2]) << "," << std::imag(psi[N/2][N/2]) << std::endl;
    }
    outfile.close();

    // Free memory
    for (int i = 0; i < N; ++i) {
        delete[] psi[i];
    }
    delete[] psi;

    return 0;
}
