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
    // return x*x +y*y;
}

// Function to initialize the wave function
void initialize_wavefunction(std::complex<double>** psi, int N, double dx) {
    double x, y;
    for (int i = 0; i < N; ++i) {
        double x = double(i)/double(N)*20-10; // Ensure x stays within [-10, 10] given in the pdf
        for (int j = 0; j < N; ++j) {
            double y = double(j)/double(N)*20-10;
            psi[i][j] = 1/M_PI * std::exp(-(x - 1) * (x - 1) - (y - 1) * (y - 1));
        }
    }
}

// Function to apply the potential operator in real space
void apply_potential_operator(std::complex<double>** psi, int N, double dx, double dt) {
    for (int i = 0; i < N; ++i) {
        double x = double(i)/double(N)*20-10;
        for (int j = 0; j < N; ++j) {
            double y = double(j)/double(N)*20-10;
            psi[i][j] *=  std::exp(-std::complex<double>(0, dt/2 * potential(x, y)));
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
// void apply_kinetic_operator(std::complex<double>** psi, int N, double dk,double dt) {
//     double kx, ky;
//     for (int i = 0; i < N; ++i) {
//         if (i <= N / 2) {
//             kx = i * dk;
//         } else {
//             kx = (i - N) * dk;
//         }
//         for (int j = 0; j < N; ++j) {
//             if (j <= N / 2) {
//                 ky = j * dk;
//             } else {
//                 ky = (j - N) * dk;
//             }
//             double k2 = kx * kx + ky * ky;
//             psi[i][j] *= std::exp(-std::complex<double>(0, dt* k2));
//         }
//     }
// }

void apply_kinetic_operator(std::complex<double> **psi, int N, double dx, double dt) {
    double kx, ky;
    double kx_factor = 2.0 * M_PI / (N * dx);
    double ky_factor = 2.0 * M_PI / (N * dx);
    
    for (int i = 0; i < N; ++i) {
        kx = i < N / 2 ? i * kx_factor : (i - N) * kx_factor;
        for (int j = 0; j < N; ++j) {
            ky = j < N / 2 ? j * ky_factor : (j - N) * ky_factor;
            double k_squared = kx*kx + ky*ky;
            psi[i][j] *= std::cos(dt * k_squared) - std::sin(dt * k_squared) * std::complex<double>(0, 1);
        }
    }
}

void transform_to_position_space(std::complex<double>** psi, int N) {
    // Cast the pointer psi[0] to fftw_complex* for input and output arrays
    fftw_complex *in = reinterpret_cast<fftw_complex*>(psi[0]);
    fftw_complex *out = reinterpret_cast<fftw_complex*>(psi[0]);

    // Create the FFTW plan for 2D backward FFT
    fftw_plan plan = fftw_plan_dft_2d(N, N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Execute the FFTW plan
    fftw_execute(plan);

    // Destroy the FFTW plan to free resources
    fftw_destroy_plan(plan);

    // Calculate the normalization factor
    double normalization_factor = 1.0 / (N * N);

    // Normalize the output array
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // Adjust the real and imaginary parts separately
            psi[i][j].real(psi[i][j].real() * normalization_factor);
            psi[i][j].imag(psi[i][j].imag() * normalization_factor);
        }
    }
}

int main() {
    const int N = 100; // Number of grid points
    const double L = 20.0; // Size of the spatial grid
    const double dx = L / N; // Spatial step
    const int tmax = 1000;

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
        double dt = 100.0 / tmax;

        // Apply the potential operator in real space
        apply_potential_operator(psi, N, dx,dt);

        // Transform the wave function to momentum space
        transform_to_momentum_space(psi, N);

        // Apply the kinetic energy operator in momentum space
        apply_kinetic_operator(psi, N, dx,dt);

        // Transform the wave function back to position space
        transform_to_position_space(psi, N);

        // Apply the potential operator in real space
        apply_potential_operator(psi, N, dx,dt);
        
        // Save the time and wave function value at point (0.1, 0) to CSV
        outfile << t << "," << std::real(psi[N/2][N/2]) << "," << std::imag(psi[N/2][N/2]) << std::endl;
    }
    outfile.close();


    std::cout << "Hello there" << std::endl;    
    return 0;

}
