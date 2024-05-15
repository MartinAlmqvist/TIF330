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
    //return x*x;
}

// Function to initialize the wave function
void initialize_wavefunction(fftw_complex* psi, int N, double dx) {
    double x, y;
    for (int i = 0; i < N; ++i) {
        x = -10.0 + i * dx; // Ensure x stays within [-10, 10] given in the pdf
        for (int j = 0; j < N; ++j) {
            y = -10.0 + j * dx;
            psi[i * N + j][0] = 1.0 / M_PI * exp(-(x - 1) * (x - 1) - (y - 1) * (y - 1));
            psi[i * N + j][1] = 0.0;
            // psi[i * N + j][0] = 1.0 / M_PI * exp(-(x - 1) * (x - 1));
            // psi[i * N + j][1] = 0;
        }
    }
}

// Function to apply the potential operator in real space
void apply_potential_operator(fftw_complex* psi, int N, double dx, double dt) {
    double x, y;
    for (int i = 0; i < N; ++i) {
        x = -10.0 + i * dx;
        for (int j = 0; j < N; ++j) {
            y = -10.0 + j * dx;
            double pot = potential(x, y);
            double exp_real = cos(dt / 2 * pot);
            double exp_imag = -sin(dt / 2 * pot);
            double new_real = psi[i * N + j][0] * exp_real - psi[i * N + j][1] * exp_imag;
            double new_imag = psi[i * N + j][0] * exp_imag + psi[i * N + j][1] * exp_real;  
            
            psi[i * N + j][0] = new_real;
            psi[i * N + j][1] = new_imag;
        }
    }
}

// Function to perform forward Fourier transform to momentum space
void transform_to_momentum_space(fftw_complex* psi, int N) {
    fftw_plan p = fftw_plan_dft_2d(N, N, psi, psi, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
}

// Function to apply the kinetic energy operator in momentum space
void apply_kinetic_operator(fftw_complex* psi, int N, double dx, double dt) {
    double kx, ky;
    double kx_factor = 2.0 * M_PI / (N * dx);
    double ky_factor = 2.0 * M_PI / (N * dx);
    for (int i = 0; i < N; ++i) {
        kx = (i < N / 2) ? i * kx_factor : (i - N) * kx_factor;
        for (int j = 0; j < N; ++j) {
            ky = (j < N / 2) ? j * ky_factor : (j - N) * ky_factor;
            double k_squared = kx * kx + ky * ky;
            double exp_real = cos(dt * k_squared);
            double exp_imag = -sin(dt * k_squared);
            double new_real = psi[i * N + j][0] * exp_real - psi[i * N + j][1] * exp_imag;
            double new_imag = psi[i * N + j][0] * exp_imag + psi[i * N + j][1] * exp_real;
            psi[i * N + j][0] = new_real;
            psi[i * N + j][1] = new_imag;
        }
    }
}

// Function to transform the wave function back to position space
void transform_to_position_space(fftw_complex* psi, int N) {
    fftw_plan p = fftw_plan_dft_2d(N, N, psi, psi, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    double normalization_factor = 1.0 / (N * N);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            psi[i * N + j][0] *= normalization_factor;
            psi[i * N + j][1] *= normalization_factor;
        }
    }
}

int main() {
    const int N = 3*512; // Number of grid points
    const double L = 20.0; // Size of the spatial grid
    const double dx = L / N; // Spatial step
    const int tmax = 5000;

    fftw_complex* psi = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);

    // Initialize the wave function
    initialize_wavefunction(psi, N, dx);

    // Open file for writing
    std::ofstream outfile("wave_function_values.csv");
    if (!outfile.is_open()) {
        std::cerr << "Error opening file for writing." << std::endl;
        return 1;
    }
    std::cout << " Start the loop version 5 " << std::endl;
    for (int t = 0; t < tmax; ++t) {
        double dt = 100.0 / tmax;
        if( t%1000 == 0) {
            std::cout << t << std::endl;
        }
        // Apply the potential operator in real space
        apply_potential_operator(psi, N, dx, dt);

        // Transform the wave function to momentum space
        transform_to_momentum_space(psi, N);

        // Apply the kinetic energy operator in momentum space
        apply_kinetic_operator(psi, N, dx, dt);

        // Transform the wave function back to position space
        transform_to_position_space(psi, N);

        // Apply the potential operator in real space
        apply_potential_operator(psi, N, dx, dt);

        int x_index = static_cast<int>((0.1+10.0) / 20.0 * (N-1)); // Convert from range [-10, 10] to [0, N-1]
        int y_index = N / 2; // Y index for y = 0
        int index_01 = y_index * N + x_index;
        outfile << t * dt << "," << psi[x_index * N + y_index][0] << "," << psi[x_index * N + y_index][1] << std::endl;

        //outfile << t*dt << "," << psi[N / 2 * N / 2 + statisssssssc_cast<int>(0.1 * N/2)][0] << "," << psi[N / 2 * N/2][1] << std::endl;
    }
    outfile.close();
 
    
    fftw_free(psi);

    std::cout << "Finished" << std::endl;
    return 0;
}
