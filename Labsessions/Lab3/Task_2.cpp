#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

double* get_y(double base[], double weights, int N) {
    double* y_res = new double[N];

    for (int i = 0; i < N; i++) {
        y_res[i] = weights * base[i];
    }
    return y_res;
}

double* get_y_dx(double k, int N) {
    double* dydx = new double[N];

    for (int i = 0; i < N; i++) {
        double x = -1 + i * (2.0 / N);
        dydx[i] = -M_PI / 2 * ((1 + 2 * k) * x) * sin(M_PI / 2 * ((1 + 2 * k) * x));
    }
    return dydx;
}

double* get_f(double dydx[], double a, int N) {
    double* f = new double[N];

    for (int i = 0; i < N; i++) {
        f[i] = sqrt(1 + a * dydx[i] * dydx[i]);
    }

    return f;
}

double* integrate(double base[], double f[], int N) {
    double* f_hat = new double[N];
    double dx = 2.0 / N;
    double sum = 0.0;

    for (int i = 0; i < N - 1; i++) {
        sum += 0.5 * dx * (base[i] * f[i] + base[i + 1] * f[i + 1]);
        f_hat[i] = sum;
    }

    f_hat[N - 1] = f_hat[N - 2];

    return f_hat;
}

double update_weights(double base[], int N, double f_hat[]) {
    double weights = 0.0;

    for (int i = 0; i < N; i++) {
        if (base[i] == 0) {
            weights = 0;// if divide with zero
        } else {
            weights += 1 / (base[i] * base[i]) * f_hat[i];
        }
    }

    return weights;
}

int main() {
    int N = 100;

    double* base = new double[N];
    double weights = -1.0;
    double* y  = new double[N];
    double* dydx  = new double[N];
    double* f  = new double[N];
    double* f_hat  = new double[N];
    double a = 1;


    std::ofstream file("output.txt");

    for (int k = 1; k < 100; k++) {
        // Basis function
        for (int i = 0; i < N; i++) {
           
            double x = -1 + i * (2.0 / N);
            base[i] = cos(M_PI / 2 * ((1 + 2 * k) * x));
        }

        y = get_y(base, weights, N);

        dydx = get_y_dx(k, N);
        f = get_f(dydx, a, N);
        f_hat = integrate(base, f, N);
        weights = update_weights(base, N, f_hat);
          
        delete[] dydx;
        delete[] f;
        delete[] f_hat;
    }

    y = get_y(base, weights, N);

    // Save values of y to CSV file
    for (int i = 0; i < N; i++) {
         file << y[i] << std::endl;
    }

    
    delete[] base;
    delete[] y;

    return 0;
}
