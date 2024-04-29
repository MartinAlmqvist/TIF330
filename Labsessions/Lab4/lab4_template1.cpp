#include <iostream>
#include <vector>
#include "plot.h"
#include <cmath>

using namespace std;

double epsilon = 1e-10;

void thomasAlgorithm(double *a, double *b, double *c, double *d, double *x, int n)
{
    if(n == 1){
        x[0] = d[0]/b[0];
    } else {
        c[0] = c[0]/b[0];
        d[0] = d[0]/b[0];
        for(int i = 1; i < n; i++){
            c[i] = c[i]/(b[i] - a[i]*c[i-1]);
            d[i] = (d[i] - a[i]*d[i-1])/(b[i] - a[i]*c[i-1]);
        }
        x[n-1] = d[n-1];
        for(int i = n-2; i >= 0; i--)
            x[i] = d[i] - c[i]*x[i+1];
    }
}

bool drum1(double x, double y) // return true for inner region
{
    if((y > 0)&&(y < 1/3.0)&&(x > 2/3.0 + epsilon)&&(y > x - 2/3.0 + epsilon)) return true;
    if((y >= 1/3.0)&&(y < 2/3.0)&&(y > -x + 2/3.0 + epsilon)&&(x < 1 - epsilon)) return true;
    if((y >= 2/3.0)&&(y < 1)&&(y < x + 2/3.0 - epsilon)&&(x < 1/3.0 - epsilon)) return true;
    return false;
}
bool drum2(double x, double y) // return true for inner region
{

    if((y > 0)&&(y < 1/3.0)&&(y > -x + 2/3.0 + epsilon)&&(x < 2/3.0 - epsilon)) return true;
    if((y >= 1/3.0)&&(y < 2/3.0)&&(y > -x + 2/3.0 + epsilon)&&(y < -x + 4/3.0 - epsilon)) return true;
    if((y >= 2/3.0)&&(y < 1)&&(x > 0)&&(x < 1/3.0 - epsilon)) return true;
    return false;
}

bool drum3(double x, double y) // return true for inner region
{
    if(((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5) < 0.24)&&((y < 1/2.0 - epsilon)||(x < 1/3.0 - epsilon))) return true;
    return false;
}

void computeNextState(double ***states, double *lambda, int &state_num, int n, double T = 1.0)
{
    double h = 1/double(n - 1);
    double tau = 0.02*h;
    double sigma = tau/(2*h*h);

    double *a = new double[n];
    double *b = new double[n];
    double *c = new double[n];
    double *d = new double[n];
    double *x = new double[n];

    double **D = new double*[n];
    for(int i = 0; i < n; i++) D[i] = new double[n];

    for(int ix = 0; ix < n; ix++)
    for(int iy = 0; iy < n; iy++){
        D[iy][ix] = drum1(ix/double(n - 1), iy/double(n - 1));
    }
    pyPlot2D(D, n, n, "drum", ", cmap='RdBu'", true);
    
    double **Ua = new double*[n];
    double **Ub = new double*[n];
    for(int iy = 0; iy < n; iy++) {
        Ua[iy] = new double[n];
        Ub[iy] = new double[n];
        for(int ix = 0; ix < n; ix++) {
            Ua[iy][ix] = -sin(M_PI*(0.5 + (ix + 0.5)/double(n)))*sin(M_PI*(ix + 0.5)/double(n))*sin(M_PI*(iy + 0.5)/double(n))*D[iy][ix];
            Ub[iy][ix] = 0;
        }
    }
    //plot2D_dots(Ua, n, n, "U0", ", cmap='RdBu'");
    int stride = 100;
    for(int it = 0; it < T/tau; it++){
        //cout << it << endl;
        for(int iy = 0; iy < n; iy++){
            int is = 0; // number of grid points within the drum for a given iy
            int i0 = n; // the first grid point inside the drum (n is the default value that means that ix is beyond the membrain)
            for(int ix = 0; ix < n; ix++){
                if((D[iy][ix] != 0)&&(i0 != n)){// internal or final point
                    b[is] = 1+2*sigma;
                    c[is] = -sigma;
                    a[is] = -sigma;
                    // if(ix = 0 && ix + 1 >= n){
                    //     d[is] = 0;
                    // }
                    d[is] = sigma*(Ua[iy-1][ix] -2*Ua[iy][ix] + Ua[iy+1][ix])/2 + Ua[iy][ix];
                    is++;
                }
                if((D[iy][ix] != 0)&&(i0 == n)){// starting point
                    i0 = ix; 
                    b[0] = 1+2*sigma;
                    c[0] = -sigma;
                    // if(ix = 0 && ix + 1 >= n){
                    //     d[is] = 0;
                    // }
                    d[0] =  sigma*(Ua[iy-1][ix] -2*Ua[iy][ix] + Ua[iy+1][ix])/2 + Ua[iy][ix];
                    is = 1;
                }
                if((is != 0)&&(i0 != 0)&&(D[iy][ix + 1] == 0)) ix = n; // termination
            }
            thomasAlgorithm(a, b, c, d, x, is);
            for(int i = 0; i < is; i++) Ub[iy][i0 + i] = x[i];
        }
        for(int ix = 0; ix < n; ix++){
            int is = 0; // number of grid points within the drum for a given iy
            int i0 = n; // the first grid point inside the drum (n is the default value that means that ix is beyond the membrain)
            for(int iy = 0; iy < n; iy++){
                if((D[iy][ix] != 0)&&(i0 != n)){// internal or final point
                    b[is] = 1+2*sigma;
                    c[is] = -sigma;
                    a[is] = -sigma;
                    // if(ix = 0 && ix + 1 >= n){
                    //     d[is] = 0;
                    // }
                    d[is] = sigma*(Ub[iy][ix-1] -2*Ub[iy][ix] + Ub[iy][ix+1])/2 + Ub[iy][ix];
                    is++;
                }
                if((D[iy][ix] != 0)&&(i0 == n)){// starting point
                    i0 = iy; 
                    b[0] = 1+2*sigma;
                    c[0] = -sigma;
                    // if(ix = 0 && ix + 1 >= n){
                    //     d[is] = 0;
                    // }
                    d[0] = sigma*(Ub[iy][ix-1] -2*Ub[iy][ix] + Ub[iy][ix+1])/2 + Ub[iy][ix];
                    is = 1;
                }
                if((is != 0)&&(i0 != 0)&&(D[iy + 1][ix] == 0)) iy = n; // termination
            }
            thomasAlgorithm(a, b, c, d, x, is);
            for(int i = 0; i < is; i++) Ua[i0 + i][ix] = x[i];
        }
        if(it%stride == 0)
        {
            for(int i = 0; i < state_num; i++){
                double C = 0;
                for(int iy = 0; iy < n; iy++)
                for(int ix = 0; ix < n; ix++)
                    if(D[iy][ix] != 0) C += Ua[iy][ix] * states[i][iy][ix] * (h*h);    
                for(int iy = 0; iy < n; iy++)
                for(int ix = 0; ix < n; ix++)
                    if(D[iy][ix] != 0) Ua[iy][ix] -= C * states[i][iy][ix];    
            }
            double S = 0;
            for(int iy = 0; iy < n; iy++)
                for(int ix = 0; ix < n; ix++)
                    if(D[iy][ix] != 0) S += Ua[iy][ix] * Ua[iy][ix] * (h*h);
            double factor = 1/sqrt(S);
            double l = -log(sqrt(S))/(stride*tau);
            lambda[state_num] = l;
            cout << it << "/" << T/tau << ", l = " << l/9.0 << endl;
            for(int iy = 0; iy < n; iy++)
                for(int ix = 0; ix < n; ix++){
                    if(D[iy][ix] != 0) Ua[iy][ix] *= factor;
                    if(D[iy][ix] != 0) states[state_num][iy][ix] = Ua[iy][ix];
                }
            //plot2D_dots(Ua, n, n, "U" + to_string(state_num) + "_" + to_string(it/stride), ", cmap='RdBu'");
        }
    }
    state_num++;
}

double max2(double **state, int n)
{
    double max = 0;
    for(int iy = 0; iy < n; iy++) 
        for(int ix = 0; ix < n; ix++) 
            if(state[iy][ix]*state[iy][ix] > max) max = state[iy][ix]*state[iy][ix];
    return max;
}

int main()
{
    int n = 3*40 + 1;
    int nS = 10;
    double *lambda = new double[nS];
    double ***states = new double**[nS];
    for(int i = 0; i < nS; i++){
        states[i] = new double*[n];
        for(int iy = 0; iy < n; iy++) {
            states[i][iy] = new double[n];
            for(int ix = 0; ix < n; ix++) states[i][iy][ix] = -1;
        }
    }
    int state_num = 0;
    for(int i = state_num; i < nS; i++) computeNextState(states, lambda, state_num, n, 0.25);
    for(int i = 0; i < state_num; i++)
    {
        cout << "lambda_" << i << " = " << lambda[i]/lambda[0]  << ", max2 = " << max2(states[i], n) << endl;
        pyPlot2D(states[i], n, n, "T" + to_string(i), ", cmap='RdBu'", true);
    }

    for(int i = 0; i < nS; i++){
        for(int iy = 0; iy < n; iy++) delete []states[i][iy];
        delete []states[i];
    }
    delete states;
    delete []lambda;
}