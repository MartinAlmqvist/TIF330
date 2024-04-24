#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

struct catenary // structure for solving the catenary problem (y'' = sqrt(1 + a*y'^2) with y(-1) = y(1) = 0) using iterative application of Galerkin method
{
    vector<double> c; // coefficients of representation
    double a; 
    catenary(int N, double a): a(a) // constructor function (function that creates an instance of the class)
    {
        c.assign(N, 0);
        c[0] = -1; // initial state
    }
    double sqr(double x) // auxilary function
    {
        return x*x;
    }
    double basisFunction(int k, double x) // k-th basis function
    {
        return cos(0.5*M_PI*(1 + 2*k)*x); // this form gives y(-1) = y(1) = 0 and also account for the symmetry of the problem
    }
    double basisFunction_x(int k, double x) // the first derivative of the k-th basis function 
    {
        return -0.5*M_PI*(1 + 2*k)*sin(0.5*M_PI*(1 + 2*k)*x);
    }
    double solution(double x) // current solution
    {
        double S = 0;
        for(int k = 0; k < c.size(); k++)
            S += c[k]*basisFunction(k, x);
        return S;
    }
    double solution_x(double x) // derivative of the current solution
    {
        double S = 0;
        for(int k = 0; k < c.size(); k++)
            S += c[k]*basisFunction_x(k, x);
        return S; 
    }
    double exact(double x) // exact for a = 1
    {
        return 0.5*(exp(x) + exp(-x)) - 0.5*(exp(1.0) + exp(-1.0)); // analytical solution
    }
    double RHS_int(int k) // scalar product of RHS and k-th basis function (we use scalar product of the form (f(x), g(x)) = \int_{-1}^1 f(x)*g(x)*dx)
    {
        double S = 0;
        int M = 1024; // number of points used for numerical integration
        for(int i = 0; i < M; i++)
        {
            double x = -1 + 2*(0.5 + i)/double(M);
            S += basisFunction(k, x)*sqrt(1 + a*sqr(solution_x(x)))*(2/double(M));
        }
        return S;
    }
    void iterate() // one iteration that updates the coefficients: new coefficients are computed using the current coefficients substituted into RHS
    {
        vector<double> c1(c); // new coefficients
        for(int i = 0; i < c.size(); i++)
            c1[i] = -RHS_int(i)/sqr(0.5*M_PI*(2*i + 1));
        for(int i = 0; i < c.size(); i++)
            c[i] = c1[i]; 
    }
    void saveSolution(string fileName)
    {
        FILE * pFile;
        pFile = fopen (fileName.c_str(), "wb");
        double x, y;
        int M = 128; // number of points
        for(int i = 0; i < M; i++)
        {            
            x = -1 + 2*(0.5 + i)/double(M);
            y = solution(x);
            fwrite (&y, sizeof(double), 1, pFile);
            y = exact(x);
            fwrite (&y, sizeof(double), 1, pFile);
        }
        fclose (pFile);
    }
    double L2Error()
    {
        double S = 0;
        int M = 1024; // number of points used for numerical integration
        for(int i = 0; i < M; i++)
        {
            double x = -1 + 2*(0.5 + i)/double(M);
            S += sqr(solution(x) - exact(x))*(2/double(M));
        }
        return S;
    }
};

int main()
{
    cout << "lab3" << endl;
    /*
    //============== 1. particular case =================
    catenary X(30, 1);
    cout << X.c.size() << " basis functions are used" << endl;
    for(int i = 0; i < 10; i++)
    {
        cout << "iteration no. " << i << ": L2 error = " << X.L2Error() << endl;
        X.iterate();
    }
    X.saveSolution("out.dat");


    //============== 2. convergence rates =================
    FILE * pFile;
    pFile = fopen ("convergence.dat", "wb");
    for(int K = 1; K < 30; K += 3)
    {
        catenary X(K, 1);
        for(int i = 0; i < 10; i++)
        {
            X.iterate();
            double error = X.L2Error();
            fwrite (&error, sizeof(double), 1, pFile);
        }
        cout << K << endl;
    }
    fclose (pFile);
    */
    //============== 3. extra point challenge ============
    catenary X(30, sqrt(3.0));
    cout << X.c.size() << " basis functions are used" << endl;
    for(int i = 0; i < 10; i++)
    {
        //cout << "iteration no. " << i << ": L2 error = " << X.L2Error() << endl;
        cout << "iteration no. " << i << ": y(0) = " << X.solution(0) << endl;
        X.iterate();
    }
};