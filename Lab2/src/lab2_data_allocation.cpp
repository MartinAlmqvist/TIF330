#include <iostream>
#include <cmath>

using namespace std;

struct grid3
{
    double **v;
    int n;
    double null;
    grid3(int n): n(n) {
        v = new double*[n];
        for(int ia = 0; ia < n; ia++) v[ia] = new double[2*(n - ia) - 1];
        null = 0;
    }
    double operator()(int ia, int ib){        
        if(ia == -1) return -v[0][ib-1]; // for boundary conditions
        if(ib == -1) return -v[ia][0]; // for boundary conditions
        if(ib == 2*(n - ia) - 1) return -v[ia][2*(n - ia) - 2]; // for boundary conditions
        return v[ia][ib];
    }
    ~grid3(){
        for(int i = 0; i < n; i++) delete []v[i];
        delete []v;
    }
};

void doStepTrapezoidalDrum(grid3 &u_p, grid3 &u, int n, double h, double tau, int a_max){
 
	//=====the following can be used to set the boundary conditions
    if(a_max < n)
        for(int ib = 0; ib < 2*(n - a_max) - 1; ib+=2)
            u.v[a_max][ib] = -u.v[a_max - 1][ib + 1];

    #pragma omp parallel for collapse(2)
    for(int ia = 0; ia < a_max; ia++)
    for(int ib = 0; ib < 2*(n - ia) - 1; ib++){
        
        double l = h;
        double A = sqrt(3.0)/4.0 * l*l;
        double dr = h /sqrt(3.0);
        if(ib%2 == 0){
            u_p.v[ia][ib] = tau*l/(dr*A)*(u.operator()(ia,ib-1) + u.operator()(ia,ib+1) + u.operator()(ia-1,ib+1) -3*u.operator()(ia,ib)) + u.operator()(ia,ib);
        }
        else{
            u_p.v[ia][ib] = tau*l/(dr*A)*(u.operator()(ia,ib-1) + u.operator()(ia,ib+1) + u.operator()(ia+1,ib-1) -3*u.operator()(ia,ib)) + u.operator()(ia,ib);   
        }
    }
		
    
}

void trapezoidalDrum(double x){
    int n = 128;
    int a_max = round(n*x);
    double t1 = 0;
    double t2 = 0;
    double u_1 = 0;
    double u_2 = 0;

    grid3 u0(n), u1(n); 

    //setting initial condition (any non-zero would work)
    for(int ia = 0; ia < n; ia++)
        for(int ib = 0; ib < 2*(n - ia) - 1; ib++)
            u0.v[ia][ib] = 1;

    double tau = 0.05/double(n*n);
    int stride = 10;
    for(int i = 0; i < 100000; i++){
        doStepTrapezoidalDrum(u1, u0, n, 1/double(n), tau, a_max);
        doStepTrapezoidalDrum(u0, u1, n, 1/double(n), tau, a_max);
        if(i == 98000){
            t1 = tau*i*2;
            u_1 =u0(int(n/3.0), int(2*n/3.0));  
        }
		if(i%1000 == 0) cout << "t = " << 2*tau*i << ": u(center) = " << u0(int(n/3.0), int(2*n/3.0)) << endl;
        if(i == 99000){
            t2 = tau*i*2;
            u_2 = u0(int(n/3.0), int(2*n/3.0));
        }

    }
  double lamdba = log(u_1/u_2)/(t2-t1);
  double omega = sqrt(lamdba);

  cout << "lamdba = " << lamdba << "omega = " << omega << endl;
};

int main(){
    trapezoidalDrum(1.0);
}
