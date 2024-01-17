#include <iostream>
#include <cmath>
#include <math.h>
#include "cpgplot.h"
#include <cstdlib>

// declare functions
double rk4(int steps);
double f_x( double x, double y, double z, double a, double K, double beta_xy, double beta_xz);
double f_y( double x, double y, double b, double delta_y);
double f_z( double x, double z, double c, double delta_z);

// main function that calls solution
int main(){
rk4(100);

return 0;
}


// Define the differential equations
double f_x( double x, double y, double z, double a, double K, double beta_xy, double beta_xz) {
    return a * x * (1 - (x / K)) - (beta_xy * x * y) - (beta_xz * x * z) ;
}

double f_y( double x, double y, double b, double delta_y) {
    return -b * y + (x * y * delta_y);
}

double f_z( double x, double z, double c, double delta_z) {
    return -c * z + (x * z * delta_z);
}



double rk4(int steps){
// parameters
    double a = 1.1;     // growth rate of prey
    double b = 0.4;     // death rate of predator (per capita) Y
    double c = 0.4;    // " Z
    double K = 1;     // environmental capacity
    double beta_xy = 0.3;  // predation rate 
    double beta_xz = 0.4;  // (often intepretted as prey death rates, as a result of a particular predator) 

    double delta_z = 0.09; // effect of prey prescence on predator growth rate/ (aka predator growth rate)
    double delta_y = 0.15;
	
// user input delta (for testing
    double delta;
    //std::cout << "delta: \n";
    //std::cin >> delta;
	//delta_y = delta; delta_z = delta;

//initial values
    double x = 0.1; double y = 0.2; double z = 0.2;

// time
    double tstart = 0.0;
    double tend = 10.0;
    double d = (tend - tstart)/steps;
    double t = tstart;

	// arrays for plot
    float* tp = (float *)calloc(steps +1, sizeof(float));
    float* yp = (float *)calloc(steps +1, sizeof(float));
    float* xp = (float *)calloc(steps +1, sizeof(float));
	float* zp = (float *)calloc(steps +1, sizeof(float));

        int n = steps;
    for (int i = 0; i < n; i++){
        tp[i] = t;
        yp[i] = y;
        xp[i] = x;
		zp[i] = z;

// iterative rk4
	double k1x = d * f_x( x, y, z, a, K, beta_xy, beta_xz);
	double k1y = d * f_y( x, y, b, delta_y);
	double k1z = d * f_z( x, z, c, delta_z);

	double k2x = d * f_x( x + 0.5 * k1x, y + 0.5 * k1y, z + 0.5 * k1z, a, K, beta_xy, beta_xz);
	double k2y = d * f_y( x + 0.5 * k1x, y + 0.5 * k1y, b, delta_y);
	double k2z = d * f_z( x + 0.5 * k1x, z + 0.5 * k1z, c, delta_z);

	double k3x = d * f_x( x + 0.5 * k2x, y + 0.5 * k2y, z + 0.5 * k2z, a, K, beta_xy, beta_xz);
	double k3y = d * f_y( x + 0.5 * k2x, y + 0.5 * k2y, b, delta_y);
	double k3z = d * f_z( x + 0.5 * k2x, z + 0.5 * k2z, c, delta_z);

	double k4x = d * f_x( x + k3x, y + k3y, z + k3z, a, K, beta_xy, beta_xz);
	double k4y = d * f_y( x + k3x, y + k3y, b, delta_y);
	double k4z = d * f_z( x + k3x, z + k3z, c, delta_z);
	

        t = t + d;
        x = x + (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
        y = y + (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
		z = z + (k1z + 2 * k2z + 2 * k3z + k4z) / 6;
        }
        tp[n] = t; yp[n] = y; xp[n] = x;  zp[n] = z; 

//plot
if (!cpgopen("/XWINDOW")) return 1;

cpgenv(0., tend, 0., 1.5, 0, 1);
cpglab("t", "x,y,z", "G:prey(z), W:pred(y), B:pred(z)");

cpgsci(1); // White: predator
cpgline(n+1, tp, yp);


cpgsci(9); // Green: prey
cpgline(n+1, tp, xp);

cpgsci(5); // blue: prey2
cpgline(n+1, tp, zp);

cpgclos();


return 0.0;
}

