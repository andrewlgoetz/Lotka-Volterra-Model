#include <iostream>
#include <cmath>
#include <math.h>
#include "cpgplot.h"
#include <cstdlib>


double rk4(int steps);
double f_x( double x, double y, double a, double K, double beta);
double f_y( double x, double y, double b, double delta);

int main(){
rk4(100); // number of steps

return 0;
}


// Define the differential equations
double f_x( double x, double y, double a, double K, double beta) {
    return a * x * (1 - (x / K)) - (beta * x * y);}


// Alternate f_x without carrying capacity, used in verification
//double f_x( double x, double y, double a, double K, double beta) {
  //  return a * x *1/* (1 - (x / K))*/ - (beta * x * y); }

double f_y( double x, double y, double b, double delta) {
    return -b * y + (x * y * delta);
}




double rk4(int steps){

    double a = 1.1;     // growth rate of prey
    double b = 0.4;     // death rate of predator (per capita) 
    double K = 1;     // environmental capacity
    double beta = 0.4;//0.4;  // predation rate 
    double delta = 0.5; // effect of prey prescence on pred/ predator effeciency
			// This variable may be the focus of my experiments
			// And demo.
    // Notes from DEMO:
	// How does changing predator effeciency (effect of prey presceence on growth rate) change the curve?
	// Consider these values: 0.4
 	// 0.65
   
//user input
 //std::cout << "delta: \n";
   // std::cin >> delta;
    
    double x = 0.1; double y = 0.2;

    double tstart = 0.0;
    double tend = 10.0;
    double d = (tend - tstart)/steps;
    double t = tstart;

    float* tp = (float *)calloc(steps +1, sizeof(float));
    float* yp = (float *)calloc(steps +1, sizeof(float));
    float* xp = (float *)calloc(steps +1, sizeof(float));

	int n = steps;
    for (int i = 0; i < n; i++){
	tp[i] = t;
	yp[i] = y;
	xp[i] = x;
 
	double k1x = d * f_x( x, y, a, K, beta);
        double k1y = d * f_y( x, y, b, delta);

        double k2x = d * f_x( x + 0.5 * k1x, y + 0.5 * k1y, a, K, beta);
        double k2y = d * f_y( x + 0.5 * k1x, y + 0.5 * k1y, b, delta);

        double k3x = d * f_x( x + 0.5 * k2x, y + 0.5 * k2y, a, K, beta);
        double k3y = d * f_y( x + 0.5 * k2x, y + 0.5 * k2y, b, delta);

        double k4x = d * f_x( x + k3x, y + k3y, a, K, beta);
        double k4y = d * f_y( x + k3x, y + k3y, b, delta);

	t = t + d;
	x = x + (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
        y = y + (k1y + 2 * k2y + 2 * k3y + k4y) / 6;

	}
	tp[n] = t; yp[n] = y; xp[n] = x;


if (!cpgopen("/XWINDOW")) return 1;

cpgenv(0., tend, 0., 1.2, 0, 1);

cpglab("t", "y", "rk4");

cpgsci(1); // White: predator
cpgline(n+1, tp, yp);

cpgsci(9); // Green: prey
cpgline(n+1, tp, xp);

cpgclos();


return 0.0;
}

