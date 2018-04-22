#ifndef ODE_H
#define ODE_H

#include <cmath>
#include <vector>

// Euler
double* euler(double (*f)(double), double xi, double yi, double xf, double h);
double* euler(double (*f)(double), double xi, double yi, double xf, int n);

// Modified Euler
double* modified_euler(double (*f)(double), double xi, double yi, double xf, double h);
double* modified_euler(double (*f)(double), double xi, double yi, double xf, int n);

// Mid-Point
double* midpoint(double (*f)(double), double xi, double yi, double xf, double h);
double* midpoint(double (*f)(double), double xi, double yi, double xf, int n);

// RK3
double* RK3(double (*f)(double), double xi, double yi, double xf, double h);
double* RK3(double (*f)(double), double xi, double yi, double xf, int n);

// Adam's-Bashforth
double* adams_bashforth(double (*f)(double), double xi, double yi, double xf, double h);
double* adams_bashforth(double (*f)(double), double xi, double yi, double xf, int n);

// Adam's-Bashforth
double* milne_simpsons(double (*f)(double), double xi, double yi, double xf, double h);
double* milne_simpsons(double (*f)(double), double xi, double yi, double xf, int n);

// RK4
double* RK4(double (*f)(double), double xi, double yi, double xf, double h);
double* RK4(double (*f)(double), double xi, double yi, double xf, int n);

int get_number_of_steps(double xi, double xf, double h);
double get_step_size(double xi, double xf, int n);

#endif
