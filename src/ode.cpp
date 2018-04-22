#include "ode.hpp"
#include <iostream>

int get_number_of_steps(double xi, double xf, double h)
{
	return (int)((xf - xi) / h);
}

double get_step_size(double xi, double xf, int n)
{
	return (xf - xi) / n;
}

double* euler(double (*f)(double), double xi, double yi, double xf, double h)
{
	int n = get_number_of_steps(xi, xf, h);
	return euler(*f, xi, yi, xf, n);
}

double* euler(double (*f)(double), double xi, double yi, double xf, int n)
{
	double* y = new double[n + 1];

	y[0] = yi;

	double h = get_step_size(xi, xf, n);
	for(int i = 0; i < n; i++) {
		y[i + 1] = y[i] + h * (*f)(y[i]);
	}

	return y;
}

double* modified_euler(double (*f)(double), double xi, double yi, double xf, double h)
{
	int n = get_number_of_steps(xi, xf, h);
	return modified_euler(*f, xi, yi, xf, n);
}

double* modified_euler(double (*f)(double), double xi, double yi, double xf, int n)
{
	double* y = new double[n + 1];

	y[0] = yi;

	double h = get_step_size(xi, xf, n);
	for(int i = 0; i < n; i++) {
		// y_i+1 = y_i + h/2 * (f(y_i) + f(y_i+1)
		y[i + 1] = y[i] + 0.5 * h * ((*f)(y[i]) + (*f)(y[i] + h * (*f)(y[i])));
	}

	return y;
}

double* midpoint(double (*f)(double), double xi, double yi, double xf, double h)
{
	int n = get_number_of_steps(xi, xf, h);
	return midpoint(*f, xi, yi, xf, n);
}

double* midpoint(double (*f)(double), double xi, double yi, double xf, int n)
{
	double* y = new double[n + 1];

	y[0] = yi;

	double h = get_step_size(xi, xf, n);
	for(int i = 0; i < n; i++) {
		// y_i+1 = y_i + h * f(y_i + h/2 * f(y_i))
		y[i + 1] = y[i] + h * (*f)(y[i] + 0.5 * h * (*f)(y[i]));
	}

	return y;
}

double* RK3(double (*f)(double), double xi, double yi, double xf, double h)
{
	int n = get_number_of_steps(xi, xf, h);
	return RK3(*f, xi, yi, xf, n);
}

double* RK3(double (*f)(double), double xi, double yi, double xf, int n)
{
	double* y = new double[n + 1];

	y[0] = yi;

	double h = get_step_size(xi, xf, n);
	for(int i = 0; i < n; i++) {
		double k1 = (*f)(y[i]);
		double k2 = (*f)(y[i] + 0.5 * h * k1);
		double k3 = (*f)(y[i] + 0.75 * h * k2);

		y[i + 1] = y[i] + h / 6 * (k1 + 4 * k2 + k3);
	}

	return y;
}

double* adams_bashforth(double (*f)(double), double xi, double yi, double xf, double h)
{
	int n = get_number_of_steps(xi, xf, h);
	return adams_bashforth(*f, xi, yi, xf, n);
}

double* adams_bashforth(double (*f)(double), double xi, double yi, double xf, int n)
{
	double* y = new double[n + 1];

	double h = get_step_size(xi, xf, n);
	y[0] = yi;
	double* temp = euler(*f, xi, yi, xi + h, 1);
	y[1] = temp[1];
	delete temp;

	for(int i = 1; i < n; i++) {
		// y_i+1 = y_i + h/2 * (3f_i - f_i-1)
		y[i + 1] = y[i] + 0.5 * h * (3 * (*f)(y[i]) - (*f)(y[i - 1]));
	}

	return y;
}

double* milne_simpsons(double (*f)(double), double xi, double yi, double xf, double h)
{
	int n = get_number_of_steps(xi, xf, h);
	return milne_simpsons(*f, xi, yi, xf, n);
}

double* milne_simpsons(double (*f)(double), double xi, double yi, double xf, int n)
{
	double* y = new double[n + 1];

	double h = get_step_size(xi, xf, n);

	y[0] = yi;
	double* temp = euler(*f, xi, yi, xi + h, 3);
	y[1] = temp[1];
	y[2] = temp[2];
	y[3] = temp[3];

	delete temp;

	for(int i = 3; i < n; i++) {
		// Predictor Equation
		// y_i+1 = y_i-3 + 4/3 h (2f_i - f_i-1 + 2f_i-2)
		y[i + 1] = y[i - 3] + 4 * h / 3 * (2 * (*f)(y[i])- (*f)(y[i - 1]) + 2 * (*f)(y[i - 2]));

		for(int j = 0; j < 5; j++) {
			// Corrector Equation
			// y_i+1 = y_i-1 + h/3 (f_i+1 + 4f_i + f_i-1)
			y[i + 1] = y[i - 1] + h / 3 * ((*f)(y[i + 1]) + 4 * (*f)(y[i]) + (*f)(y[i - 1]));
		}
	}

	return y;
}

double* RK4(double (*f)(double), double xi, double yi, double xf, double h)
{
	int n = get_number_of_steps(xi, xf, h);
	return RK4(*f, xi, yi, xf, n);
}

double* RK4(double (*f)(double), double xi, double yi, double xf, int n)
{
	double* y = new double[n + 1];

	y[0] = yi;

	double h = get_step_size(xi, xf, n);
	for(int i = 0; i < n; i++) {
		double k1 = (*f)(y[i]);
		double k2 = (*f)(y[i] + 0.5 * h * k1);
		double k3 = (*f)(y[i] + 0.5 * h * k2);
		double k4 = (*f)(y[i] + h * k3);

		y[i + 1] = y[i] + h / 6 * (k1 + 2 * (k2 + k3) + k4);
	}

	return y;
}
