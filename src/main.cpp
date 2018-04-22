#include "ode.hpp"
#include <iostream>

double dv_dt(double v)
{
	double dvdt = 9.81 - 12.5 / 68.1 * v;
	return dvdt;
}

int main()
{
	double v0 = 44.9189;
	double t0 = 10;
	double tf = 12;
	double h = 0.25;
	double* y = RK4(&dv_dt, t0, v0, tf, h);

	for(int i = 0; i < (int)((tf - t0) / h) + 1; i++) {
		std::cout << y[i] << std::endl;
	}

	delete y;
	return 0;
}
