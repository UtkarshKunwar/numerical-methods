#include "gauss_jordan.hpp"

// Returns identity matrix of size n
double** eye(int n)
{
	double** identity = new double*[n];
	for(int i = 0; i < n; i++) {
		identity[i] = new double[n]();	// Initialises all values to 0 simultaneously.
		identity[i][i] = 1.0;
	}

	return identity;
}

// Pivots and swaps the rows depending on the max element in that row.
int pivot_and_swap(double** A, double** I, int n)
{
	for(int j = 0; j < n; j++) {
		int current_index = j;
		int max_row_index = j;
		double max_val = abs(A[j][j]);

		// Find new max element in current column.
		for(int i = j; i < n; i++) {
			if(abs(A[i][j]) > max_val) {
				max_val = A[i][j];
				max_row_index = i;
			}
		}

		// Swap if new max element is found.
		if(current_index != max_row_index) {
			// Exchanging matrix rows.
			double* temp_mat = A[current_index];
			A[current_index] = A[max_row_index];
			A[max_row_index] = temp_mat;

			// Exchanging I matrix rows.
			double* temp_I = I[current_index];
			I[current_index] = I[max_row_index];
			I[max_row_index] = temp_I;
		}
	}

	return 0;
}

// Eliminates the elements above and below the diagonal elements.
int elimination(double** A, double** I, int n)
{
	for(int j = 0; j < n; j++) {
		if(abs(A[j][j]) < 1e-9) {
			return 1;
		}
		for(int i = 0; i < n; i++) {
			if(i != j) {
				double factor = A[i][j] / A[j][j];
				for(int k = 0; k < n; k++) {
					A[i][k] -= factor * A[j][k];
					I[i][k] -= factor * I[j][k];
				}
			}
		}
	}

	return 0;
}

// Scales matrix depending on the diagonal element of reduced A.
int scale(double** A, double** I, int n)
{
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			I[i][j] /= A[i][i];
			A[i][j] /= A[i][i];
		}
	}

	return 0;
}

// Returns the inverse of a matrix.
double** matrixInverse(double** A, int n)
{
	double** I = eye(n);
	pivot_and_swap(A, I, n);
	if(elimination(A, I, n) == 1)
		return NULL;
	scale(A, I, n);

	return I;
}
