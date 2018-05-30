#include <iostream>
#include <iomanip>
#include <fstream>
#include "gauss_jordan.hpp"

void printUsage(char* prog_name)
{
	std::cout << "./" << prog_name << " <input file>\n" << std::endl;
	std::cout << "Example <input file>\n3\n1 2 3\n4 5 6\n7 8 9" << std::endl;
}

int main(int argc, char* argv[])
{
	if(argc < 2) {
		printUsage(argv[0]);
		return 1;
	}

	// Reading file.
	std::ifstream fin;
	fin.open(argv[1]);
	if(fin.fail()) {
		std::cout << "Error reading file." << std::endl;
		std::exit(0);
	}

	// Size of the matrix (assuming matrix is a square matrix)
	int n;
	fin >> n;

	// Reading matrix from file.
	double** A= new double*[n];
	for(int i = 0; i < n; i++) {
		A[i] = new double[n];
		for(int j = 0; j < n; j++) {
			fin >> A[i][j];
		}
	}
	fin.close();

	double** A_inverse = matrixInverse(A, n);

	if(A_inverse == NULL) {
		std::cout << "Inverse for the matrix does not exist." << std::endl;
		return 1;
	}

	// Opening file for writing.
	std::ofstream fout;
	fout.open("results.txt");
	if(fout.fail()) {
		std::cout << "Error writing to file." << std::endl;
		std::exit(0);
	}

	// Writing the inverse matrix to the output file.
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			fout << std::fixed << A_inverse[i][j] << "\t";
		}
		fout << std::endl;
	}
	fout.close();

	return 0;
}
