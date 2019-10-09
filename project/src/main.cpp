#include "../../rowMatrix.h"
#include <iostream>
#include <omp.h>

#include "../include/LUDecomposition.h"



using namespace std;

int main() {

	// Initial block
	int m = 1024;  // Number of rows
	int n = 1024;  // Number of cols
	int b = 32;
	double* initialMatrix = createRandomRowMatrix(m, n);
	double* A = getCopy(initialMatrix, m, n);
	double* B = getCopy(initialMatrix, m, n);
	double* C = getCopy(initialMatrix, m, n);
	double* D = getCopy(initialMatrix, m, n);


	cout << "Usual LU-Decomposition" << endl;
	double inTime = omp_get_wtime();
	double* newA1 = LUDecomposition(A, m, n);
	double outTime = omp_get_wtime();
	cout << "Difference with A = " << checkLU(initialMatrix, newA1, m) << endl;
	cout << "Time spent: " << outTime - inTime << endl;
	cout << endl << endl;


	cout << "Parallel usual LU-Decomposition" << endl;
	inTime = omp_get_wtime();
	double* newA2 = LUDecompositionParal(B, m, n);
	outTime = omp_get_wtime();
	cout << "Difference with A = " << checkLU(initialMatrix, newA2, m) << endl;
	cout << "Time spent: " << outTime - inTime << endl;
	cout << endl << endl;


	cout << "Block LU-Decomposition" << endl;
	inTime = omp_get_wtime();
	double* newA3 = LUBlockDecomposition(C, m, b);
	outTime = omp_get_wtime();
	cout << "Difference with A = " << checkLU(initialMatrix, newA3, m) << endl;
	cout << "Time spent: " << outTime - inTime << endl;
	cout << endl << endl;

	
	cout << "Parallel block LU-Decomposition" << endl;
	inTime = omp_get_wtime();
	double* newA4 = LUBlockDecompositionParal(D, m, b);
	outTime = omp_get_wtime();
	cout << "Difference with A = " << checkLU(initialMatrix, newA4, m) << endl;
	cout << "Time spent: " << outTime - inTime << endl;
	cout << endl << endl;


	deletePointMatr(A, m);
	deletePointMatr(B, m);
	deletePointMatr(C, m);
	deletePointMatr(D, m);
	deletePointMatr(newA3, m);
	deletePointMatr(newA4, m);
	system("pause");
}