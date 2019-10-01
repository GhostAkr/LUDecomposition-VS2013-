#include <iostream>
#include <omp.h>

#include "../include/LUDecomposition.h"


using namespace std;

int main() {
	int m = 4;  // Number of rows
	int n = 4;  // Number of cols
	double** A = createRandomMatrix(m, n);
	//cout << "Initial matrix" << endl;
	//matrixPrint(A, m, n);
	cout << "Usual LU-Decomposition" << endl;
	double inTime = omp_get_wtime();
	double** newA1 = LUDecomposition(A, m, n);
	double** res = 
	double outTime = omp_get_wtime();
	cout << "Time spent: " << outTime - inTime << endl;
	cout << "Parallel LU-Decomposition" << endl;
	inTime = omp_get_wtime();
	double** newA2 = LUDecompositionParal(A, m, n);
	outTime = omp_get_wtime();
	cout << "Time spent: " << outTime - inTime << endl;
	if (compareMatrices(newA1, m, n, newA2, m, n)) {
		cout << "Matrixes are equal" << endl;
	}
	else {
		cout << "Matrixes are NOT equal" << endl;
	}
	system("pause");
}