#include <iostream>
#include <omp.h>

#include "../include/LUDecomposition.h"


using namespace std;

int main() {
	// Initial block
	int m = 1000;  // Number of rows
	int n = 1000;  // Number of cols
	double** A = createRandomMatrix(m, n);
	/*cout << "Matrix A is" << endl;
	matrixPrint(A, m, n);*/
	double** B = getCopy(A, m, n);
	/*cout << "Matrix B is" << endl;
	matrixPrint(B, m, n);*/
	double** C = getCopy(A, m, n);
	/*cout << "Matrix C is" << endl;
	matrixPrint(C, m, n);*/
	cout << endl << endl;


	cout << "Usual LU-Decomposition" << endl;
	double inTime = omp_get_wtime();
	double** newA1 = LUDecomposition(A, m, n);
	double outTime = omp_get_wtime();
	/*cout << "Result is" << endl;
	matrixPrint(newA1, m, n);*/
	cout << "Time spent: " << outTime - inTime << endl;
	cout << endl << endl;


	cout << "Parallel usual LU-Decomposition" << endl;
	inTime = omp_get_wtime();
	double** newA2 = LUDecompositionParal(B, m, n);
	outTime = omp_get_wtime();
	/*cout << "Result is" << endl;
	matrixPrint(newA2, m, n);*/
	cout << "Time spent: " << outTime - inTime << endl;
	cout << endl << endl;


	cout << "Block LU-Decomposition" << endl;
	inTime = omp_get_wtime();
	double** newA3 = LUBlockDecomposition(C, m);
	outTime = omp_get_wtime();
	/*cout << "Result is" << endl;
	matrixPrint(newA3, m, n);*/
	cout << "Time spent: " << outTime - inTime << endl;
	cout << endl << endl;


	// Comparing
	if (compareMatrices(newA1, m, n, newA3, m, n)) {
		cout << "Matrixes are equal" << endl;
	}
	else {
		cout << "Matrixes are NOT equal" << endl;
	}
	system("pause");
}