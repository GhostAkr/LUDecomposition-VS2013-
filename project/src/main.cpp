#include <iostream>
#include <omp.h>

#include "../include/LUDecomposition.h"


using namespace std;

int main() {
	// Initial block
	int m = 1024;  // Number of rows
	int n = 1024;  // Number of cols
	double** initialMatrix = createRandomMatrix(m, n);
	/*cout << endl;
	cout << "Initial matrix is" << endl;
	matrixPrint(initialMatrix, m, m);
	cout << endl;*/
	double** A = getCopy(initialMatrix, m, n);
	/*cout << "Matrix A is" << endl;
	matrixPrint(A, m, n);*/
	double** B = getCopy(A, m, n);
	/*cout << "Matrix B is" << endl;
	matrixPrint(B, m, n);*/
	double** C = getCopy(A, m, n);
	/*cout << "Matrix C is" << endl;
	matrixPrint(C, m, n);*/
	double** D = getCopy(A, m, n);
	//cout << endl << endl;


	cout << "Usual LU-Decomposition" << endl;
	double inTime = omp_get_wtime();
	double** newA1 = LUDecomposition(A, m, n);
	double outTime = omp_get_wtime();
	/*cout << "Result is" << endl;
	matrixPrint(newA1, m, n);*/
	/*cout << endl;
	partPrint(newA1);
	cout << endl;*/
	cout << "Difference with A = " << checkLU(initialMatrix, newA1, m) << endl;
	cout << "Time spent: " << outTime - inTime << endl;
	cout << endl << endl;

	


	cout << "Parallel usual LU-Decomposition" << endl;
	inTime = omp_get_wtime();
	double** newA2 = LUDecompositionParal(B, m, n);
	outTime = omp_get_wtime();
	/*cout << "Result is" << endl;
	matrixPrint(newA2, m, n);*/
	cout << "Difference with A = " << checkLU(initialMatrix, newA2, m) << endl;
	cout << "Time spent: " << outTime - inTime << endl;
	cout << endl << endl;


	cout << "Block LU-Decomposition" << endl;
	inTime = omp_get_wtime();
	double** newA3 = LUBlockDecomposition(C, m);
	outTime = omp_get_wtime();
	/*cout << "Result is" << endl;
	matrixPrint(newA3, m, n);*/
	cout << "Difference with A = " << checkLU(initialMatrix, newA3, m) << endl;
	//printf("%f", checkLU(initialMatrix, newA3, m));
	cout << "Time spent: " << outTime - inTime << endl;
	cout << endl << endl;

	
	//cout << "Parallel block LU-Decomposition" << endl;
	//inTime = omp_get_wtime();
	//double** newA4 = LUBlockDecompositionParal(D, m);
	//outTime = omp_get_wtime();
	///*cout << "Result is" << endl;
	//matrixPrint(newA3, m, n);*/
	//cout << "Difference with A = " << checkLU(initialMatrix, newA4, m) << endl;
	//cout << "Time spent: " << outTime - inTime << endl;
	//cout << endl << endl;


	deletePointMatr(A, m);
	deletePointMatr(B, m);
	deletePointMatr(C, m);
	deletePointMatr(D, m);
	deletePointMatr(newA3, m);
	//deletePointMatr(newA4, m);
	system("pause");
}