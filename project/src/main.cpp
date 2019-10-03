#include <iostream>
#include <omp.h>

#include "../include/LUDecomposition.h"


using namespace std;

int main() {
	// Initial block
	int m = 1024;  // Number of rows
	int n = 1024;  // Number of cols
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
	cout << endl;
	partPrint(newA1);
	cout << endl;
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








	/*double** R = createRandomMatrix(3, 3);
	R[0][0] = 1;
	R[1][1] = 1;
	R[2][2] = 1;
	R[0][1] = 0;
	R[0][2] = 0;
	R[1][2] = 0;
	cout << "R is" << endl;
	matrixPrint(R, 3, 3);
	double** b = createRandomMatrix(3, 4);
	cout << "b is" << endl;
	matrixPrint(b, 3, 4);
	matrixPrint(linSolveDown(R, b, 3, 4), 3, 4);*/





	//double** R = createRandomMatrix(3, 3);
	//R[1][0] = 0;
	//R[2][0] = 0;
	//R[2][1] = 0;
	//cout << "R is" << endl;
	//matrixPrint(R, 3, 3);
	//double** b = createRandomMatrix(4, 3);
	//cout << "b is" << endl;
	//matrixPrint(b, 4, 3);
	////linSolveUp(R, b, 5, 4);
	//matrixPrint(linSolveUp(R, b, 4, 3), 4, 3);
	system("pause");
}