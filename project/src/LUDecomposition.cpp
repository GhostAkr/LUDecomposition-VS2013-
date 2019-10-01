#include <iostream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <omp.h>

#include "../include/LUDecomposition.h"

using namespace std;

// Main methods

double** LUDecomposition(double** _A, int m, int n) {
	if (m < n) {
		cout << "Wrong sizes in LUBlockDecomposition" << endl;
		return nullptr;
	}
	int iLim = min(m - 1, n);
	for (int i = 0; i < iLim; ++i) {
		for (int k = i + 1; k < m; ++k) {
			_A[k][i] /= _A[i][i];
		}
		if (i < n) {
			for (int k = i + 1; k < m; ++k) {
				for (int l = i + 1; l < n; ++l) {
					_A[k][l] -= _A[k][i] * _A[i][l];
				}
			}
		}
		//cout << "Matrix A on " << i << " iteration is" << endl;
		//matrixPrint(_A, m, n);
	}
	return _A;
}

double** LUDecompositionParal(double** _A, int m, int n) {
	if (m < n) {
		cout << "Wrong sizes in LUBlockDecomposition" << endl;
		return nullptr;
	}
	int iLim = min(m - 1, n);
#pragma omp parallel num_threads(4)
	{
#pragma omp for
		for (int i = 0; i < iLim; ++i) {
			for (int k = i + 1; k < m; ++k) {
				_A[k][i] /= _A[i][i];
			}
			if (i < n) {
				for (int k = i + 1; k < m; ++k) {
					for (int l = i + 1; l < n; ++l) {
						_A[k][l] -= _A[k][i] * _A[i][l];
					}
				}
			}
			//cout << "Matrix A on " << i << " iteration is" << endl;
			//matrixPrint(_A, m, n);
		}
		cout << "Number of threads is " << omp_get_num_threads() << endl;
	}
	return _A;
}

double** LUBlockDecomposition(double** _A, int n) {
	int b = 30;
	double** res = new double*[b];
	for (int p = 0; p < n; ++p) {
		res[p] = new double[n];
	}

	//int b = 2;
	int _i = 0;
	int lastsize = 0;
	int lenght; // TODO: написать правильно, если написано неправильно, но а если правильно, то можно вообще не трогать, но тут на вкус и цвет, как говориться.
		double** block = new double* [b];
		for (int p = 0; p < b; ++p) {
			block[p] = new double[b];
		}

		for (int i = 0; i < n - 1 - b; i += b) {
		lenght = n - b * (i + 1);
		double** block1 = new double*[b];
		for (int p = 0; p < b; ++p) {
			block1[p] = new double[lenght];
		}
		double**block2 = new double*[b];
		for (int p = 0; p < lenght; ++p) {
			block2[p] = new double[b];
		}


		for (int k = 0; k < b; ++k) {
			for (int l = 0; l < b; ++l) {
				block[k][l] = _A[k + i][l + i];
			}
		}
		LUDecompositionParal(block, b, b);  // TODO: make it void
		linSolveDown(block, block1, b, lenght);// сначала считам U
		linSolveUp(block, block2, b, lenght);  // потом L. матрица block портится 
		
		for (int k = 0; k < b; ++k) {
			for (int l = 0; l < n - i * b; ++l) {
				res[l][k] = block2[k][l];  //сначала записываем L тк на диагонали 1 и их не жалко перекрыть значениями из U
				res[k][l] = block1[k][l];
			}
		}
		_i = i;
		lastsize = n - _i * b;
		double** last = new double*[lastsize]; // для общего случая. если последняя матрица такого же размера как первая, можно использовать block
		for (int p = 0; p < lastsize; ++p) {
			last[p] = new double[lastsize];
		}
		LUDecomposition(last, lastsize, lastsize);
		for (int k = 0; k < lastsize; ++k) {
			for (int l = 0; l < lastsize; ++l) {
				res[k][l] = last[k][l];
			}
		}
	}
		
	return res;
}

// Other methods

void matrixPrint(double** _source, int m, int n) {
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			cout << _source[i][j] << " ";
		}
		cout << endl;
	}
}

double** createRandomMatrix(int m, int n) {
	double** result = new double* [m];
	for (int i = 0; i < m; ++i) {
		result[i] = new double[n];
	}
	//srand(time(NULL));
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			result[i][j] = rand() / 1000 + 5;
		}
	}
	return result;
}

bool compareMatrices(double** _source1, int m1, int n1, double** _source2, int m2, int n2) {
	if (m1 != m2 || n1 != n2) {
		return false;
	}
	for (int i = 0; i < m1; ++i) {
		for (int j = 0; j < n1; ++j) {
			if (_source1[i][j] != _source2[i][j]) {
				return false;
			}
		}
	}
	return true;
}

double** getU22(double** _source, int m, int n) {
	double** result = new double* [m];
	for (int i = 0; i < m; ++i) {
		result[i] = new double[n];
	}
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			result[i][j] = 0;
		}
	}
	for (int i = 0; i < n; ++i) {
		for (int j = i; j < n; ++j) {
			result[i][j] = _source[i][j];
		}
	}
	return result;
}

double** getL(double** _source, int m, int n) {
	double** result = new double* [m];
	for (int i = 0; i < m; ++i) {
		result[i] = new double[m];
	}
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			if (i == j) {
				result[i][j] = 1.0;
			}
			else {
				result[i][j] = 0;
			}
		}
	}
	for (int i = 1; i < m; ++i) {
		for (int j = 0; j < i; ++j) {
			result[i][j] = _source[i][j];
		}
	}
	return result;
}

double** getL22(double** _source, int m, int n) {
	double** L = getL(_source, m, n);
	double** result = new double* [m / 2];
	for (int i = 0; i < m / 2; ++i) {
		result[i] = new double[m];
	}
	for (int i = 0; i < m / 2; ++i) {
		for (int j = 0; j < m; ++j) {
			result[i][j] = L[i][j];
		}
	}
	return result;
}

double** getL32(double** _source, int m, int n) {
	double** L = getL(_source, m, n);
	double** result = new double* [m / 2];
	for (int i = 0; i < m / 2; ++i) {
		result[i] = new double[m];
	}
	for (int i = m / 2; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			result[i][j] = L[i][j];
		}
	}
	return result;
}

void linSolveDown(double** _A, double** _b, int n, int m) {
	double** res = new double*[n];
	for (int i = 0; i < n; ++i) {
		res[i] = new double[m];
	}
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < m; ++j){
			res[i][j] = _b[i][j];
			for (int k = 0; k < i; k++)
				res[i][j] -= _A[i-1][k] * res[k][j];
		}
	}
	_b = res;
}

void linSolveUp(double** _A, double** _b, int n, int m) {
	double e;
	double** res = new double*[n];
	for (int i = 0; i < n; ++i) {
		res[i] = new double[m];
	}
	for (int i = 0; i < n; ++i){
			e = 1 / _A[i][i];
		for (int j = 0; j < i; ++j){
			_A[i][j] = -e * _A[i][j];
		}
		_A[i][i] = 1;
	}
	//res = matrixMult(_A, _b, );
	_b = res;
}

double** matrixMult(double** _source1, double** _source2, int m, int n, int s) {
	double** result = new double*[m];
	for (int i = 0; i < m; ++i) {
		result[i] = new double[s];
		for (int j = 0; j < s; ++j) {
			result[i][j] = 0.0;
		}
	}

	for (int k = 0; k < n; k++) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < s; j++) {
				result[i][j] += _source1[i][k] * _source2[k][j];
			}
		}
	}
	return result;
}