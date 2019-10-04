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
	}
	return _A;
}

double** LUDecompositionParal(double** _A, int m, int n) {
	if (m < n) {
		cout << "Wrong sizes in LUBlockDecomposition" << endl;
		return nullptr;
	}
	int iLim = min(m - 1, n);
//#pragma omp parallel //for num_threads(4)
	//{
//#pragma omp for //schedule(static, 4)
		for (int i = 0; i < iLim; ++i)
		{
			//printf("Thread on %d is %d\n", i, omp_get_thread_num());
			//cout << "Thread on " << i << " is " << omp_get_thread_num() << endl;
//#pragma omp for
#pragma omp parallel for
			for (int k = i + 1; k < m; ++k) {
				//printf("Thread1 on [%d][%d] is %d\n", k, i, omp_get_thread_num());
				_A[k][i] /= _A[i][i];
			}
			if (i < n) {
#pragma omp parallel for
				for (int k = i + 1; k < m; ++k) {
					for (int l = i + 1; l < n; ++l) {
						//printf("Thread2 on [%d][%d] is %d\n", k, l, omp_get_thread_num());
						_A[k][l] -= _A[k][i] * _A[i][l];
					}
				}
			}
		}
		//cout << "Number of threads is " << omp_get_num_threads() << endl;
	//}
	return _A;
}

double** LUBlockDecomposition(double** _A, int n) {
	int b = 32;
	double** res = new double*[n];
	for (int p = 0; p < n; ++p) {
		res[p] = new double[n];
	}
	int _i = 0;
	int lastsize = 0;
	int length;
	double** block = new double* [b];
	for (int p = 0; p < b; ++p) {
		block[p] = new double[b];
	}
	for (int i = 0; i < n; i += b) {
		length = n - i - b;
		double** block1 = new double*[b];
		for (int p = 0; p < b; ++p) {
			block1[p] = new double[length];
		}
		double** block2 = new double*[length];
		for (int p = 0; p < length; ++p) {
			block2[p] = new double[b];
		}
		for (int k = 0; k < b; ++k) {
			for (int l = 0; l < b; ++l) {
				block[k][l] = _A[k + i][l + i];
			}
		}
		for (int k = 0; k < b; ++k) {
			for (int l = 0; l < length; ++l) {
				block1[k][l] = _A[k + i][l + i + b];
			}
		}

		for (int k = 0; k < length; ++k) {
			for (int l = 0; l < b; ++l) {
				block2[k][l] = _A[k + i + b][l + i];
			}
		}
		LUDecomposition(block, b, b);
		for (int k = i; k < i + b; ++k) {
			for (int l = i; l < i + b; ++l) {
				res[k][l] = block[k - i][l - i];
			}
		}
		block1 = linSolveDown(block, block1, b, length);  // Counting U
		block2 = linSolveUp(block, block2, length, b);   // Counting L
		for (int k = i; k < b + i; ++k) {
			for (int l = i + b; l < n; ++l) {
				res[l][k] = block2[l - i - b][k - i];
				res[k][l] = block1[k - i][l - i - b];
			}
		}
		double** mm = matrixMult(block2, block1, length, b, length);
		for (int k = 0; k < length; ++k) {
			for (int l = 0; l < length; ++l) {
				_A[i + b + k][i + b + l] -= mm[k][l];
			}
		}
		
		_i = i;
		if (i == n - b) {
			for (int k = 0; k < b; ++k) {
				for (int l = 0; l < b; ++l) {
					block[k][l] = _A[n - b + k][n - b + l];
				}
			}
			LUDecomposition(block, b, b);
			for (int i = n - b; i < n; ++i) {
				for (int j = n - b; j < n; ++j) {
					res[i][j] = block[i - n + b][j - n + b];
				}
			}
		}
		deletePointMatr(block1, b);
		deletePointMatr(block2, length);
		deletePointMatr(mm, length);
	}
	deletePointMatr(block, b);
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
	double epsNull = 1e-4;
	if (m1 != m2 || n1 != n2) {
		return false;
	}
	for (int i = 0; i < m1; ++i) {
		for (int j = 0; j < n1; ++j) {
			if (abs(_source1[i][j] - _source2[i][j]) > epsNull) {
				//cout << "i = " << i << "; j = " << j << endl;
				//printf("_source1[i][j] = %f; _source2[i][j] = %f\n", _source1[i][j], _source2[i][j]);
				/*cout << "_source1[i][j] = " << _source1[i][j] << endl;
				cout << "_source2[i][j] = " << _source2[i][j] << endl;*/
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

double** linSolveDown(double** _A, double** _b, int n, int m) {
	double** res = new double*[n];
	for (int i = 0; i < n; ++i) {
		res[i] = new double[m];
	}
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < m; ++j){
			res[i][j] = _b[i][j];
			for (int k = 0; k < i; k++) {
				res[i][j] -= _A[i][k] * res[k][j];
			}
		}
	}
	return res;
}

//double** linSolveUp(double** _A, double** _b, int n, int m) {
//	double e;
//	double** res = new double*[n];
//	for (int i = 0; i < n; ++i) {
//		res[i] = new double[m];
//	}
//	for (int i = 0; i < m; ++i){
//			e = 1 / _A[i][i];
//		for (int j = 0; j < m; ++j){
//			_A[i][j] *= e;
//
//		}
//		for (int k = 0; k < i; ++k) {
//			_A[k + 1][k] = -e * _A[k + 1][k];
//			_A[i - 1][j] -= _A[i][j] * _A[i - 1][j];
//		}
//		_A[i][i] = e;
//	}
//	cout << "Inv" << endl;
//	matrixPrint(_A, m, m);
//	res = matrixMult(_b, _A, n, m, m);
//	//_b = res;
//	return res;
//}

double** linSolveUp(double** _A, double** _b, int n, int m) {
	double** res = new double*[n];
	for (int i = 0; i < n; ++i) {
		res[i] = new double[m];
	}
	double** ident = new double*[m];
	for (int i = 0; i < m; ++i) {
		ident[i] = new double[m];
		for (int j = 0; j < m; ++j) {
			ident[i][j] = 0;
		}
		ident[i][i] = 1.0;
	}
	for (int i = 0; i < m - 1; ++i) {
		double e = 1.0 / _A[i][i];
		for (int j = 0; j < m; ++j) {
			_A[i][j] *= e;
			ident[i][j] *= e;
		}
		for (int k = i + 1; k < m; ++k) {
			double coeff = 1.0 / _A[k][k];
			coeff *= _A[i][k];
			for (int j = i + 1; j < m; ++j) {
				
				ident[i][j] -= ident[k][j] * coeff;
				_A[i][j] -= _A[k][j] * coeff;
			}
		}
	}
	double eLast = 1.0 / _A[m - 1][m - 1];
	ident[m - 1][m - 1] *= eLast;
	res = matrixMult(_b, ident, n, m, m);
	return res;
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

double** getCopy(double** _source, int m, int n) {
	double** res = new double*[m];
	for (int i = 0; i < m; ++i) {
		res[i] = new double[n];
		for (int j = 0; j < n; ++j) {
			res[i][j] = _source[i][j];
		}
	}
	return res;
}

void partPrint(double** _source) {
	for (int i = 508; i < 512; ++i) {
		for (int j = 508; j < 512; ++j) {
			cout << _source[i][j] << " ";
		}
		cout << endl;
	}
}

void deletePointMatr(double** _source, int m) {
	//cout << "Deleting matrix" << endl;
	//matrixPrint(_source, m, 32);
	for (int i = 0; i < m; ++i) {
		//cout << "i = " << i << endl;
		delete[] _source[i];
	}
	delete[] _source;
}
