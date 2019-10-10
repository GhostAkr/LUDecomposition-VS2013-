#include <iostream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <omp.h>

#include "../include/LUDecomposition.h"

using namespace std;

// Main methods

double* LUDecomposition(double* _A, int m, int n) {
	if (m < n) {
		cout << "Wrong sizes in LUBlockDecomposition" << endl;
		return nullptr;
	}
	int iLim = min(m - 1, n);
	for (int i = 0; i < iLim; ++i) {
		for (int k = i + 1; k < m; ++k) {
			_A[k * n + i] /= _A[i * n + i];
		}
		if (i < n) {
			for (int k = i + 1; k < m; ++k) {
				for (int l = i + 1; l < n; ++l) {
					_A[k * n + l] -= _A[k * n + i] * _A[i * n + l];
				}
			}
		}
	}
	return _A;
}

double* LUDecompositionParal(double* _A, int m, int n) {
	if (m < n) {
		cout << "Wrong sizes in LUBlockDecomposition" << endl;
		return nullptr;
	}
	int iLim = min(m - 1, n);
	for (int i = 0; i < iLim; ++i) {
#pragma omp parallel for
		for (int k = i + 1; k < m; ++k) {
			_A[k * n + i] /= _A[i * n + i];
		}
		if (i < n) {
#pragma omp parallel for
			for (int k = i + 1; k < m; ++k) {
				for (int l = i + 1; l < n; ++l) {
					_A[k * n + l] -= _A[k * n + i] * _A[i * n + l];
				}
			}
		}
	}
	return _A;
}

double* LUBlockDecomposition(double* _A, int n, int b) {
	double* res = new double[n * n];
	int _i = 0;
	int lastsize = 0;
	int length;
	double* block = new double [b * b];
	double t1 = 0.0, t2 = 0.0;
	for (int i = 0; i < n; i += b) {
		if (n - i < b) {
			b = n - i;
		}
		length = n - i - b;
		double* block1 = new double[b * length];
		double* block2 = new double[length * b];
		for (int k = 0; k < b; ++k) {
			for (int l = 0; l < b; ++l) {
				block[k * b + l] = _A[(k + i) * n + (l + i)];
			}
		}
		for (int k = 0; k < b; ++k) {
			for (int l = 0; l < length; ++l) {
				block1[k * length + l] = _A[(k + i) * n + (l + i + b)];
			}
		}

		for (int k = 0; k < length; ++k) {
			for (int l = 0; l < b; ++l) {
				block2[k * b + l] = _A[(k + i + b) * n + (l + i)];
			}
		}
		LUDecomposition(block, b, b);
		for (int k = i; k < i + b; ++k) {
			for (int l = i; l < i + b; ++l) {
				res[k * n + l] = block[(k - i) * b + (l - i)];
			}
		}
		block1 = linSolveDown(block, block1, b, length);  // Counting U
		block2 = linSolveUp(block, block2, length, b);   // Counting L
		for (int k = i; k < b + i; ++k) {
			for (int l = i + b; l < n; ++l) {
				res[l * n + k] = block2[(l - i - b) * b + (k - i)];
				res[k * n + l] = block1[(k - i) * length + (l - i - b)];
			}
		}
		for (int k = i + b; k < n; ++k) {
			for (int l = i + b; l < n; ++l) {
				double sum = 0.0;
				for (int p = 0; p < b; ++p) {
					sum += block2[(k - i - b) * b + p] * block1[p * length + (l - i - b)];
				}
				_A[k * n + l] -= sum;
			}
		}
		_i = i;
		if (i == n - b) {
			for (int k = 0; k < b; ++k) {
				for (int l = 0; l < b; ++l) {
					block[k * b + l] = _A[(n - b + k) * n + (n - b + l)];
				}
			}
			LUDecomposition(block, b, b);
			for (int i = n - b; i < n; ++i) {
				for (int j = n - b; j < n; ++j) {
					res[i * n + j] = block[(i - n + b) * b + (j - n + b)];
				}
			}
		}
		delete[] block1;
		delete[] block2;
	}
	delete[] block;
	return res;
}

double* LUBlockDecompositionParal(double* _A, int n, int b) {
	double* res = new double[n * n];
	int _i = 0;
	int lastsize = 0;
	int length;
	double* block = new double[b * b];
	for (int i = 0; i < n; i += b) {
		if (n - i < b) {
			b = n - i;
		}
		length = n - i - b;
		double* block1 = new double[b * length];
		double* block2 = new double[length * b];
		for (int k = 0; k < b; ++k) {
			for (int l = 0; l < b; ++l) {
				block[k * b + l] = _A[(k + i) * n + (l + i)];
			}
		}
		for (int k = 0; k < b; ++k) {
			for (int l = 0; l < length; ++l) {
				block1[k * length + l] = _A[(k + i) * n + (l + i + b)];
			}
		}

		for (int k = 0; k < length; ++k) {
			for (int l = 0; l < b; ++l) {
				block2[k * b + l] = _A[(k + i + b) * n + (l + i)];
			}
		}
		LUDecomposition(block, b, b);
		for (int k = i; k < i + b; ++k) {
			for (int l = i; l < i + b; ++l) {
				res[k * n + l] = block[(k - i) * b + (l - i)];
			}
		}
			block1 = linSolveDownParal(block, block1, b, length);  // Counting U
			block2 = linSolveUpParal(block, block2, length, b);   // Counting L
		for (int k = i; k < b + i; ++k) {
			for (int l = i + b; l < n; ++l) {
				res[l * n + k] = block2[(l - i - b) * b + (k - i)];
				res[k * n + l] = block1[(k - i) * length + (l - i - b)];
			}
		}
		omp_set_num_threads(4);
#pragma omp parallel for
		for (int k = i + b; k < n; ++k) {
			for (int l = i + b; l < n; ++l) {
				double sum = 0.0;
					for (int p = 0; p < b; ++p) {
						sum += block2[(k - i - b) * b + p] * block1[p * length + (l - i - b)];
					}
				_A[k * n + l] -= sum;
			}
		}
		_i = i;
		if (i == n - b) {
			for (int k = 0; k < b; ++k) {
				for (int l = 0; l < b; ++l) {
					block[k * b + l] = _A[(n - b + k) * n + (n - b + l)];
				}
			}
			LUDecomposition(block, b, b);
			for (int i = n - b; i < n; ++i) {
				for (int j = n - b; j < n; ++j) {
					res[i * n + j] = block[(i - n + b) * b + (j - n + b)];
				}
			}
		}
		delete[] block1;
		delete[] block2;
	}
	delete[] block;
	return res;
}


// Other methods


void matrixPrint(double* _source, int m, int n) {
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			cout << _source[i * n + j] << " ";
		}
		cout << endl;
	}
}

double* createRandomRowMatrix(int m, int n) {
	double* result = new double[m * n];
	//srand(time(NULL));
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			result[i * n + j] = rand() / 1000 + 5;
		}
	}
	return result;
}

bool compareMatrices(double** _source1, int m1, int n1, double** _source2, int m2, int n2) {
	double epsNull = 1e-2;
	if (m1 != m2 || n1 != n2) {
		return false;
	}
	for (int i = 0; i < m1; ++i) {
		for (int j = 0; j < n1; ++j) {
			if (abs(_source1[i][j] - _source2[i][j]) > epsNull) {
				//printf("_source1[i][j] = %f; _source2[i][j] = %f\n", _source1[i][j], _source2[i][j]);
				return false;
			}
		}
	}
	return true;
}

double* linSolveDown(double* _A, double* _b, int n, int m) {
	double* res = new double[n * m];
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < m; ++j){
			res[i * m + j] = _b[i * m + j];
			for (int k = 0; k < i; k++) {
				res[i * m + j] -= _A[i * n + k] * res[k * m + j];
			}
		}
	}
	return res;
}

double* linSolveDownParal(double* _A, double* _b, int n, int m) {
	double* res = new double[n * m];
	omp_set_num_threads(4);
	for (int i = 0; i < n; ++i){
#pragma omp parallel for
		for (int j = 0; j < m; ++j){
			res[i * m + j] = _b[i * m + j];
			double tmp = 0.0;
				for (int k = 0; k < i; k++) {
					res[i * m + j] -= _A[i * n + k] * res[k * m + j];
				}
		}
	}
	return res;
}

double* linSolveUp(double* _A, double* _b, int n, int m) {
	double* res = new double[n * m];
	double* ident = new double[m * m];
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			ident[i * m + j] = 0;
		}
		ident[i * m + i] = 1.0;
	}
	for (int i = 0; i < m - 1; ++i) {
		double e = 1.0 / _A[i * m + i];
		for (int j = 0; j < m; ++j) {
			_A[i * m + j] *= e;
			ident[i * m + j] *= e;
		}
		for (int k = i + 1; k < m; ++k) {
			double coeff = 1.0 / _A[k * m + k];
			coeff *= _A[i * m + k];
			for (int j = i + 1; j < m; ++j) {
				
				ident[i * m + j] -= ident[k * m + j] * coeff;
				_A[i * m + j] -= _A[k * m + j] * coeff;
			}
		}
	}
	double eLast = 1.0 / _A[(m - 1) * m + (m - 1)];
	ident[(m - 1) * m + (m - 1)] *= eLast;
	res = matrixMult(_b, ident, n, m, m);
	delete[] ident;
	return res;
}

double* linSolveUpParal(double* _A, double* _b, int n, int m) {
	double* res = new double[n * m];
	double* ident = new double[m * m];
	omp_set_num_threads(4);
#pragma omp parallel for
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			ident[i * m + j] = 0;
		}
		ident[i * m + i] = 1.0;
	}
	for (int i = 0; i < m - 1; ++i) {
		double e = 1.0 / _A[i * m + i];
#pragma omp parallel for
		for (int j = 0; j < m; ++j) {
			_A[i * m + j] *= e;
			ident[i * m + j] *= e;
		}
		for (int k = i + 1; k < m; ++k) {
			double coeff = 1.0 / _A[k * m + k];
			coeff *= _A[i * m + k];
			for (int j = i + 1; j < m; ++j) {
				ident[i * m + j] -= ident[k * m + j] * coeff;
				_A[i * m + j] -= _A[k * m + j] * coeff;
			}
		}
	}
	double eLast = 1.0 / _A[(m - 1) * m + (m - 1)];
	ident[(m - 1) * m + (m - 1)] *= eLast;
	res = matrixMult(_b, ident, n, m, m);
	delete[] ident;
	return res;
}

double* matrixMult(double* _source1, double* _source2, int m, int n, int s) {
	double* result = new double[m * s];
	for (int i = 0; i < m * s; ++i) {
			result[i] = 0.0;
	}
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < s; j++) {
				result[i * s + j] += _source1[i * n + k] * _source2[k * s + j];
			}
		}
	}
	return result;
}

double* getCopy(double* _source, int m, int n) {
	double* res = new double[m * n];
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			res[i * n + j] = _source[i * n + j];
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

//void deletePointMatr(double* _source, int m) {
//	delete[] _source;
//}

double norm(double* _source, int m, int n) {
	double max = 0.0;
	for (int i = 0; i < m; ++i) {
		double sum = 0.0;
		for (int j = 0; j < n; ++j) {
			sum += fabs(_source[i * n + j]);
		}
		if (sum > max) {
			max = sum;
		}
	}
	return max;
}

double* matrixDiff(double* _source1, double* _source2, int m, int n) {
	double* result = new double[m * n];
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			result[i * n + j] = _source1[i * n + j] - _source2[i * n + j];
		}
	}
	return result;
}

double* getL(double* _source, int m) {
	double* result = new double[m * m];
	for (int i = 0; i < m * m; ++i) {
		result[i] = 0.0;
	}
	for (int i = 0; i < m; ++i) {
		result[i * m + i] = 1.0;
		for (int j = 0; j < i; ++j) {
			result[i * m + j] = _source[i * m + j];
		}
	}
	return result;
}

double* getU(double* _source, int m) {
	double* result = new double[m * m];
	for (int i = 0; i < m * m; ++i) {
		result[i] = 0.0;
	}
	for (int i = 0; i < m; ++i) {
		for (int j = i; j < m; ++j) {
			result[i * m + j] = _source[i * m + j];
		}
	}
	return result;
}

double checkLU(double* _initial, double* _final, int m) {
	double* L = getL(_final, m);
	double* U = getU(_final, m);
	double* LU = matrixMult(L, U, m, m, m);
	double* diff = matrixDiff(_initial, LU, m, m);
	delete[] L;
	delete[] U;
	delete[] LU;
	double result = norm(diff, m, m);
	delete[] diff;
	return result;
}
