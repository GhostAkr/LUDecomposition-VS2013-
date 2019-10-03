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
		//cout << "Number of threads is " << omp_get_num_threads() << endl;
	}
	return _A;
}

double** LUBlockDecomposition(double** _A, int n) {
	int b = 2;
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
	double** bl1 = new double*[b];
	for (int p = 0; p < b; ++p) {
		bl1[p] = new double[b];
	}
	double** bl2 = new double*[b];
	for (int p = 0; p < b; ++p) {
		bl2[p] = new double[b];
	}
	for (int i = 0; i < n - 1 - b; i += b) {
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
		/*if (i == b) {
			cout << "Blockkkkkk" << endl;
			partPrint(block);
		}*/
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
		LUDecomposition(block, b, b);  // TODO: make it void
		/*cout << "block is" << endl;
		matrixPrint(block, b, b);*/
		for (int k = i; k < i + b; ++k) {
			for (int l = i; l < i + b; ++l) {
				res[k][l] = block[k - i][l - i];
			}
		}



		//partPrint(res);
		block1 = linSolveDown(block, block1, b, length);// сначала считам U
		//if (i == 2) {
		//	cout << "block1 = " << block1[0][3] << endl;
		//	//
		//}
		//cout << "block1 = " << block1[0][5] << endl;
		block2 = linSolveUp(block, block2, length, b);  // потом L. матрица block портится 
		for (int k = i; k < b + i; ++k) {
			for (int l = i + b; l < n; ++l) {
				res[l][k] = block2[l - i - b][k - i];  //сначала записываем L тк на диагонали 1 и их не жалко перекрыть значениями из U
				res[k][l] = block1[k - i][l - i - b];
			}
		}


		//for (int k = i + b; k < n; ++k) {
		//	for (int l = i + b; l < b; ++l) {
		//		res[k][l] = block2[l - i - b][k - i];  //сначала записываем L тк на диагонали 1 и их не жалко перекрыть значениями из U
		//		//res[k][l] = block1[k - i][l - i - b];
		//	}
		//}
		//for (int k = i; k < b; ++k) {
		//	for (int l = i + b; l < n; ++l) {
		//		//res[l][k] = block2[l - i - b][k - i];  //сначала записываем L тк на диагонали 1 и их не жалко перекрыть значениями из U
		//		res[k][l] = block1[k - i][l - i - b];
		//	}
		//}

		//cout << "res[0][7] = " << res[0][7] << endl;
		if (i == 2) {
			cout << "res[2][7] = " << res[2][7] << endl;
		}
		double** mm = matrixMult(block2, block1, length, b, length);
		for (int k = 0; k < length; ++k) {
			for (int l = 0; l < length; ++l) {
				_A[i + b + k][i + b + l] -= mm[k][l];
			}
		}
		
		//deletePointMatr(block, b);
		_i = i;
		//_A = getCopy(res, n, n);
		//_A = res;
		if (i == n - 2 - b) {
			/*cout << "length = " << length << endl;
			cout << "block2[0][0] = " << block2[0][0] << endl;
			cout << "block2[0][1] = " << block2[0][1] << endl;
			cout << "block2[1][0] = " << block2[1][0] << endl;
			cout << "block2[1][1] = " << block2[1][1] << endl;*/
			/*cout << "Block11111" << endl;
			matrixPrint(block1, b, b);
			cout << "Block22222" << endl;
			matrixPrint(block2, b, b);*/
			/*cout << "A is" << endl;
			matrixPrint(_A, n, n);*/
			for (int k = 0; k < b; ++k) {
				for (int l = 0; l < b; ++l) {
					block[k][l] = _A[_i + b + k][_i + b + l];
				}
			}
			LUDecomposition(block, b, b);
			for (int i = n - b; i < n; ++i) {
				for (int j = n - b; j < n; ++j) {
					res[i][j] = block[i - n + b][j - n + b];
				}
			}
			bl1 = getCopy(block1, b, b);
			bl2 = getCopy(block2, b, b);
		}
		deletePointMatr(block1, b);
		deletePointMatr(block2, length);
	}
		//lastsize = n - _i;
		//double** last = new double*[b]; // для общего случая. если последняя матрица такого же размера как первая, можно использовать block
		//for (int p = 0; p < b; ++p) {
		//	last[p] = new double[b];
		//}
		//LUDecomposition(last, b, b);
		/*cout << "bl1111" << endl;
		matrixPrint(bl1, b, b);
		cout << "bl2222" << endl;
		matrixPrint(bl2, b, b);*/




	/*double** lastBlock = new double*[b];
	for (int i = 0; i < b; ++i) {
		lastBlock[i] = new double[b];
	}
		double** mm = matrixMult(bl2, bl1, length, b, length);
		for (int k = 0; k < b; ++k) {
			for (int l = 0; l < b; ++l) {
				lastBlock[k][l] = _A[_i + b + k][_i + b + l] - mm[k][l];
			}
		}
		cout << "Last block" << endl;
		matrixPrint(lastBlock, b, b);
		cout << endl;
		LUDecomposition(lastBlock, b, b);
		for (int i = n - b; i < n; ++i) {
			for (int j = n - b; j < n; ++j) {
				res[i][j] = lastBlock[i - n + b][i - n + b];
			}
		}*/


		/*for (int i = n - b; i < n; ++i) {
			for (int j = n - b; j < n; ++j) {
				res[i][j] = _A[i - n + b][i - n + b];
			}
		}*/





		/*for (int k = 0; k < length; ++k) {
			for (int l = 0; l < length; ++l) {
				_A[_i + b + k][_i + b + l] -= mm[k][l];
			}
		}*/
		//LUDecomposition(block, b, b);
		/*for (int k = _i; k < b; ++k) {
			for (int l = _i; l < b; ++l) {
				res[k][l] = block[k][l];
			}
		}*/
		//_A = res;
		//_A = getCopy(res, n, n);
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
	double epsNull = 1e-5;
	if (m1 != m2 || n1 != n2) {
		return false;
	}
	for (int i = 0; i < m1; ++i) {
		for (int j = 0; j < n1; ++j) {
			if (abs(_source1[i][j] - _source2[i][j]) > epsNull) {
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
	//cout << "_b[0][0] = " << _b[0][0] << endl;
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < m; ++j){
			res[i][j] = _b[i][j];
			for (int k = 0; k < i; k++) {
				//cout << "_A[i][k] = " << _A[i][k] << endl;
				//cout << "res[k][j] = " << res[k][j] << endl;
				res[i][j] -= _A[i][k] * res[k][j];
			}
		}
	}
	//cout << "res = " << res[1][0] << endl;
	return res;
	//_b = res;
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
	//double e;
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
			//cout << ident[i][j] << endl;
		}
		//cout << "_A[0][2] = " << _A[0][2] << endl;
		for (int k = i + 1; k < m; ++k) {
			double coeff = 1.0 / _A[k][k];
			coeff *= _A[i][k];
			for (int j = i + 1; j < m; ++j) {
				
				ident[i][j] -= ident[k][j] * coeff;
				_A[i][j] -= _A[k][j] * coeff;
			}
			//cout << "_A[0][2] = " << _A[0][2] << endl;
			/*cout << endl;
			matrixPrint(_A, m, m);
			cout << endl;*/
		}
	}
	//_A[m - 1][m - 1] /= _A[m - 1][m - 1];
	double eLast = 1.0 / _A[m - 1][m - 1];
	ident[m - 1][m - 1] *= eLast;
	res = matrixMult(_b, ident, n, m, m);
	/*cout << "QQQB" << endl;
	partPrint(res);*/
	/*cout << "Inv" << endl;
	matrixPrint(ident, m, m);*/
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
	cout << "source[0][0] = " << _source[0][0] << endl;
	cout << "source[0][1] = " << _source[0][1] << endl;
	cout << "source[1][0] = " << _source[1][0] << endl;
	cout << "source[1][1] = " << _source[1][1] << endl;
	cout << "m = " << m << endl;
	cout << "n = " << n << endl;
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
	for (int i = 0; i < 20; ++i) {
		for (int j = 0; j < 20; ++j) {
			cout << _source[i][j] << " ";
		}
		cout << endl;
	}
}

void deletePointMatr(double** _source, int m) {
	for (int i = 0; i < m; ++i) {
		delete[] _source[i];
	}
	delete[] _source;
}
