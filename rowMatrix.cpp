#include "rowMatrix.h"


rowMatrix::rowMatrix() {
	rows = 0;
	cols = 0;
	data = nullptr;
}


rowMatrix::~rowMatrix() {
	delete[] data;
}

rowMatrix::rowMatrix(int _newRows, int _newCols) {
	rows = _newRows;
	cols = _newCols;
	data = new double[rows * cols];
	for (int i = 0; i < rows * cols; ++i) {
		data[i] = 0.0;
	}
}

double rowMatrix::operator()(int _indexRow, int _indexCol) {
	return data[_indexRow * cols + _indexCol];
}

void rowMatrix::print() {
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			cout << data[i * cols + j] << endl;
		}
		cout << endl;
	}
}
