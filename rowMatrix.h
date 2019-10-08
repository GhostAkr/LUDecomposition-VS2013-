#pragma once

#include <iostream>
using namespace std;

class rowMatrix
{
private:
	double* data;
	int rows;
	int cols;
public:
	rowMatrix();
	rowMatrix(int _newRows, int _newCols);
	double operator()(int _indexRow, int _indexCol);
	void print();
	~rowMatrix();
};

