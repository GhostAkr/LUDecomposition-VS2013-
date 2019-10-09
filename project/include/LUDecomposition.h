#ifndef PROJECT_INCLUDE_LUDECOMPOSITION_H_
#define PROJECT_INCLUDE_LUDECOMPOSITION_H_

// Main methods
double* LUDecomposition(double* _A, int m, int n);
double* LUDecompositionParal(double* _A, int m, int n);
double* LUBlockDecomposition(double* _A, int n, int b);
double* LUBlockDecompositionParal(double* _A, int n, int b);

// Other methods
double norm(double* _source, int m, int n);
void matrixPrint(double** _source, int m, int n);
double** createRandomMatrix(int m, int n);
double* createRandomRowMatrix(int m, int n);
bool compareMatrices(double** _source1, int m1, int n1, double** _source2, int m2, int n2);
double* linSolveDown(double* _A, double* _b, int n, int m);
double* linSolveUp(double* _A, double* _b, int n, int m);
double* linSolveDownParal(double* _A, double* _b, int n, int m);
double* linSolveUpParal(double* _A, double* _b, int n, int m);
double* getL(double* _source, int m);
double* getU(double* _source, int m);
double* matrixMult(double* _source1, double* _source2, int m, int n, int s);  // m, n --- size of _source1; s --- nOfCols in _source2
double* matrixDiff(double* _source1, double* _source2, int m, int n);
double* getCopy(double* _source, int m, int n);
double checkLU(double* _initial, double* _final, int m);


void deletePointMatr(double* _source, int m);

void partPrint(double** _source);

#endif  // PROJECT_INCLUDE_LUDECOMPOSITION_H_