#ifndef PROJECT_INCLUDE_LUDECOMPOSITION_H_
#define PROJECT_INCLUDE_LUDECOMPOSITION_H_

// Main methods
double** LUDecomposition(double** _A, int m, int n);
double** LUDecompositionParal(double** _A, int m, int n);
double** LUBlockDecomposition(double** _A, int n);

// Other methods
void matrixPrint(double** _source, int m, int n);
double** createRandomMatrix(int m, int n);
bool compareMatrices(double** _source1, int m1, int n1, double** _source2, int m2, int n2);
void linSolveDown(double** _A, double** _b, int n, int m);
void linSolveUp(double** _A, double** _b, int n, int m);
double** getU22(double** _source, int m, int n);
double** getL(double** _source, int m, int n);
double** getL22(double** _source, int m, int n);
double** getL32(double** _source, int m, int n);

#endif  // PROJECT_INCLUDE_LUDECOMPOSITION_H_