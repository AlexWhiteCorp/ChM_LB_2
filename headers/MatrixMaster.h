#ifndef LB2_MATRIXMASTER_H
#define LB2_MATRIXMASTER_H

#endif //LB2_MATRIXMASTER_H

#include <iostream>
#include <cmath>

using namespace std;

class MatrixMaster
{
public:
    static double determinant(int n,  double** matrix);
    static double conditionality(int n,  double** matrix);
    static double** swapMatrix(int n,  double** matrix);
    static double* tridiagonalSolve (int n, double** matrix, double *f);
    static void printMatrix(int n, int m, double** matrix, double* f);
    static void printMatrixWithValues(int n, int m, double** matrix, double* f);
    static double mNorm(int n, int m, double** matrix);
    static double lNorm(int n, int m, double** matrix);
    static double sqrtNorm(int n, int m, double** matrix);
    static double max(int n,  double* matrix);
private:
    static double** minor(int n, int head_i, int head_j, double** matrix);

};