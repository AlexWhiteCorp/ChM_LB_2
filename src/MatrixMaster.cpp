#include "../headers/MatrixMaster.h"

double MatrixMaster::determinant(int n,  double** matrix)
{
    double det = 0;
    for(int i = 0; i <= n; i++)
    {
        if(n == 1) return matrix[0][0];

        det += pow(-1, 1 + (i + 1)) * matrix[0][i] * determinant(n - 1, minor(n, 0, i, matrix));
    }

    return det;
}

double MatrixMaster::conditionality(int n,  double** matrix)
{
    return lNorm(n, n, matrix) * lNorm(n, n, swapMatrix(n, matrix));
}

double* MatrixMaster::tridiagonalSolve (int n, double** matrix, double *f)
{
    //діагоналі
    auto top = new double[n]{0};
    auto main = new double[n]{0};
    auto bottom = new double[n]{0};
    auto result = new double[n];

    //зчитування трьох діагоналей з вхідної матриці
    top[0] = matrix[0][1];
    main[0] = matrix[0][0];
    main[n - 1] = matrix[n - 1][n - 1];
    bottom[n - 1] = matrix[n - 1][n - 2];
    bottom[n - 1] = matrix[n - 1][n - 2];

    for(int i = 1; i < n - 1; i++)
    {
        top[i] = matrix[i][i + 1];
        main[i] = matrix[i][i];
        bottom[i] = matrix[i][i - 1];
    }

    //прямий хід
    top[0] = top[0] / main[0];
    for (int i = 1; i < n - 1; i++){
        double div = (main[i] - top[i - 1] * bottom[i]);
        top[i] = top[i] / div;
        if(div == 0) throw -3;
    }

    f[0] = f[0] / main[0];//
    for (int i = 1; i < n; i++){
        double div = main[i] - top[i - 1] * bottom[i];
        f[i] = (f[i] - f[i - 1] * bottom[i]) / div;
        if(div == 0) throw -3;
    }

    //зворотній хід
    result[n - 1] = f[n - 1];
    for (int i = n - 2; i >= 0; i--)
        result[i] = f[i] - top[i] * result[i + 1];

    return result;
}

double** MatrixMaster::swapMatrix(int n,  double** matrix)
{
    double det = determinant(n, matrix);
    if(det == 0){
        cout << "Error! Swaped Matrix: det = 0!" << endl;
        throw -1;
    }

    auto swapedMatrix = new double*[n];
    for(int i = 0; i < n; i++) swapedMatrix[i] = new double[n];

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            swapedMatrix[i][j] = pow(-1, (i + 1) + (j + 1)) * determinant(n - 1, minor(n, i, j, matrix));
            swapedMatrix[i][j] /= det;
        }
    }

    return swapedMatrix;
}

double** MatrixMaster::minor(int n, int head_i, int head_j, double** matrix)
{
    auto minor = new double*[n - 1];
    for(int j = 0; j < n - 1; j++) minor[j] = new double[n - 1];

    for(int j = 0, currPos_i = 0; j < n; j++)
    {
        if(j == head_i) continue;

        for(int k = 0, currPos_j = 0; k < n; k++)
        {
            if(k != head_j){
                minor[currPos_i][currPos_j] = matrix[j][k];
                currPos_j++;
            }
        }

        currPos_i++;
    }

    return minor;
}

double MatrixMaster::mNorm(int n, int m, double** matrix)
{
    double sums[n];

    for(int i = 0; i < n; i++)
    {
        double sum = 0;
        for(int j = 0; j < m; j++)
        {
            sum += abs(matrix[i][j]);
        }
        sums[i] = sum;
    }

    return max(n, sums);
}

double MatrixMaster::lNorm(int n, int m, double** matrix)
{
    double sums[m];

    for(int i = 0; i < m; i++)
    {
        double sum = 0;
        for(int j = 0; j < n; j++)
        {
            sum += abs(matrix[j][i]);
        }
        sums[i] = sum;
    }

    return max(m, sums);
}

double MatrixMaster::sqrtNorm(int n, int m, double** matrix)
{
    double sum = 0;
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < n; j++)
        {
            sum += pow(abs(matrix[j][i]), 2);
        }
    }

    return sqrt(sum);
}

void MatrixMaster::printMatrix(int n, int m, double** matrix, double* f)
{
    for(int i = 0; i < n; i++)
    {
        cout << "|| ";
        for(int j = 0; j < n; j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << "| " << f[i] << " ||" << endl;
    }
}

void MatrixMaster::printMatrixWithValues(int n, int m, double** matrix, double* f)
{
    for(int i = 0; i < n; i++)
    {
        double f_i = 0;
        cout << matrix[i][0] << " * x" << 1;
        f_i += matrix[i][0] * f[0];
        for(int j = 1; j < n; j++)
        {
            f_i += matrix[i][j] * f[j];
            cout << " + " << matrix[i][j] << " * x" << (j + 1);
        }
        cout << " = " << f_i << endl;
    }
}

double MatrixMaster::max(int n,  double* matrix)
{
    double max = matrix[0];

    for(int i = 1; i < n; i++)
    {
        if(matrix[i] > max) max = matrix[i];
    }

    return max;
}