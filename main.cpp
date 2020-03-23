#include <iostream>
#include <cmath>

#include "headers/MatrixMaster.h"

using namespace std;

const double eps = pow(10, -6);

void copy(int n, double** arr1, double** arr2)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            arr2[i][j] = arr1[i][j];
        }
    }
}

void runningMethod(int n, double** matrix, double* f)
{
    MatrixMaster::printMatrix(n, n, matrix, f);
    cout << endl;

    cout << "Determinant: " << MatrixMaster::determinant(n, matrix) << endl;
    cout << "Conditionality: " << MatrixMaster::conditionality(n, matrix) << endl << endl;

    try
    {
        auto result = MatrixMaster::tridiagonalSolve(n, matrix, f);
        for(int i = 0; i < n; i++) cout << "x" << i << " = " << result[i] << ", ";
        cout << endl << endl;

        //перевірка знайдених коренів
        MatrixMaster::printMatrixWithValues(n, n, matrix, result);
    }
    catch (int e)
    {
        cout << "Для даної матриці метод Прогонки застосовувати не можна(ділення на 0)!" << endl;
    }
}

double* seidelMethod(int n, double** matrix_, double* f, double eps)
{
    auto matrix = new double*[n];
    for(int i = 0; i < n; i++) matrix[i] = new double[n];
    copy(n, matrix_, matrix);
    MatrixMaster::printMatrix(n, n, matrix, f);
    cout << endl;

    //перевірка умови збіжності
    for(int i = 0; i < n; i++)
    {
        double sum = 0;
        for(int j = 0; j < n; j++)
        {
            if(i != j) sum += abs(matrix_[i][j]);
        }

        if(sum > abs(matrix_[i][i])){
            cout << "Для даної матриці метод Зейделя не сходиться!" << endl;
            throw -2;
        }
    }

    cout << "Determinant: " << MatrixMaster::determinant(n, matrix) << endl;
    cout << "Conditionality: " << MatrixMaster::conditionality(n, matrix) << endl << endl;

    //підготовка вхідних данних
    auto result = new double[n];
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if(i != j) matrix[i][j] = matrix[i][j] / matrix[i][i] * -1;
        }
        f[i] = f[i] / matrix[i][i];
        matrix[i][i] = 0;
        result[i] = f[i];
    }

    //ітераційний процес
    auto deltas = new double[n];
    for (int iter = 0; iter >= 0; iter++)
    {
        cout << "Iteration: " << iter << endl;
        for(int i = 0; i < n; i++)
        {
            cout << "x" << i << " = " << result[i] << ", ";
        }
        cout << endl;

        for(int i = 0; i < n; i++)
        {
            double newX = f[i];

            for(int j = 0; j < n; j++)
            {
                newX += matrix[i][j] * result[j];
            }

            deltas[i] = abs(result[i] - newX);
            result[i] = newX;
        }

        double max = MatrixMaster::max(n, deltas);
        if(max <= eps || iter >= 1000) break;
    }

    //перевірка знайденик коренів
    MatrixMaster::printMatrixWithValues(n, n, matrix_, result);

    return result;
}

int main() {
    int n = 3, m = 3;

    auto matrix = new double*[n]{
        new double[m]{ 5, 2, 0},
        new double[m]{ 2, 8, 3},
        new double[m]{ 0, 2, 2}
    };
    auto f = new double[n]{-1, 13, 15};

    cout << "Метод 'Прогонки'" << endl;
    runningMethod(n, matrix, f);

    n = 3, m = 3;
    matrix = new double*[n]{
            new double[m]{ 12,  4, -7},
            new double[m]{ -2,  8, -5},
            new double[m]{ -1,  6, -8}
    };
    f = new double[n]{12, 13, 14};

    cout << endl << "Метод Зейделя" << endl;
    seidelMethod(n, matrix, f, eps);

    return 0;
}