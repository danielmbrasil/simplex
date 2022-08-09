#include <iostream>
#include <vector>
#include <cmath>

std::vector<std::vector<int>> multiply(std::vector<std::vector<int>> &A, std::vector<std::vector<int>> &B)
{
    int columnsA = A[0].size(); // columns A
    int rowsB = B.size(); // rows B

    if (columnsA != rowsB) exit(1); // check if it's possible to multiply

    // new matrix dimensions
    int rowsA = A.size(); // rows A
    columnsA = A[0].size(); // columns A
    int columnsB = B[0].size(); // columns B

    std::vector<std::vector<int>> C(rowsA, std::vector<int>(columnsB, 0));

    for (int i = 0; i < rowsA; ++i)
        for (int j = 0; j < columnsB; ++j)
            for (int k = 0; k < columnsA; ++k)
                C[i][j] += A[i][k] * B[k][j];

    return C;
}

long long determinant(std::vector<std::vector<int>> &A) {
    int rows = A.size();
    int columns = A[0].size();

    if (rows != columns) exit(1); // non-square matrix

    if (columns == 1) return A[0][0];

    // base case, 2x2 matrix
    if (columns == 2)
        return A[0][0]*A[1][1] - A[0][1]*A[1][0];

    std::vector<std::vector<int>> subMatrix(rows-1, std::vector<int>(columns-1, 0));
    long long det = 0;

    for (int x = 0; x < columns; x++) {
        int subi = 0;
        for (int i = 1; i < columns; i++) {
            int subj = 0;
            for (int j = 0; j < columns; j++) {
                if (x == j) continue;

                subMatrix[subi][subj] = A[i][j];
                subj++;
            }
            subi++;
        }
        det += (pow(-1, x) * A[0][x] * determinant(subMatrix));
    }

    return det;
}

std::vector<std::vector<int>> transpose(std::vector<std::vector<int>> &A) {
    int rows = A.size();
    int columns = A[0].size();

    std::vector<std::vector<int>> B(columns, std::vector<int>(rows, 0));

    for (int i = 0; i < columns; i++)
        for (int j = 0; j < rows; j++)
            B[i][j] = A[j][i];

    return B;
}

int main(int argc, char const *argv[])
{
    int n1, m1;
    std::cin >> n1 >> m1;
    std::vector<std::vector<int>> A(n1, std::vector<int>(m1, 0));

    for (int i = 0; i < n1; i++)
        for (int j = 0; j < m1; j++)
            std::cin >> A[i][j];
/*
    int n2, m2;
    std::cin >> n2 >> m2;
    std::vector<std::vector<int>> B(n2, std::vector<int>(m2, 0));
    for (int i = 0; i < n2; i++)
        for (int j = 0; j < m2; j++)
            std::cin >> B[i][j];

    std::vector<std::vector<int>> C = multiply(A, B);

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < m2; j++)
            std::cout << C[i][j] << " ";
        std::cout << std::endl;
    }*/

    //std::cout << "determinant " << determinant(A) << std::endl;

    auto B = transpose(A);

    for (int i = 0; i < m1; i++) {
        for (int j = 0; j < n1; j++)
            std::cout <<  B[i][j] << " ";
        std::cout << std::endl;
    }
    return 0;
}
