#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

template<typename T, typename D>
std::vector<std::vector<D>> multiply(std::vector<std::vector<T>> &A, std::vector<std::vector<D>> &B)
{
    int columnsA = A[0].size(); // columns A
    int rowsB = B.size(); // rows B

    if (columnsA != rowsB) exit(1); // check if it's possible to multiply

    // new matrix dimensions
    int rowsA = A.size(); // rows A
    columnsA = A[0].size(); // columns A
    int columnsB = B[0].size(); // columns B

    std::vector<std::vector<D>> C(rowsA, std::vector<D>(columnsB, 0));

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

std::vector<std::vector<int>> getCofactor(std::vector<std::vector<int>> &A, int x, int y) {
    int n = A.size();
    std::vector<std::vector<int>> B(n-1, std::vector<int>(n-1, 0));

    int i = 0, j = 0;

    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            if (row != x && y != col) {
                B[i][j++] = A[row][col];

                if (j == n-1) {
                    j = 0;
                    i++;
                }
            }
        }
    }

    return B;
}

std::vector<std::vector<int>> adjoint(std::vector<std::vector<int>> &A) {
    int n = A.size();
    std::vector<std::vector<int>> adj(n, std::vector<int>(n, 0));

    if (n==1) {
        adj[0][0] = 1;
        return adj;
    }

    int sign = 1;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            auto temp = getCofactor(A, i, j);

            // if i+j is even, sign is positive
            sign = ((i+j) % 2 == 0) ? 1 : -1;

            adj[i][j] = sign * determinant(temp);
        }
    }

    return transpose(adj);
}

std::vector<std::vector<double>> inverse(std::vector<std::vector<int>> &A) {
    int det = determinant(A);

    if (det == 0) exit(1);

    int n = A.size();

    std::vector<std::vector<double>> B(n, std::vector<double>(n, 0.f));

    /*
    A^(-1) = (1/det(A))*adj(A)
    */

   double in_det = 1.0/det;
   auto adj = adjoint(A);

   for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        B[i][j] = (double)(in_det*adj[i][j]);

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

    auto B = inverse(A);
    auto C = multiply(A, B);

    for (int i = 0; i < m1; i++) {
        for (int j = 0; j < n1; j++)
            std::cout <<  B[i][j] << " ";
        std::cout << std::endl;
    }

        for (int i = 0; i < m1; i++) {
        for (int j = 0; j < n1; j++)
            std::cout <<  C[i][j] << " ";
        std::cout << std::endl;
    }
    return 0;
}
