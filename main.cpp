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
    int columnsB = B[0].size(); // columns B

    std::vector<std::vector<D>> C(rowsA, std::vector<D>(columnsB, 0));

    for (int i = 0; i < rowsA; ++i)
        for (int j = 0; j < columnsB; ++j)
            for (int k = 0; k < columnsA; ++k)
                C[i][j] += A[i][k] * B[k][j];

    return C;
}


std::vector<double> multiplyMatrixByVector(std::vector<std::vector<double>> const &A, std::vector<double> const &b) {
    int n = A.size(); // A is a square matrix

    if (n != b.size()) exit(25);

    std::vector<double> B(n, 0.f);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            B[i] += A[i][j]*b[j];

    return B;
}

double determinant(std::vector<std::vector<double>> &A) {
    int rows = A.size();
    int columns = A[0].size();

    if (rows != columns) exit(1); // non-square matrix

    if (columns == 1) return A[0][0];

    // base case, 2x2 matrix
    if (columns == 2)
        return A[0][0]*A[1][1] - A[0][1]*A[1][0];

    std::vector<std::vector<double>> subMatrix(rows-1, std::vector<double>(columns-1, 0.f));
    double det = 0;

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

std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> &A) {
    int rows = A.size();
    int columns = A[0].size();

    std::vector<std::vector<double>> B(columns, std::vector<double>(rows, 0.f));

    for (int i = 0; i < columns; i++)
        for (int j = 0; j < rows; j++)
            B[i][j] = A[j][i];

    return B;
}

std::vector<std::vector<double>> getCofactor(std::vector<std::vector<double>> &A, int x, int y) {
    int n = A.size();
    std::vector<std::vector<double>> B(n-1, std::vector<double>(n-1, 0.f));

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

std::vector<std::vector<double>> adjoint(std::vector<std::vector<double>> &A) {
    int n = A.size();
    std::vector<std::vector<double>> adj(n, std::vector<double>(n, 0.f));

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

std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> &A) {
    int det = determinant(A);
    std::cout << "det " << det << std::endl;
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

std::vector<std::vector<double>> createNewMatrix(std::vector<std::vector<double>> const &A, std::vector<int> const &indexes) {
    int rows = A.size();
    int cols = indexes.size();

    std::vector<std::vector<double>> newMatrix(rows, std::vector<double>(cols, 0.f));

    int j = 0;
    for (auto col : indexes) {
        for (int i = 0; i < rows; i++) {
            newMatrix[i][j] = A[i][col];
        }
        j++;
    }

    return newMatrix;
}

std::vector<std::vector<double>> createBasicCoeficientsVector(std::vector<std::vector<double>> const &c, std::vector<int> const &B) {
    int n = B.size();

    std::vector<std::vector<double>> cB(n, std::vector<double>(1, 0.f));
 
    for (int i = 0; i < n; i++)
        cB[i][0] = c[B[i]][0];

    return cB;
}

std::vector<std::vector<double>> getColumnFromMatrix(std::vector<std::vector<double>> const &A, int col) {
    int rows = A.size();

    std::vector<std::vector<double>> columnMatrix(rows, std::vector<double>(1, 0.f));

    for (int i = 0; i < rows; i++)
        columnMatrix[i][0] = A[i][col];

    return columnMatrix;
}

void minimum(std::vector<double> const &c_hat, int *min, int *k) {
    unsigned int c_hat_size = c_hat.size();

    *min = c_hat[0];
    for (int i = 1; i < c_hat_size; i++) {
        if (*min > c_hat[i])
            *min = c_hat[i], *k = i;
    }
}

int main(int argc, char const *argv[])
{
    int nE, nVO, nVF; // num equações, num var originais, num v folgas
    std::cin >> nE >> nVO >> nVF;

    std::vector<std::vector<double>> A(nE, std::vector<double>(nVO+nVF, 0.f)); // matriz original
    std::vector<std::vector<double>> b(nE, std::vector<double>(1, 0.f)); // termos independentes, matriz nEx1
    std::vector<std::vector<double>> c(nVO+nVF, std::vector<double>(1, 0.f)); // coeficients from Z, matriz (nVO+nVF)x1
    std::vector<double> c_hat((nVO+nVF)-nE, 0.f);

    // le vetor de coeficientes de Z
    for (int i = 0; i < nVO; i++)
        std::cin >> c[i][0];

    // le matriz original
    for (int i = 0; i < nE; i++)
        for (int j = 0; j < (nVO+nVF); j++)
            std::cin >> A[i][j];

    // lee el arreglo de los tiermos independientes
    for (int i = 0; i < nE; i++)
        std::cin >> b[i][0];

    std::vector<int> basicIndexes(nE, 0); // vetor dos indices basicos
    std::vector<int> nonBasicIndexes((nVO+nVF)-nE, 0); // vetor dos indices não básicos

    // read basic indexes vector
    for (int i = 0; i < nE; i++)
        std::cin >> basicIndexes[i];

    // read non-basic indexes vector
    for (int i = 0; i < ((nVO+nVF)-nE); i++)
        std::cin >> nonBasicIndexes[i];

    auto B = createNewMatrix(A, basicIndexes); // matrix Basica
    auto N = createNewMatrix(A, nonBasicIndexes); //non-basic matrix

    // print basic matrix for testing purposes
    std::cout << "Matriz basica" << std::endl;
    for (int i = 0; i < nE; i++) {
        for (int j = 0; j < nE; j++)
            std::cout << B[i][j] << " ";
        std::cout << std::endl;
    }

    auto B_inverse = inverse(B); // B^(-1)

    // printa matriz B inversa para teste
    std::cout << "B^(-1) \n";
    for (int i = 0; i < nE; i++) {
        for (int j = 0; j < nE; j++)
            std::cout << B_inverse[i][j] << " ";
        std::cout << std::endl;
    }

    /* ATEÇÃO - ATENCIÓN - BEWARE 
    * AQUI COMEÇA O SIMPLEX - AQUÍ EMPIEZA EL SIMPLEX - HERE THE SIMPLEX STARTS
    * PASSO 1: TIO(Xb) <- B^(-1)*b ----> tio(xb) é um vetor das sol. iniciais basicas
    *          TIO(Xn) <- 0        ----> tio(Xn) é um vetor inicialmente nulo das soluções n basicas
    */


    auto Xb = multiply(B_inverse, b); // basic solutions --- TIO(Xb) <- B^(-1)*b
    std::vector<double> Xn((nVO+nVF)-nE, 0.f); // non-basic solutions --- TIO(Xn) <- 0 

    /* ATENÇÃO - ATENCIÓN - BEWARE
    * PASSO 2: Cálculo dos custos relativos
    *          2.1: vetor multiplicador simplex
    *               λ ← cB B^(−1) ---> λ vetor mult, cB custos básicos, B_inverse
    */

    auto cB = createBasicCoeficientsVector(c, basicIndexes);
    auto ctb = transpose(cB);
    auto lambida = multiply(ctb, B_inverse);
    auto lambida_t = transpose(lambida);

    std::cout << "cb\n";
    auto sla = cB.size();
    for (int i = 0; i < sla; i++)
        std::cout << cB[i][0] << " ";
    std::cout << std::endl << "lambida " << lambida.size() << " " << lambida[0].size() << std::endl;
    auto lambida_size = lambida_t.size();
    for (int i = 0; i < lambida_size; i++)
        std::cout << lambida_t[i][0] << " ";
    std::cout << std::endl;


    /* ATENÇÃO - ATENCIÓN - BEWARE
    * PASSO 2: Cálculo dos custos relativos
    *          2.2: custos relativos
    *               cChapeuNj <- cNj - lambida_t*aNj j = 1, 2, ..., n-m
    */

    std::cout << "dimensoes lambida transposta " << lambida_t.size() << " " << lambida_t[0].size() << std::endl;
    for (int j = 0; j < (nVO+nVF)-nE; j++) {
        auto a = getColumnFromMatrix(N, nonBasicIndexes[j]);
        auto a_t = transpose(a);
        std::cout << "dimensoes a " << a.size() << " " << a[0].size() << std::endl;
        auto lambida_nao_basico = multiply(lambida_t, a_t);
        c_hat[nonBasicIndexes[j]] = c[nonBasicIndexes[j]][0] - lambida_nao_basico[0][0]; 
    }

    std::cout << "c hat\n";

    for (auto c : c_hat) {
        std::cout << c << " ";
    }
    std::cout << std::endl;

    /* ATENÇÃO - ATENCIÓN - BEWARE
    * PASO 2: determinación de la variable a entrar en la matriz básica
    *         cSombreroNk <- min {cSombreroNj, j=1,2,..,n-m} (la variable xNk entra en la matriz básica)
    *         k es el valor del índice j del valor minimo de cSombrero 
    */

    int min, k;
    minimum(c_hat, &min, &k);
    c_hat[nonBasicIndexes[k]] = min;

        std::cout << "c hat\n";

    for (auto c : c_hat) {
        std::cout << c << " ";
    }
    std::cout << std::endl;

    return 0;
}
