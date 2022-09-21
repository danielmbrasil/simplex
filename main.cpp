#include <iostream>
#include <vector>
#include "Matrix.h"

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

    Matrix A(nE, nVO+nVF); // matriz original
    Matrix b(nE, 1); // termos independentes, matriz nEx1
    Matrix c(nVO+nVF, 1); // coeficients from Z, matriz (nVO+nVF)x1
    std::vector<double> c_hat((nVO+nVF)-nE, 0.f);

    // le vetor de coeficientes de Z
    c.readMatrix();

    // le matriz original
    A.readMatrix();

    // lee el arreglo de los tiermos independientes
    b.readMatrix();

    std::vector<int> basicIndexes(nE, 0); // vetor dos indices basicos
    std::vector<int> nonBasicIndexes((nVO+nVF)-nE, 0); // vetor dos indices não básicos

    // read basic indexes vector
    for (int i = 0; i < nE; i++)
        std::cin >> basicIndexes[i];

    // read non-basic indexes vector
    for (int i = 0; i < ((nVO+nVF)-nE); i++)
        std::cin >> nonBasicIndexes[i];

    Matrix B(A, basicIndexes); // matrix Basica
    Matrix N(A, nonBasicIndexes); //non-basic matrix

    // print basic matrix for testing purposes
    std::cout << "Matriz basica\n";
    B.printMatrix();

    Matrix B_inverse = B.inverse(); // B^(-1)

    // printa matriz B inversa para teste
    B_inverse.printMatrix();

    /* ATEÇÃO - ATENCIÓN - BEWARE 
    * AQUI COMEÇA O SIMPLEX - AQUÍ EMPIEZA EL SIMPLEX - HERE THE SIMPLEX STARTS
    * PASSO 1: TIO(Xb) <- B^(-1)*b ----> tio(xb) é um vetor das sol. iniciais basicas
    *          TIO(Xn) <- 0        ----> tio(Xn) é um vetor inicialmente nulo das soluções n basicas
    */

    Matrix Xb = B_inverse * b; // basic solutions --- TIO(Xb) <- B^(-1)*b
    std::vector<double> Xn((nVO+nVF)-nE, 0.f); // non-basic solutions --- TIO(Xn) <- 0 

    /* ATENÇÃO - ATENCIÓN - BEWARE
    * PASSO 2: Cálculo dos custos relativos
    *          2.1: vetor multiplicador simplex
    *               λ ← cB B^(−1) ---> λ vetor mult, cB custos básicos, B_inverse
    */

    Matrix cB = c.getBasicCoeficientsMatrix(basicIndexes);

    Matrix lambida = cB.transpose() * B_inverse;
    Matrix lambida_t = lambida.transpose();

    std::cout << "cb\n";
    cB.printMatrix();

    std::cout << std::endl << "dimensoes lambida " << lambida.getRows() << " " << lambida.getCols() << std::endl;
    std::cout << "dimensoes lambida transposta " << lambida_t.getRows() << " " << lambida_t.getCols() << std::endl;
    std::cout << "lambida transposta\n";
    lambida_t.printMatrix();

    /* ATENÇÃO - ATENCIÓN - BEWARE
    * PASSO 2: Cálculo dos custos relativos
    *          2.2: custos relativos
    *               cChapeuNj <- cNj - lambida_t*aNj j = 1, 2, ..., n-m
    */
    
    for (unsigned j = 0; j < (nVO+nVF)-nE; j++) {
        Matrix a = N.getColumnMatrix(nonBasicIndexes[j]);
        std::cout << "dimensoes a " << a.getRows() << " " << a.getCols() << std::endl;
        Matrix lambida_nao_basico = lambida_t * a.transpose();
        c_hat[nonBasicIndexes[j]] = c.getMatrix()[nonBasicIndexes[j]][0] - lambida_nao_basico.getMatrix()[0][0]; 
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

    return 0;
}
