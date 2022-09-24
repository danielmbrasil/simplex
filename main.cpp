#include <iostream>
#include <vector>
#include "Matrix.h"

#define MAX_ITERATIONS 10 //18446744073709551615

using ulooooooooong = unsigned long long;

template<typename T>
void minimum(std::vector<double> const &vec, T *min, unsigned *index) {
    unsigned int c_hat_size = vec.size();

    *min = vec[0];
    for (int i = 1; i < c_hat_size; i++) {
        if (*min > vec[i])
            *min = vec[i], *index = i;
    }
}

template<typename T>
void epsilon(std::vector<double> const &vec, T *min, int *index) {
    unsigned vec_size = vec.size();

    *min = 4294967295;
    for (int i = 0; i < vec_size; i++) {
        if (vec[i] > 0 && *min > vec[i])
            *min = vec[i], *index = i;
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

    ulooooooooong it = 0;

    while (it++ < MAX_ITERATIONS) {
        /* ATEÇÃO - ATENCIÓN - BEWARE 
        * AQUI COMEÇA O SIMPLEX - AQUÍ EMPIEZA EL SIMPLEX - HERE THE SIMPLEX STARTS
        * PASSO 1: TIO(Xb) <- B^(-1)*b ----> tio(xb) é um vetor das sol. iniciais basicas
        *          TIO(Xn) <- 0        ----> tio(Xn) é um vetor inicialmente nulo das soluções n basicas
        */

       Matrix B_inverse = B.inverse(); // B^(-1)

        // printa matriz B inversa para teste
        std::cout << "Matriz B inversa\n";
        B_inverse.printMatrix();

        std::cout << "dim " << b.getRows() << " " << b.getCols() << std::endl;
        std::cout << "  aa" << b.test() << std::endl;
        std::cout << std::endl;
        Matrix Xb = B_inverse * b; // basic solutions --- TIO(Xb) <- B^(-1)*b
        std::cout << "n acredito\n";
        std::vector<double> Xn((nVO+nVF)-nE, 0.f); // non-basic solutions --- TIO(Xn) <- 0 

        /* ATENÇÃO - ATENCIÓN - BEWARE
        * PASSO 2: Cálculo dos custos relativos
        *          2.1: vetor multiplicador simplex
        *               λ ← cB B^(−1) ---> λ vetor mult, cB custos básicos, B_inverse
        */
        std::cout << "ola\n";
        Matrix cB = c.getBasicCoeficientsMatrix(basicIndexes);
        std::cout << "ola tudo bem\n";
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
            Matrix lambida_nao_basico = lambida_t * a.transpose();
            c_hat[nonBasicIndexes[j]] = c.getMatrix()[nonBasicIndexes[j]][0] - lambida_nao_basico.getMatrix()[0][0]; 
        }

        std::cout << "c hat\n";

        for (auto c : c_hat) {
            std::cout << c << " ";
        }
        std::cout << std::endl;

        /* ATENÇÃO - ATENCIÓN - BEWARE
        * PASO 2: cálculo de los costos relativos
        *         determinación de la variable a entrar en la matriz básica
        *         cSombreroNk <- min {cSombreroNj, j=1,2,..,n-m} (la variable xNk entra en la matriz básica)
        *         k es el valor del índice j del valor minimo de cSombrero 
        */

        int min;
        unsigned k;
        minimum(c_hat, &min, &k);
        c_hat[nonBasicIndexes[k]] = min;

        /* ATENÇÃO - ATENCIÓN - BEWARE
        *  3rd STEP: optimality test
                    if c_hatNk >= 0: STOP (the solution at the currect iteration is optimal) 🙏
        */

       if (c_hat[nonBasicIndexes[k]] >= 0) break;


       /* ATENÇÃO - ATENCIÓN - BEWARE
       * PASSO 4: Cálculo da direção simplex
       *          y <- B^(-1) * aNk (equivalente a resolva o sistema By=aNk)
       */

        Matrix y = B_inverse * N.getColumnMatrix(k);

        /* ATENÇÃO - ATENCIÓN - BEWARE
        * PASO 5: determinación del paso y variable que va salir de la matriz base
        *         I) si y <= 0, entonces el problema no tiene solución finita óptima, BREAK
        *         II) de lo contrario:
        *               sombrero del evandro <- xSombreroBl/yl = min {xSombreroBi / yi>0, i=1,..., m}
        *               xBl sale de la base
        */

        if (y.isNegative()) exit(4);

        std::vector<double> temp_min;

        for (unsigned i = 0; i < nE; i++) {
            if (y.getMatrix()[i][0] > 0) {
                temp_min.push_back(Xb.getMatrix()[i][0] / y.getMatrix()[i][0]);
            } else {
                temp_min.push_back(0.f); // pegar o indice l correto
            }
        }

        int l = -1;
        double epsilon_hat;
        epsilon(temp_min, &epsilon_hat, &l);

        if (l == -1) break; // y <= 0

        std::cout << "chapeu evandro " << epsilon_hat << " l " << l << std::endl;

        /* ATENÇÃO - ATENCIÓN - BEWARE
        * PASSO 6: atualização, nova partição básica
        *           troque a l-ésima coluna de B pela k-ésima coluna de N 
        */

        B.swapColumns(N, k, (unsigned)l);

        // update basic and non-basic indees
        auto temp = basicIndexes[l];
        basicIndexes[l] = nonBasicIndexes[k];
        nonBasicIndexes[k] = temp;
    }

    std::cout << "final\n";

    for (auto i : c_hat)
        std::cout << i << " ";
    std:: cout << std::endl;

    return 0;
}
