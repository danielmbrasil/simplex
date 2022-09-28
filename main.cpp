#include <iostream>
#include <vector>
#include "Matrix.h"

#define MAX_ITERATIONS 10 //18446744073709551615

using ulooooooooong = unsigned long long;

template<typename T>
void minimum(std::vector<double> const &vec, T *min, unsigned *index) {
    unsigned int c_hat_size = vec.size();

    *min = vec[0];
    *index = 0;
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

Matrix createBasicCoeficientsVector(std::vector<double> const &c, std::vector<unsigned> const &B) {
    unsigned n = B.size();

    Matrix cB(n, 1);
 
    for (unsigned i = 0; i < n; i++)
        cB.setNewValueAtSpecificPositionOnMatrix(i, 0, c[B[i]]);

    return cB;
}

int main(int argc, char const *argv[])
{
    unsigned contraints_number, variables_number; // num equações, num var originais, num v folgas
    std::cin >> contraints_number >> variables_number;
    unsigned total_variables = variables_number + contraints_number;

    Matrix *A = new Matrix(contraints_number, total_variables); // matriz original
    Matrix *b = new Matrix(contraints_number, 1); // termos independentes, matriz contraints_numberx1
    std::vector<double> c(total_variables, 0.f);

    // le vetor de coeficientes de Z
    for (unsigned i = 0; i < variables_number; i++)
        std::cin >> c[i];

    // le matriz original
    A->readMatrix();

    // lee la matriz de los tiermos independientes
    b->readMatrix();

    std::vector<unsigned> basicIndexes(contraints_number, 0); // vetor dos indices basicos
    std::vector<unsigned> nonBasicIndexes(total_variables-contraints_number, 0); // vetor dos indices não básicos

    // read basic indexes vector
    for (int i = 0; i < contraints_number; i++)
        std::cin >> basicIndexes[i];

    // read non-basic indexes vector
    for (int i = 0; i < variables_number; i++)
        std::cin >> nonBasicIndexes[i];

    Matrix *B = new Matrix(A, basicIndexes); // matrix Basica
    Matrix *N = new Matrix(A, nonBasicIndexes); //non-basic matrix

    ulooooooooong it = 0;

    while (it++ < MAX_ITERATIONS) {
        std::cout << "\nITERAÇÃO " << it << std::endl;

        /* ATEÇÃO - ATENCIÓN - BEWARE 
        * AQUI COMEÇA O SIMPLEX - AQUÍ EMPIEZA EL SIMPLEX - HERE THE SIMPLEX STARTS
        * PASSO 1: TIO(Xb) <- B^(-1)*b ----> tio(xb) é um vetor das sol. iniciais basicas
        *          TIO(Xn) <- 0        ----> tio(Xn) é um vetor inicialmente nulo das soluções n basicas
        */

        // print basic matrix for testing purposes
        std::cout << "Matriz basica\n";
        B->printMatrix();

        std::cout << "Matriz nao basica\n";
        N->printMatrix();

        std::cout << "basic indexes\n";
        for (auto i : basicIndexes)
            std::cout << i << " ";
        std::cout << std::endl;

        std::cout << "non-basic indexes\n";
        for (auto i : nonBasicIndexes)
            std::cout << i << " ";
        std::cout << std::endl;

        Matrix *B_inverse = new Matrix(B->inverse()); // B^(-1)

        // printa matriz B inversa para teste
        std::cout << "Matriz B inversa\n";
        B_inverse->printMatrix();

        Matrix *Xb = new Matrix(B_inverse->multiply(b)); // basic solutions --- TIO(Xb) <- B^(-1)*b
        std::vector<double> Xn(total_variables-contraints_number, 0.f); // non-basic solutions --- TIO(Xn) <- 0 

        std::cout << "Xb" << std::endl;
        Xb->printMatrix();

        /* ATENÇÃO - ATENCIÓN - BEWARE
        * PASSO 2: Cálculo dos custos relativos
        *          2.1: vetor multiplicador simplex
        *               λ ← cB B^(−1) ---> λ vetor mult, cB custos básicos, B_inverse
        */

        Matrix *cB = new Matrix(createBasicCoeficientsVector(c, basicIndexes));

        Matrix *lambida = new Matrix(cB->transpose().multiply(B_inverse));

        std::cout << "cb\n";
        cB->printMatrix();

        std::cout << "lambida\n";
        lambida->printMatrix();

        /* ATENÇÃO - ATENCIÓN - BEWARE
        * PASSO 2: Cálculo dos custos relativos
        *          2.2: custos relativos
        *               cChapeuNj <- cNj - lambida*aNj j = 1, 2, ..., n-m
        */

       std::vector<double> c_hat;
        
        for (unsigned j = 0; j < variables_number; j++) {
            Matrix *a = new Matrix(A->getColumnMatrix(nonBasicIndexes[j]));
            Matrix *lambida_nao_basico = new Matrix(lambida->multiply(a));
            c_hat.push_back(c[nonBasicIndexes[j]] - lambida_nao_basico->getMatrix()[0][0]);

            delete a;
            delete lambida_nao_basico;
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
        std::cout << "min " << min << " k " << k << std::endl;
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

        Matrix *temp = new Matrix(N->getColumnMatrix(k));
        Matrix *y = new Matrix(B_inverse->multiply(temp));
        delete temp;

        /* ATENÇÃO - ATENCIÓN - BEWARE
        * PASO 5: determinación del paso y variable que va salir de la matriz base
        *         I) si y <= 0, entonces el problema no tiene solución finita óptima, BREAK
        *         II) de lo contrario:
        *               sombrero del evandro <- xSombreroBl/yl = min {xSombreroBi / yi>0, i=1,..., m}
        *               xBl sale de la base
        */

        if (y->isNegative()) exit(4);

        std::vector<double> temp_min;

        for (unsigned i = 0; i < contraints_number; i++) {
            if (y->getMatrix()[i][0] > 0) {
                temp_min.push_back(Xb->getMatrix()[i][0] / y->getMatrix()[i][0]);
            } else {
                temp_min.push_back(0.f); // pegar o indice l correto
            }
        }

        int l = -1;
        double epsilon_hat;
        epsilon(temp_min, &epsilon_hat, &l);

        if (l == -1) break; // y <= 0

        std::cout << "y\n";
        y->printMatrix();
        std::cout << "epsilson " << epsilon_hat << " l=" << l << std::endl;

        /* ATENÇÃO - ATENCIÓN - BEWARE
        * PASSO 6: atualização, nova partição básica
        *           troque a l-ésima coluna de B pela k-ésima coluna de N 
        */

        B->swapColumns(N, k, (unsigned)l);

        // update basic and non-basic indees
        auto temp0 = basicIndexes[l];
        basicIndexes[l] = nonBasicIndexes[k];
        nonBasicIndexes[k] = temp0;

        delete B_inverse;
        delete Xb;
        delete cB;
        delete lambida;
        delete y;
    }

    std::cout << "final\n";

    // free memory
    delete A;
    delete B;
    delete N;
    delete b;
    return 0;
}
