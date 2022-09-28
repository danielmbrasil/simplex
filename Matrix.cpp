#include "Matrix.h"
#include <iostream>
#include <cmath>

Matrix::Matrix(unsigned rows, unsigned cols) {
    this->m_rows = rows;
    this->m_cols = cols;
    this->m_matrix.resize(rows, std::vector<double>(cols, 0.f));
}

Matrix::Matrix(Matrix *A, std::vector<unsigned> const &indexes) {
    this->m_rows = A->getRows();
    this->m_cols = indexes.size();
    this->m_matrix.resize(this->m_rows, std::vector<double>(this->m_cols, 0.f));

    unsigned j = 0;
    for (auto col : indexes) {
        for (unsigned i = 0; i < this->getRows(); i++) {
            this->m_matrix[i][j] = A->m_matrix[i][col];
        }
        j++;
    }
}

Matrix Matrix::multiply(Matrix *B) {
    if (this->getCols() != B->getRows()) exit(1); // check if it's possible to multiply

    Matrix result(this->getRows(), B->getCols());

    for (unsigned i = 0; i < this->getRows(); i++)
        for (unsigned j = 0; j < B->getCols(); j++)
            for (unsigned k = 0; k < this->getCols(); k++)
                result.m_matrix[i][j] += this->m_matrix[i][k] * B->m_matrix[k][j];

    return result;
}

Matrix Matrix::transpose() {
    Matrix result(this->getCols(), this->getRows());


    for (unsigned i = 0; i < this->getCols(); i++)
        for (unsigned j = 0; j < this->getRows(); j++)
            result.m_matrix[i][j] = this->m_matrix[j][i];

    return result;
}

double Matrix::determinant() {
    if (this->getRows() != this->getCols()) exit(2); // non-square matrix

    if (this->getCols() == 1) return this->m_matrix[0][0]; // 1-element matrix

    // base case, 2x2 matrix
    if (this->getCols() == 2)
        return this->m_matrix[0][0]*this->m_matrix[1][1] - this->m_matrix[0][1]*this->m_matrix[1][0];

    double det = 0.f;
    Matrix subMatrix(this->getRows()-1, this->getCols()-1);

    for (unsigned x = 0; x < this->getCols(); x++) {
        unsigned subi = 0;
        for (unsigned i = 1; i < this->getCols(); i++) {
            unsigned subj = 0;
            for (unsigned j = 0; j < this->getCols(); j++) {
                if (x == j) continue;

                subMatrix.m_matrix[subi][subj] = this->m_matrix[i][j];
                subj++;
            }
            subi++;
        }
        det += (pow(-1, x) * this->m_matrix[0][x] * subMatrix.determinant());
    }

    return det;
}

Matrix Matrix::cofactor(unsigned x, unsigned y) {
    Matrix result(this->getRows()-1, this->getRows()-1);

    unsigned i = 0, j = 0;

    for (unsigned row = 0; row < this->getRows(); row++) {
        for (unsigned col = 0; col < this->getRows(); col++) {
            if (row != x && y != col) {
                result.m_matrix[i][j++] = this->m_matrix[row][col];

                if (j == this->getRows()-1) {
                    j = 0;
                    i++;
                }
            }
        }
    }

    return result;
}

Matrix Matrix::adjoint() {
    Matrix result(this->getRows(), this->getRows());

    if (this->getRows() == 1) {
        result.m_matrix[0][0] = 1;
        return result;
    }

    char sign = 1;

    for (unsigned i = 0; i < this->getRows(); i++) {
        for (unsigned j = 0; j < this->getRows(); j++) {
            Matrix temp = this->cofactor(i, j);

            // if i+j is even, sign is positive
            sign = ((i+j) % 2 == 0) ? 1 : -1;

            result.m_matrix[i][j] = sign * temp.determinant();
        }
    }

    return result.transpose();
}

Matrix Matrix::inverse() {
    int det = (int) this->determinant();

    if (det == 0) exit(3);

    Matrix result(this->getRows(), this->getCols());

    /*
    A^(-1) = (1/det(A))*adj(A)
    */

   double in_det = 1.0/det;
   Matrix adj = this->adjoint();

   for (unsigned i = 0; i < this->getRows(); i++)
    for (unsigned j = 0; j < this->getCols(); j++)
        result.m_matrix[i][j] = (double)(in_det * adj.m_matrix[i][j]);

   return result;
}

std::vector<std::vector<double>> &Matrix::getMatrix() {
    return m_matrix;
}

void Matrix::readMatrix() {
    for (unsigned i = 0; i < this->getRows(); i++)
        for (unsigned j = 0; j < this->getCols(); j++)
            std::cin >> this->m_matrix[i][j];
}

void Matrix::printMatrix() {
    for (unsigned i = 0; i < this->getRows(); i++) {
        for (unsigned j = 0; j < this->getCols(); j++) {
            std::cout << this->m_matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

Matrix Matrix::getColumnMatrix(unsigned col) {
    Matrix columnMatrix(this->getRows(), 1);

    for (unsigned i = 0; i < this->getRows(); i++)
        columnMatrix.m_matrix[i][0] = this->m_matrix[i][col];


    return columnMatrix;
}

bool Matrix::isNegative() {
    for (unsigned i = 0; i < this->getRows(); i++)
        for (unsigned j = 0; j < this->getCols(); j++)
            if (this->m_matrix[i][j] > 0)
                return false;

    return true;
}

void Matrix::setNewValueAtSpecificPositionOnMatrix(unsigned i, unsigned j, double newValue) {
    this->m_matrix[i][j] = newValue;
}

void Matrix::swapColumns(Matrix *N, unsigned k, unsigned l) {
    for (unsigned i = 0; i < this->getRows(); i++) {
        double temp = this->m_matrix[i][l];
        this->m_matrix[i][l] = N->getMatrix()[i][k];
        N->setNewValueAtSpecificPositionOnMatrix(i, k, temp);
    }
}
