#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <vector>

class Matrix
{
public:
    // constructor
    Matrix(unsigned, unsigned);
    Matrix(Matrix &, std::vector<int> const&);

    // getters
    inline unsigned getRows() const { return m_rows; } // get number of rows in theta(1)
    inline unsigned getCols() const { return m_cols; } // get number of cols in theta(1)
    inline std::vector<std::vector<double>> getMatrix() const { return m_matrix; }
    Matrix getColumnMatrix(unsigned);
    Matrix getBasicCoeficientsMatrix(std::vector<int> const &);

    // operators overload
    Matrix operator*(Matrix const &);

    // matrix transformations
    Matrix transpose();
    Matrix inverse();

    // i/o
    void readMatrix();
    void printMatrix();

private:
    unsigned m_rows;
    unsigned m_cols;
    std::vector<std::vector<double>> m_matrix;

    // methods used internally
    double determinant();

    Matrix cofactor(unsigned, unsigned);

    Matrix adjoint();
};

#endif // ends __MATRIX_H__