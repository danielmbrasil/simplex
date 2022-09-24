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
    inline double test() { return m_matrix[0][0]; }
    Matrix getColumnMatrix(unsigned);
    Matrix getBasicCoeficientsMatrix(std::vector<int> const &);

    // setters
    void setNewValueAtSpecificPositionOnMatrix(unsigned, unsigned, double);

    // operators overload
    Matrix operator*(Matrix const &);

    // matrix transformations
    Matrix transpose();
    Matrix inverse();
    void swapColumns(Matrix &, unsigned, unsigned);

    // i/o
    void readMatrix();
    void printMatrix();

    // verification
    bool isNegative();

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