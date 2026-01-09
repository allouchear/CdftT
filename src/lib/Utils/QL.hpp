#ifndef CDFTT_QL_HPP_INCLUDED
#define CDFTT_QL_HPP_INCLUDED

#include <vector>


/**
 * @brief Reduces a real, symmetric matrix to tridiagonal form using Householder transformations.
 */
void reductionToTridiagonal(std::vector<std::vector<double>>& matrix, int dimension, std::vector<double>& diagonal, std::vector<double>& subdiagonal);

int diagonalisationOfATridiagonalMatrix(double* D, double* E, int n, double** V);

#endif // CDFTT_QL_HPP_INCLUDED