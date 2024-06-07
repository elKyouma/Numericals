#pragma once
#include <span>
#include <stack>
#include <utility>

#include "matrix.h"
#include "vector.h"

typedef float real;
//typedef double real;

typedef std::stack<std::pair<size_t, size_t>> permutation_stack;

enum MatrixFlag
{
    NORMAL,
    PARTIAL_SELECT,
    FULL_SELECT
};

//w = x^n-1*a_n + x^n-2*a_n-1 + a_0 
real solve_polynomial(std::span<real>,const real);
real solve_polynomial_horner(std::span<real>,const real);

vector<real> solve_high_trian_matrix_eq(const matrix<real>& a, const vector<real>& b, bool assumeDiagonalOnes = false);
vector<real> solve_low_trian_matrix_eq(const matrix<real>& a, const vector<real>& b, bool assumeDiagonalOnes = false);
vector<real> solve_matrix_eq_gauss( matrix<real> a, vector<real> b, MatrixFlag flag = NORMAL);
vector<real> solve_matrix_eq_jordan( matrix<real> a, vector<real> b, MatrixFlag flag = NORMAL);
vector<real> solve_matrix_eq_with_lu_decomposition(const matrix<real>& a, const vector<real>& b, MatrixFlag flag = NORMAL);

matrix<real> lu_decomposition(matrix<real> a, MatrixFlag flag = NORMAL);
matrix<real> ldl_decomposition(matrix<real> a, MatrixFlag flag = NORMAL);
