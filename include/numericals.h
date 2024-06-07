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

//returns permutation stack
void revertPermutations(matrix<real>& a, vector<real>& b, permutation_stack& permutations);

vector<real> solve_triangular_matrix_equation(const matrix<real>& a, const vector<real>& b);
vector<real> solve_matrix_equation_gauss( matrix<real> a, vector<real> b, MatrixFlag flag = NORMAL);
vector<real> solve_matrix_equation_jordan( matrix<real> a, vector<real> b, MatrixFlag flag = NORMAL);

matrix<real> lu_decomposition(matrix<real> a);
