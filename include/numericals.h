#pragma once
#include <span>

typedef float real;
//typedef double real;


//w = x^n-1*a_n + x^n-2*a_n-1 + a_0 
real solve_polynomial(std::span<real>,const real);
real solve_polynomial_horner(std::span<real>,const real);
real solve_matrix_equation();
