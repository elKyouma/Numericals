#pragma once
#include <span>

//w = x^n-1*a_n + x^n-2*a_n-1 + a_0 
double solve_polynomial(std::span<double>);
