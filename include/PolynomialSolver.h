#pragma once
#include "numerical_types.h"
#include <functional>
#include <span>

real solve_polynomial(std::span<real>,const real);
real solve_polynomial_horner(std::span<real>,const real);

real find_function_zero_with_bisection(std::function<real(real)>, real, real);
real find_function_zero_with_falsi(std::function<real(real)>, real, real);
