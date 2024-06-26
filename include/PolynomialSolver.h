#pragma once
#include "numerical_types.h"
#include <functional>
#include <span>
#include <vector>
#include "vector.h"

using MFunc = std::function<real(real)>;

real find_derivative(MFunc, real, real);

real solve_polynomial(std::span<real>,const real);
real solve_polynomial_horner(std::span<real>,const real);

std::vector<real> get_chebyshev_polynomial_zeros(size_t n, real a, real b);

real find_function_zero_with_bisection(MFunc, real, real);
real find_function_zero_with_falsi(MFunc, real, real);
real find_function_zero_with_secant(MFunc, real, real);
real find_function_zero_with_newton_raphson(MFunc, real, real);

vector<real> get_polynomial_approximation(std::span<real> x, std::span<real> y, size_t n, std::vector<std::function<real(real)>> base = {});
MFunc get_lagrange_interpolation(std::span<real> x, std::span<real> y);
MFunc get_newton_interpolation(std::span<real> x, std::span<real> y);
