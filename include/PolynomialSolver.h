#pragma once
#include "numerical_types.h"
#include <span>

real solve_polynomial(std::span<real>,const real);
real solve_polynomial_horner(std::span<real>,const real);


