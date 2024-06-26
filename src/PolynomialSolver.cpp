#include "PolynomialSolver.h"
#include <cmath>
#include <ranges>
#include "matrix.h"
#include "vector.h"
#include "MatrixSolver.h"

real find_derivative(MFunc func, real x, real delta)
{
    return (func(x + delta) - func(x - delta))/(2 * delta);
}

real solve_polynomial(std::span<real> coefficients, const real x)
{
    real result{0.0};
    for (size_t i = 0; i < coefficients.size(); i++)
        result += pow(x, i) * coefficients[i];
    
    return result;
}

real solve_polynomial_horner(std::span<real> coefficients, const real x)
{
    real result{0.0};
    for(const auto& a : std::ranges::views::reverse(coefficients))
        result = result * x + a;
    return result;
}

real find_function_zero_with_bisection(MFunc func, real a, real b)
{
    if(fabs(a - b) < 0.00001)
        return (a + b)/2;
    else
    {
        auto x = (a + b)/2;
        if(func(x) * func(a) > 0)
            return find_function_zero_with_bisection(func, x, b);
        else
            return find_function_zero_with_bisection(func, a, x);
    }
}

real find_function_zero_with_falsi(MFunc func, real a, real b)
{
    if(fabs(func(a)) < 0.00001)
        return a;
    return find_function_zero_with_falsi(func, a - func(a) / (func(b) - func(a)) * (b - a), b);
}

real find_function_zero_with_secant(MFunc func, real a, real b)
{
    if(fabs(func(a)) < 0.00001)
        return a;
    return find_function_zero_with_secant(func, a - func(a) / (func(a) - func(b)) * (a - b), a);  
}

real find_function_zero_with_newton_raphson(MFunc func, real a, real b)
{
    real x;
    if((find_derivative(func, a-1e-6, 1e-6) - find_derivative(func, a+1e-6, 1e-6))/2e-6 * a > 0)
        x = a;
    else
        x = b;

    while (fabs(func(x)) > 1e-6 )
        x -= func(x)/find_derivative(func, x, 1e6);

    return x;
}


vector<real> get_polynomial_approximation(std::span<real> input, std::span<real> output, size_t n)
{
    size_t m = input.size();
    vector<real> polynomial = vector<real>(m);
    matrix<real> d = matrix<real>{n + 1, m};
    for(size_t y = 0; y < m; y++)
        for(size_t x = 0; x <= n; x++)
            d.GetElement(x, y) = pow(input[y], x);
    
    polynomial = numericals::solve_matrix_eq_jordan(d.Transposed() * d, d.Transposed() * vector<real>(output));

    return polynomial;
}
