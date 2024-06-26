#include "PolynomialSolver.h"
#include <cmath>
#include <iostream>
#include <ranges>

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

real find_function_zero_with_bisection(std::function<real(real)> func, real a, real b)
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

real find_function_zero_with_falsi(std::function<real(real)> func, real a, real b)
{
    if(fabs(func(a)) < 0.00001)
        return a;
    return find_function_zero_with_falsi(func, a - func(a) / (func(b) - func(a)) * (b - a), b);
}

real find_function_zero_with_secant(std::function<real(real)> func, real a, real b)
{
    if(fabs(func(a)) < 0.00001)
        return a;
    if(!std::isnan(a))
        std::cout << a << " " << b << " " << func(a) << " " << func(b)  << "\n";
    return find_function_zero_with_secant(func, a - func(a) / (func(a) - func(b)) * (a - b), a);  
}
