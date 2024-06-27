#include "PolynomialSolver.h"
#include <cmath>
#include <ranges>
#include <utility>
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

std::vector<real> get_chebyshev_polynomial_zeros(size_t n, real a, real b)
{
    std::vector<real> result(n);
   
    for(size_t i = 0; i < n; i++)
        result[i] = ((a - b) * cos((2*i + 1.0) / (2*(n-1) + 2.0) * M_PI) + (a + b))*0.5f;
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


vector<real> get_polynomial_approximation(std::span<real> input, std::span<real> output, size_t n, std::vector<std::function<real(real)>> base)
{
    size_t m = input.size();
    vector<real> polynomial(m);
    matrix<real> d(n + 1, m);
    for(size_t y = 0; y < m; y++)
        for(size_t x = 0; x <= n; x++)
            if(base.empty())
                d.GetElement(x, y) = pow(input[y], x);
            else
                d.GetElement(x, y) = base[x](input[y]);

    polynomial = numericals::solve_matrix_eq_jordan(d.Transposed() * d, d.Transposed() * vector<real>(output));

    return polynomial;
}

MFunc get_lagrange_interpolation(std::span<real> input, std::span<real> output)
{
    size_t m = input.size();

    return [m, input, output](real x)
    {
        real sum = 0.0;
        for(size_t i = 0; i < m; i++)
        {
            real mult = 1.0;
            for(size_t j = 0; j < m; j++)
                if(i != j)
                    mult *= (x - input[j])/(input[i] - input[j]);
            sum += mult * output[i];
        }
        
        return sum;
    };
}

MFunc get_newton_interpolation(std::span<real> input, std::span<real> output)
{
    std::vector<std::vector<real>> coeff; // [a[], f[], f'[], f''[],...]

    size_t n = input.size();
    coeff.push_back({});
    coeff.push_back({});
    for(size_t i = 0; i < n; i++)
    {   
        coeff[0].push_back(input[i]);
        coeff[1].push_back(output[i]);
    }

    for(size_t i = 2; i <= n; i++)
    {   
        coeff.push_back(std::vector<real>(n - i + 1));
        for(size_t j = 0; j < n - i + 1; j++)
            coeff[i][j] = (coeff[i - 1][j + 1] - coeff[i - 1][j]) / (coeff[0][i - 1 + j] - coeff[0][j]);
    }
    return [coeff = std::move(coeff), n](real x)
    {
        real sum = 0.0;
        for(size_t i = 0; i < n; i++)
        {
            real mul = coeff[i + 1][0];
            for(size_t j = 0; j < i; j++)
                mul *= (x - coeff[0][j]);
           
            sum += mul;
        }
        return sum;
    };
}
