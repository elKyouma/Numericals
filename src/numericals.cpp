#include "numericals.h"
#include "vector.h"
#include <cmath>
#include <cstddef>
#include <ranges>
#include <stdexcept>

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

vector<real> solve_triangular_matrix_equation(const matrix<real>& a, const vector<real>& b)
{
    if(a.GetSizeY() != a.GetSizeY() || a.GetSizeX() != b.size()) [[unlikely]] std::runtime_error("Wrong matrix-vector sizes in solver");

    vector<real> x(a.GetSizeX());
    size_t last = x.size() - 1;
    for(size_t i = 0; i < x.size(); i++)
    {    
        size_t index = last - i;
        x[index] = b[index];
        for(size_t j = 1; j <= i; j++)
            x[index] -= a.GetElement(index + j, index) * x[index + j]; 
        x[index] /= a.GetElement(index, index);
    }
    
    return x;
}
