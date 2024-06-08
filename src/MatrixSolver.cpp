#include "MatrixSolver.h"
#include "PivotingStrategy.h"
#include "vector.h"
#include "MatrixDecomposer.h"

#include <cmath>
#include <cstddef>
#include <ranges>
#include <stdexcept>

namespace numericals{

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

vector<real> solve_high_trian_matrix_eq(const matrix<real>& a, const vector<real>& b, bool assumeDiagonalOnes)
{
    if(a.GetSizeY() != a.GetSizeY() || a.GetSizeX() != b.GetSize()) [[unlikely]] std::runtime_error("Wrong matrix-vector sizes in solver");

    vector<real> x(a.GetSizeX());
    size_t last = x.GetSize() - 1;
    for(size_t i = 0; i < x.GetSize(); i++)
    {    
        size_t index = last - i;
        x[index] = b[index];
        for(size_t j = 1; j <= i; j++)
            x[index] -= a.GetElement(index + j, index) * x[index + j]; 
        if(!assumeDiagonalOnes)
            x[index] /= a.GetElement(index, index);
    }
    
    return x;
}

vector<real> solve_low_trian_matrix_eq(const matrix<real>& a, const vector<real>& b, bool assumeDiagonalOnes)
{
    if(a.GetSizeY() != a.GetSizeY() || a.GetSizeX() != b.GetSize()) [[unlikely]] std::runtime_error("Wrong matrix-vector sizes in solver");

    vector<real> x(a.GetSizeX());
    for(size_t i = 0; i < x.GetSize(); i++)
    {    
        x[i] = b[i];
        for(size_t j = 1; j <= i; j++)
            x[i] -= a.GetElement(i - j, i) * x[i - j]; 
        if(!assumeDiagonalOnes)
            x[i] /= a.GetElement(i, i);
    }
    
    return x;
}

vector<real> solve_matrix_eq_gauss( matrix<real> a, vector<real> b, PivotingStrategy&& strategy)
{
    if(a.GetSizeY() != a.GetSizeY() || a.GetSizeX() != b.GetSize()) [[unlikely]] std::runtime_error("Wrong matrix-vector sizes in solver");
    
    permutation_stack stack;
    size_t size_y = a.GetSizeX(); 
    for(size_t d = 0; d < size_y; d++)
    {
        strategy.PreIteration(a, b, d);
        b[d] /= a.GetElement(d, d);
        a.GetRowSlice(d) = a.GetRow(d) / a.GetElement(d, d);
        for (size_t i = d + 1; i < size_y; i++)
        {   
            b[i] -= b[d] * a.GetElement(d, i);
            a.GetRowSlice(i) -= a.GetElement(d, i) * a.GetRow(d);
        }    
    }

    auto x = solve_high_trian_matrix_eq(a, b);
    strategy.CleanUp(x); 
    return x;
}

vector<real> solve_matrix_eq_jordan( matrix<real> a, vector<real> b, PivotingStrategy&& strategy)
{
    if(a.GetSizeY() != a.GetSizeY() || a.GetSizeX() != b.GetSize()) [[unlikely]] std::runtime_error("Wrong matrix-vector sizes in solver");
    
    size_t size_y = a.GetSizeX(); 
    for(size_t d = 0; d < size_y; d++)
    {
        strategy.PreIteration(a, b, d);
        b[d] /= a.GetElement(d, d);
        a.GetRowSlice(d) = a.GetRow(d) / a.GetElement(d, d);
        for (size_t i = 0; i < size_y; i++)
        {   
            if(i == d) continue;

            b[i] -= b[d] * a.GetElement(d, i);
            a.GetRowSlice(i) -= a.GetElement(d, i) * a.GetRow(d);
        }
    }
    
    strategy.CleanUp(b);
    
    return b;
}

vector<real> solve_matrix_eq_with_lu_decomposition(const matrix<real>& a, const vector<real>& b, PivotingStrategy&& strategy)
{
    constexpr bool assume_diagonal_ones = true;
    matrix<real> lu = lu_decomposition(a, std::move(strategy));
    vector<real> y = solve_low_trian_matrix_eq(lu, b, assume_diagonal_ones);
    return solve_high_trian_matrix_eq(lu, std::move(y));
}

vector<real> solve_matrix_eq_with_ldlt_decomposition(const matrix<real>& a, const vector<real>& b, PivotingStrategy&& strategy)
{
    constexpr bool assume_diagonal_ones = true;
    matrix<real> ldlt = ldlt_decomposition(a, std::move(strategy));
    vector<real> z = solve_low_trian_matrix_eq(ldlt, b, assume_diagonal_ones);
    size_t size = z.GetSize();
    vector<real> y = vector<real>(size);
    
    for(size_t i = 0; i < size; i++)
        y[i] = z[i] / ldlt.GetElement(i, i);
    
    return solve_high_trian_matrix_eq(ldlt, y, assume_diagonal_ones);
}
}
