#include "numericals.h"
#include "vector.h"
#include "utils.h"

#include <cmath>
#include <cstddef>
#include <iterator>
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
    if(a.GetSizeY() != a.GetSizeY() || a.GetSizeX() != b.GetSize()) [[unlikely]] std::runtime_error("Wrong matrix-vector sizes in solver");

    vector<real> x(a.GetSizeX());
    size_t last = x.GetSize() - 1;
    for(size_t i = 0; i < x.GetSize(); i++)
    {    
        size_t index = last - i;
        x[index] = b[index];
        for(size_t j = 1; j <= i; j++)
            x[index] -= a.GetElement(index + j, index) * x[index + j]; 
        x[index] /= a.GetElement(index, index);
    }
    
    return x;
}


vector<real> solve_matrix_equation_gauss( matrix<real> a, vector<real> b, MatrixFlag flag)
{
    if(a.GetSizeY() != a.GetSizeY() || a.GetSizeX() != b.GetSize()) [[unlikely]] std::runtime_error("Wrong matrix-vector sizes in solver");

    size_t size_y = a.GetSizeX(); 
    for(size_t y = 0; y < size_y; y++)
    {
        switch (flag) {
            case PARTIAL_SELECT:
                size_t maxInd = find_index_of_column_max<real>(a.GetColumnSlice(y), y, size_y);               
                if(maxInd == y) break;
                swap_slices(a.GetRowSlice(y), a.GetRowSlice(maxInd));
                break;
            //case FULL_SELECT:
            //    break;
            //case NORMAL:
            //    break;
        }

        b[y] /= a.GetElement(y, y);
        a.GetRowSlice(y) = a.GetRow(y) / a.GetElement(y, y);
        for (size_t i = y + 1; i < size_y; i++)
        {   
            b[i] -= b[y] * a.GetElement(y, i);
            a.GetRowSlice(i) -= a.GetElement(y, i) * a.GetRow(y);
        }
    }

    return solve_triangular_matrix_equation(a, b);
}

vector<real> solve_matrix_equation_jordan( matrix<real> a, vector<real> b, [[maybe_unused]] MatrixFlag flag)
{
    if(a.GetSizeY() != a.GetSizeY() || a.GetSizeX() != b.GetSize()) [[unlikely]] std::runtime_error("Wrong matrix-vector sizes in solver");

    size_t size_y = a.GetSizeX(); 
    for(size_t y = 0; y < size_y; y++)
    {
        b[y] /= a.GetElement(y, y);
        a.GetRowSlice(y) = a.GetRow(y) / a.GetElement(y, y);
        for (size_t i = 0; i < size_y; i++)
        {   
            if(i == y) continue;

            b[i] -= b[y] * a.GetElement(y, i);
            a.GetRowSlice(i) -= a.GetElement(y, i) * a.GetRow(y);
        }
    }

    return b;
}

