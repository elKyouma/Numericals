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
    
    permutation_stack stack;
    size_t size_y = a.GetSizeX(); 
    for(size_t d = 0; d < size_y; d++)
    {
        switch (flag) {
            case PARTIAL_SELECT:
                {
                    size_t maxInd = find_index_of_valarray_max<real>(a.GetColumnSlice(d), d, size_y);               
                    if(maxInd == d) break;
                    swap_slices(a.GetRowSlice(d), a.GetRowSlice(maxInd));
                    std::swap(b[d], b[maxInd]);
                }
                break;
            case FULL_SELECT:
                {
                    auto [maxIndx, maxIndy] = find_index_of_matrix_max(a, d, d);
                    
                    if(maxIndy == d && maxIndx == d) break;
                    
                    swap_slices(a.GetRowSlice(d), a.GetRowSlice(maxIndy));
                    std::swap(b[d], b[maxIndy]);
                    swap_slices(a.GetColumnSlice(d), a.GetColumnSlice(maxIndx));
                    stack.push({d, maxIndx});
                }
                break;
            case NORMAL:
                break;
        }

        b[d] /= a.GetElement(d, d);
        a.GetRowSlice(d) = a.GetRow(d) / a.GetElement(d, d);
        for (size_t i = d + 1; i < size_y; i++)
        {   
            b[i] -= b[d] * a.GetElement(d, i);
            a.GetRowSlice(i) -= a.GetElement(d, i) * a.GetRow(d);
        }    
    }
    
    auto x = solve_triangular_matrix_equation(a, b);
    while(!stack.empty())
    {
        std::swap(x[stack.top().first], x[stack.top().second]);
        stack.pop();
    }
    
    return x;
}

vector<real> solve_matrix_equation_jordan( matrix<real> a, vector<real> b, MatrixFlag flag)
{
    if(a.GetSizeY() != a.GetSizeY() || a.GetSizeX() != b.GetSize()) [[unlikely]] std::runtime_error("Wrong matrix-vector sizes in solver");

    permutation_stack stack;
    
    size_t size_y = a.GetSizeX(); 
    for(size_t d = 0; d < size_y; d++)
    {
         switch (flag) {
            case PARTIAL_SELECT:
                {
                    size_t maxInd = find_index_of_valarray_max<real>(a.GetColumnSlice(d), d, size_y);               
                    if(maxInd == d) break;
                    swap_slices(a.GetRowSlice(d), a.GetRowSlice(maxInd));
                    std::swap(b[d], b[maxInd]);
                }
                break;
            case FULL_SELECT:
                {
                    auto [maxIndx, maxIndy] = find_index_of_matrix_max(a, d, d);
                    
                    if(maxIndy == d && maxIndx == d) break;
                    
                    swap_slices(a.GetRowSlice(d), a.GetRowSlice(maxIndy));
                    std::swap(b[d], b[maxIndy]);
                    swap_slices(a.GetColumnSlice(d), a.GetColumnSlice(maxIndx));
                    stack.push({d, maxIndx});
                }
                break;
            case NORMAL:
                break; 
        }

       b[d] /= a.GetElement(d, d);
        a.GetRowSlice(d) = a.GetRow(d) / a.GetElement(d, d);
        for (size_t i = 0; i < size_y; i++)
        {   
            if(i == d) continue;

            b[i] -= b[d] * a.GetElement(d, i);
            a.GetRowSlice(i) -= a.GetElement(d, i) * a.GetRow(d);
        }
    }

    while(!stack.empty())
    {
        std::swap(b[stack.top().first], b[stack.top().second]);
        stack.pop();
    }
    
    return b;
}


matrix<real> lu_decomposition(matrix<real> a)
{

    if(a.GetSizeY() != a.GetSizeY()) [[unlikely]] std::runtime_error("Wrong matrix-vector sizes in solver");
    size_t size = a.GetSizeX();
    
    for(size_t j = 0; j < size; j++)
    {    
        for(size_t i = j + 1; i < size; i++)
        {    
            const real multiplier = a.GetElement(j, i) / a.GetElement(j, j); 
            a.GetRowSlice(i, j) -= multiplier * a.GetRow(j, j);
            a.GetElement(j, i) = multiplier;
        }
    }
    return a;
}
