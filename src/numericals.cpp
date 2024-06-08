#include "numericals.h"
#include "vector.h"
#include "utils.h"

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

vector<real> solve_matrix_eq_gauss( matrix<real> a, vector<real> b, MatrixFlag flag)
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
    
    auto x = solve_high_trian_matrix_eq(a, b);
    while(!stack.empty())
    {
        std::swap(x[stack.top().first], x[stack.top().second]);
        stack.pop();
    }
    
    return x;
}

vector<real> solve_matrix_eq_jordan( matrix<real> a, vector<real> b, MatrixFlag flag)
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


matrix<real> lu_decomposition(matrix<real> a, [[maybe_unused]]MatrixFlag flag)
{
    if(a.GetSizeY() != a.GetSizeY()) [[unlikely]] std::runtime_error("Wrong matrix-vector sizes in solver");
    size_t size = a.GetSizeX();
    
    for(size_t j = 0; j < size; j++)
    {    
        for(size_t i = j+1; i < size; i++)
        {
            const real multiplier = a.GetElement(j, i) / a.GetElement(j, j); 
            a.GetRowSlice(i, j) -= multiplier * a.GetRow(j, j);
            a.GetElement(j, i) = multiplier;
        }
    }
    return a;
}

vector<real> solve_matrix_eq_with_lu_decomposition(const matrix<real>& a, const vector<real>& b, MatrixFlag flag)
{
    constexpr bool assume_diagonal_ones = true;
    matrix<real> lu = lu_decomposition(a, flag);
    vector<real> y = solve_low_trian_matrix_eq(lu, b, assume_diagonal_ones);
    return solve_high_trian_matrix_eq(lu, std::move(y));
}

matrix<real> ldlt_decomposition(const matrix<real>& a, [[maybe_unused]] MatrixFlag flag)
{
    if(a.GetSizeY() != a.GetSizeY()) [[unlikely]] std::runtime_error("Wrong matrix-vector sizes in ldlt decomposition");
    
    size_t size = a.GetSizeX();
    matrix<real> result {size, size};
    result.GetElement(0, 0) = a.GetElement(0, 0); 
    for(size_t i = 1; i < size; i++)
    {
        for(size_t j = 0; j < i; j++)
        {
            auto& l = result.GetElement(j, i);
            auto& c = result.GetElement(i, j);
            
            l = a.GetElement(j, i);
            for(size_t k = 0; k < j; k++)
                l -= result.GetElement(i, k) * result.GetElement(k, j);
            l /= result.GetElement(j, j);

            c = result.GetElement(j, j) * l;
            result.GetElement(i, i) = a.GetElement(i, i);
            
            for(size_t k = 0; k < i; k++)
                result.GetElement(i, i) -= result.GetElement(i, k) * result.GetElement(k, i);       
        }   
    }
    
    for(size_t i = 0; i < size; i++)
        for(size_t j = i + 1; j < size; j++)
            result.GetElement(j, i) = result.GetElement(i, j);
    return result;
}

vector<real> solve_matrix_eq_with_ldlt_decomposition(const matrix<real>& a, const vector<real>& b, MatrixFlag flag)
{
    constexpr bool assume_diagonal_ones = true;
    matrix<real> ldlt = ldlt_decomposition(a);
    vector<real> z = solve_low_trian_matrix_eq(ldlt, b, assume_diagonal_ones);
    size_t size = z.GetSize();
    vector<real> y = vector<real>(size);
    
    for(size_t i = 0; i < size; i++)
        y[i] = z[i] / ldlt.GetElement(i, i);
    
    return solve_high_trian_matrix_eq(ldlt, y, assume_diagonal_ones);
}


matrix<real> ldl_decomposition(const matrix<real>& a, MatrixFlag flag) 
{
    if(a.GetSizeY() != a.GetSizeY()) [[unlikely]] std::runtime_error("Wrong matrix-vector sizes in ldl decomposition");
 
    size_t size = a.GetSizeX();
    matrix<real> result{size, size};
    for(size_t i = 0; i < size; i++)
    {    
        auto& d = result.GetElement(i, i); 
        d = a.GetElement(i, i);
        for(size_t k = 0; k < i - 1; k++)
            d -= pow(result.GetElement(k, i), 2);
        d = sqrt(d);
        
        for(size_t j = 0; j < i - 1; j++)
        {    
            result.GetElement(j , i) = a.GetElement(j, i);
            for(size_t k = 0; k < i - 1; k++)
                result.GetElement(j, i) -= result.GetElement(j, k) * result.GetElement(i, k);
        }
    }

    for(size_t i = 0; i < size; i++)
        for(size_t j = i + 1; j < size; j++)
            result.GetElement(j, i) = result.GetElement(i, j);

    return result;
}
