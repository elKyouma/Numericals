#pragma once
#include "numerical_types.h"
#include "matrix.h"
#include "vector.h"
#include <utility>
#include "utils.h"

class PivotingStrategy
{
public:
    virtual void PreIteration(matrix<real>& A, vector<real>& b, const size_t i) = 0;
    virtual void CleanUp(vector<real>& x) = 0;
};

class NoPivotingStragegy : public PivotingStrategy
{
public:
    void PreIteration([[maybe_unused]]matrix<real>& A, 
                      [[maybe_unused]]vector<real>& b, 
                      [[maybe_unused]]const size_t i) override {}
    void CleanUp([[maybe_unused]]vector<real>& x) override {}
};

class PartialPivotingStragegy : public PivotingStrategy
{
public:
    void PreIteration(matrix<real>& A, vector<real>& b, const size_t i) override
    {
        size_t maxInd = find_index_of_valarray_max<real>(A.GetColumnSlice(i), i, A.GetSizeY());               
        if(maxInd == i) return;
        
        swap_slices(A.GetRowSlice(i), A.GetRowSlice(maxInd));
        std::swap(b[i], b[maxInd]);
    }
    void CleanUp([[maybe_unused]] vector<real>& x) override{}
};

class FullPivotingStragegy : public PivotingStrategy
{
public:
    void PreIteration(matrix<real>& A, vector<real>& b, const size_t i) override
    {
        auto [maxIndx, maxIndy] = find_index_of_matrix_max(A, i, i);
                    
        if(maxIndy == i && maxIndx == i) return;
                    
        swap_slices(A.GetRowSlice(i), A.GetRowSlice(maxIndy));
        std::swap(b[i], b[maxIndy]);
        swap_slices(A.GetColumnSlice(i), A.GetColumnSlice(maxIndx));
        stack.push({i, maxIndx});
 
    }

    void CleanUp(vector<real>& x) override
    {
        while(!stack.empty())
        {
            std::swap(x[stack.top().first], x[stack.top().second]);
            stack.pop();
        }
    }

private:
    permutation_stack stack;
};
