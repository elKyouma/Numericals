#include "MatrixDecomposer.h"
#include "PivotingStrategy.h"

namespace numericals
{

matrix<real> llt_decomposition(const matrix<real>& a) 
{
    if(a.GetSizeY() != a.GetSizeY()) [[unlikely]] std::runtime_error("Wrong matrix-vector sizes in ldl decomposition");
 
    size_t size = a.GetSizeX();
    matrix<real> result{size, size};
    for(int i = 0; i < (int)size; i++)
    {
        for(int j = 0; j < i; j++)
        {
            result.GetElement(j, i) = a.GetElement(j, i);
            for(int k = 0; k < j; k++)
               result.GetElement(j, i) -= result.GetElement(k, j) * result.GetElement(k, i);
            result.GetElement(j, i) /= result.GetElement(j, j);
        }
        result.GetElement(i, i) = a.GetElement(i, i);
        for(int k = 0; k < i; k++)
            result.GetElement(i, i) -= pow(result.GetElement(k, i), 2);
        result.GetElement(i, i) = sqrt(result.GetElement(i, i));
    } 

    for(size_t i = 0; i < size; i++)
        for(size_t j = i + 1; j < size; j++)
            result.GetElement(j, i) = result.GetElement(i, j);
    return result;
}

matrix<real> ldlt_decomposition(const matrix<real>& a)
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

matrix<real> lu_decomposition(matrix<real> a, PivotingStrategy&& strategy)
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

}
