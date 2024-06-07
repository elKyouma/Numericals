#pragma once
#include <valarray>
#include <matrix.h>

template <typename T> requires std::is_arithmetic_v<T>
size_t find_index_of_valarray_max(const std::valarray<T>& vals, size_t start, size_t end)
{
    size_t maxInd = start;
    T max = vals[start];
    for(size_t i = start + 1; i < end; i++)
        if(vals[i] > max)
        {
            max = vals[i];
            maxInd = i;
        }

    return maxInd; 
}

template <typename T> requires std::is_arithmetic_v<T>
size_t find_index_of_matrix_max(const matrix<T>& mat)
{
    T max = mat.GetElement(0);
    size_t ind = 0;
    size_t size = mat.GetSizeX() * mat.GetSizeY();
    for(size_t i = 1; i < size; i++)
        if(mat.GetElement(i) > max)
        {
            ind = i;
            max = mat.GetElement(i);
        }

    return ind;
}

template <typename T> requires std::is_arithmetic_v<T>
void swap_slices(std::slice_array<T> first, std::slice_array<T> second)
{
    std::valarray<T> tmp = first;
    first = second;
    second = tmp; 
}

template <typename T> requires std::is_arithmetic_v<T>
void swap_slices(std::valarray<T>& first, std::valarray<T>& second)
{
    std::valarray<T> tmp = first;
    first = second;
    second = tmp; 
}
