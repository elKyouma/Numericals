#pragma once
#include <utility>
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
std::pair<size_t, size_t> find_index_of_matrix_max(const matrix<T>& mat, const size_t startX = 0, const size_t startY = 0,
                                                   size_t endX = 0, size_t endY = 0)
{
    if(endX == 0) endX = mat.GetSizeX();
    if(endY == 0) endY = mat.GetSizeY();

    T max = mat.GetElement(0);
    size_t indX = startX;
    size_t indY = startY;
    for(size_t y = startY; y < endY; y++)
        for(size_t x = startX; x < endX; x++)
            if(mat.GetElement(x, y) > max)
            {
                indX = x;
                indY = y;
                max = mat.GetElement(x, y);
            }

    return {indX, indY};
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
