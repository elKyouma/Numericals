#pragma once
#include <valarray>

template <typename T> requires std::is_arithmetic_v<T>
size_t find_index_of_column_max(const std::valarray<T>& vals, size_t start, size_t end)
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
void swap_slices(std::slice_array<T> first, std::slice_array<T> second)
{
    auto tmp = first;
    first = second;
    second = tmp; 
}
