#pragma once

#include <initializer_list>
#include <stdexcept>
#include <valarray>

template <typename T> requires std::is_integral_v<T>
class matrix
{
public:
    matrix(size_t cols, size_t rows) : data(rows * cols), rows(rows), cols(cols) {}
    matrix(size_t cols, size_t rows, std::initializer_list<T> list) : data(list), rows(rows), cols(cols) 
    {
        if(rows * cols != list.size()) [[unlikely]] std::runtime_error("Wrong matrix dimentions");
    }

    std::valarray<T>& GetRow(size_t row) const { return data[std::slice(row * cols, cols, 1)]; }
    std::valarray<T>& GetColumn(size_t col) const { return data[std::slice(col, rows, cols)]; }

    std::slice_array<T>& GetColumn(size_t col) { return data[std::slice(col, rows, cols)]; }
    std::slice_array<T>& GetRow(size_t row) { return data[std::slice(row * cols, cols, 1)]; }

private:
    std::valarray<T> data;
    size_t rows;
    size_t cols;

};
