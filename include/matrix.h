#pragma once

#include <initializer_list>
#include <stdexcept>
#include <valarray>

template <typename T> requires std::is_arithmetic_v<T>
class matrix
{
public:
    matrix(const size_t size_x, const size_t size_y) : data(size_y * size_x), size_y(size_y), size_x(size_x) {}
    matrix(const size_t size_x, const size_t size_y, std::initializer_list<T> list) : data(list), size_y(size_y), size_x(size_x) 
    {
        if(size_y * size_x != list.size()) [[unlikely]] std::runtime_error("Wrong matrix dimentions");
    }

    std::valarray<T> GetRow(const size_t row) const { return data[std::slice(row * size_x, size_x, 1)]; }
    std::valarray<T> GetColumn(const size_t col) const { return data[std::slice(col, size_y, size_x)]; }

    std::slice_array<T> GetColumnSlice(const size_t col) { return data[std::slice(col, size_y, size_x)]; }
    std::slice_array<T> GetRowSlice(const size_t row) { return data[std::slice(row * size_x, size_x, 1)]; }

    T GetElement(const size_t x, const size_t y) const { return data[x + y * size_x]; }
    T& GetElement(const size_t x, const size_t y) { return data[x + y * size_x]; }

    size_t GetSizeY() const { return size_y; }
    size_t GetSizeX() const { return size_x; }
private:
    std::valarray<T> data;
    size_t size_y;
    size_t size_x;

};
