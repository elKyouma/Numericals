#pragma once

#include <initializer_list>
#include <ostream>
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
    T GetElement(const size_t i) const { return data[i]; }
    T& GetElement(const size_t i) { return data[i]; }


    size_t GetSizeY() const { return size_y; }
    size_t GetSizeX() const { return size_x; }
   
    auto begin() { return data.begin(); }
    auto end() { return data.end(); }
    auto begin() const { return data.begin(); }
    auto end() const { return data.end(); }

private:
    std::valarray<T> data;
    size_t size_y;
    size_t size_x;

};

template <typename T> requires std::is_arithmetic_v<T>
std::ostream& operator<< (std::ostream& stream, matrix<T> toPrint)
{
    for(size_t y = 0; y < toPrint.GetSizeY(); y++)
    {
        for(size_t x = 0; x < toPrint.GetSizeX(); x++)
            stream << toPrint.GetElement(x, y) << '\t';
        stream << '\n';
    }

    return stream;
}

