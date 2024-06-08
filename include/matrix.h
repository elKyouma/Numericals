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

    std::valarray<T> GetRow(const size_t row, const size_t offset = 0) const        
    { 
        return data[std::slice(row * size_x + offset, size_x - offset, 1)]; 
    }

    std::valarray<T> GetColumn(const size_t col, const size_t offset = 0) const     
    { 
        return data[std::slice(col + offset * size_x, size_y - offset, size_x)]; 
    }

    std::slice_array<T> GetColumnSlice(const size_t col, const size_t offset = 0)   
    { 
        return data[std::slice(col + offset * size_x, size_y - offset, size_x)]; 
    }

    std::slice_array<T> GetRowSlice(const size_t row, const size_t offset = 0)      
    { 
        return data[std::slice(row * size_x + offset, size_x - offset, 1)]; 
    }

    const std::slice_array<T> GetColumnSlice(const size_t col, const size_t offset = 0) const  
    { 
        return data[std::slice(col + offset * size_x, size_y - offset, size_x)]; 
    }

    const std::slice_array<T> GetRowSlice(const size_t row, const size_t offset = 0) const      
    { 
        return data[std::slice(row * size_x + offset, size_x - offset, 1)]; 
    }

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
    
    matrix<T> operator*(const matrix<T>& other)
    {
        if(GetSizeX() != other.GetSizeY()) [[unlikely]] std::runtime_error("Wrong matrices dimensions on multiplication operator");
        
        matrix<T> result{other.GetSizeX(), GetSizeY()};

        for(size_t i = 0; i < GetSizeX(); i++)
            for(size_t x = 0; x < GetSizeX(); x++)
                for(size_t y = 0; y < GetSizeY(); y++)
                     result.GetElement(x, y) += GetElement(i, y) * other.GetElement(x, i);
        
        return result;
    }

    matrix<T>& operator+= (const matrix<T>& other)
    {
        if(GetSizeX() != other.GetSizeX() || GetSizeY() != other.GetSizeY()) [[unlikely]] std::runtime_error("Wrong matrices dimensions on multiplication operator");
    
        for(size_t x = 0; x < GetSizeX(); x++)
            for(size_t y = 0; y < GetSizeY(); y++)
                GetElement(x, y) += other.GetElement(x, y);
        return *this;
    }

    matrix<T> operator-()
    {
        matrix<T> result{GetSizeX(), GetSizeY()};
        for(size_t x = 0; x < GetSizeX(); x++)
            for(size_t y = 0; y < GetSizeY(); y++)
                result.GetElement(x, y) = -GetElement(x, y);
        
        return result;
    }

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

