#pragma once

#include <initializer_list>
#include <iostream>
#include <valarray>
#include <iostream>

template <typename T> requires std::is_arithmetic_v<T>
class vector{
public:
    vector(size_t size) : data(size){};
    vector(std::initializer_list<T> list) : data(list){}
    vector(std::valarray<T> array) : data(array){}
    size_t GetSize() const {return data.size();}

    T operator*(const vector<T>& other)
    {
        T result = 0.0;
        for(size_t i = 0; i < data.size(); i++)
            result += data[i] * other[i];
        return result;
    }

    T& operator[](const size_t index) { return data[index]; }
    T operator[](const size_t index) const { return data[index]; }
    
    operator std::valarray<T> () const {return data;}
private:
    std::valarray<T> data;
};


template <typename F> requires std::is_arithmetic_v<F>
std::ostream& operator<< (std::ostream& stream, vector<F> toPrint)
{
    for(size_t i = 0; i < toPrint.GetSize(); i++)
        stream << toPrint[i] << "\t";
    stream << std::endl; 
    return stream;
}
