#pragma once

#include <initializer_list>
#include <valarray>

template <typename T> requires std::is_integral_v<T>
class vector{
public:
    vector(size_t size) : data(size){};
    vector(std::initializer_list<T> list) : data(list){}
    size_t size() {return data.size();}

    T& operator[](const size_t index) { return data[index]; }
    T operator[](const size_t index) const { return data[index]; }
    
private:
    std::valarray<T> data;
};

