#pragma once
#include <array>
#include <cmath>
#include <ostream>

using std::array;

template<typename T,size_t N>
inline T operator*(const array<T,N>& a, const array<T,N>& b){
    T sum = 0;
    for (size_t i = 0; i < N; ++i)
        sum += a[i] * b[i];
    return sum;
}

template<typename T, size_t N>
inline array<T,N> operator*(const array<T,N>& a, T b){
    array<T,N> c;
    for (size_t i = 0; i < N; ++i)
        c[i] = a[i] * b;
    return c;
}

template<typename T, size_t N>
inline array<T,N> operator*(T b, const array<T,N>& a){ return a * b; }

template<typename T, size_t N>
inline array<T,N> operator-(const array<T,N>& a, const array<T,N>& b){
    array<T,N> c;
    for(size_t i = 0; i < N; ++i)
        c[i] = a[i] - b[i];
    
    return c; 
}


template<typename T, size_t N>
inline array<T,N>& operator+=(array<T,N>& a, const array<T,N>& b){
    for (size_t i = 0; i < N; ++i)
        a[i] += b[i];
    return a;
}

template<typename T, size_t N>
inline std::ostream& operator<<(std::ostream& osstream, const array<T,N>& a){
    for (size_t i = 0; i < N; ++i)
        osstream << a[i] << " ";
    return osstream;
}

template<size_t N>
inline double norm(const array<double,N>& v) {return std::sqrt(v*v); }

template<size_t N>
inline double norm(const array<int,N>& v) {return std::sqrt(v*v); }

