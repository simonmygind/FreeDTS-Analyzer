#include <vector>
#include <iostream>
#pragma once

using std::vector, std::cout;


template<typename T> class NDArray
{
private:
    vector<T> m_data; // ex: [A*B*C]
    vector<int> m_shape; //  Contain shapes ie. ranks ex: [A,B,C]
    vector<int> m_strides; // ex: [A*B*C,B*C,C]
public:
    NDArray(vector<int> shape);

    //coords = [i,j,k] ex: LinearIndex(3,4,2) returns (3*4+4+2)
    int LinearIndex(vector<int> coords) const;

    T& At(vector<int> coords) 
        {return m_data[LinearIndex(coords)];} 
    const T& At(vector<int> coords) const 
        {return m_data[LinearIndex(coords)];}

    int Shape(int dim) const { return m_shape[dim]; } //returns dimension size ex: Shape(1) = B
};

#include "NDArray.tpp"
