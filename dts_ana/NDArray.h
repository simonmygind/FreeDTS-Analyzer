#include <vector>
#include <iostream>
#pragma once

using std::vector, std::array;

template<typename T> class NDArray2D
{
    private:
        vector<T> m_data;
        array<int,2> m_shape;
        int m_strides;
    public:
        NDArray2D(int A, int B);
        NDArray2D(int A, int B, const T& initial_value);

        inline T& operator()(int i, int j) 
            {return m_data[i * m_strides + j];}
        inline const T& operator()(int i, int j) const 
            {return m_data[i * m_strides + j];}

        int Shape(int dim) const {return m_shape[dim];}
};

template<typename T> class NDArray3D
{
    private:
        vector<T> m_data;
        array<int,3> m_shape;
        array<int,2> m_strides;
    public:
        NDArray3D(int A, int B, int C);
        NDArray3D(int A, int B, int C, const T& initial_value);

        inline T& operator()(int i, int j, int k) 
            {return m_data[i * m_strides[0] + j * m_strides[1] + k];}
        inline const T& operator()(int i, int j, int k) const 
            {return m_data[i * m_strides[0] + j * m_strides[1] + k];}

        inline int Shape(int dim) const {return m_shape[dim];}
};

template<typename T> class NDArray
{
private:
    vector<T> m_data; // ex: [A*B*C]
    vector<int> m_shape; //  Contain shapes ie. ranks ex: [A,B,C]
    vector<int> m_strides; // ex: [A*B*C,B*C,C]
public:
    NDArray(const vector<int>& shape);
    NDArray(const vector<int>& shape, const T& initial_value);

    //coords = [i,j,k] ex: LinearIndex(3,4,2) returns (3*4+4+2)
    int LinearIndex(const vector<int>& coords) const;
    
    T& AtDirect(int i, int j) 
        {return m_data[i * m_strides[0] + j];}

    const T& AtDirect(int i, int j) const 
        {return m_data[i * m_strides[0] + j];}

    int Shape(int dim) const { return m_shape[dim]; } //returns dimension size ex: Shape(1) = B
};

#include "NDArray.tpp"
