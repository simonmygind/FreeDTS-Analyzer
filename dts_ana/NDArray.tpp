template<typename T>
NDArray2D<T>::NDArray2D(int A, int B)
    : m_data(A * B),
      m_shape{A, B},
      m_strides(B)
{}

template<typename T>
NDArray2D<T>::NDArray2D(int A, int B, const T& initial_values)
    : m_data(A * B,initial_values),
      m_shape{A, B},
      m_strides(B)
{}

template<typename T>
NDArray3D<T>::NDArray3D(int A, int B, int C)
    : m_data(A*B*C),
      m_shape{A,B,C},
      m_strides{B*C,C}
{}

template<typename T>
NDArray3D<T>::NDArray3D(int A, int B, int C, const T& initial_values) :
    m_data(A*B*C, initial_values),
    m_shape{A,B,C},
    m_strides{B*C,C}
{}

template <typename T>
NDArray<T>::NDArray(const vector<int>& shape) : m_shape(shape){
    m_strides.resize(m_shape.size());
    m_strides.back() = 1;
    // ex: idx =  2 -> 1 -> 0
    for (int i = m_shape.size()-2; i >= 0; --i){
        m_strides[i] = m_strides[i+1] * m_shape[i+1];
    }

    int total_size = 1;
    for (int s  : m_shape){
        total_size *= s;
    }

    m_data.resize(total_size);
}

template <typename T>
NDArray<T>::NDArray(const vector<int>& shape,const T& initial_value) : m_shape(shape)
{
    m_strides.resize(m_shape.size());
    m_strides.back() = 1;
    // ex: idx =  2 -> 1 -> 0
    for (int i = m_shape.size()-2; i >= 0; --i){
        m_strides[i] = m_strides[i+1] * m_shape[i+1];
    }

    int total_size = 1;
    for (int s  : m_shape){
        total_size *= s;
    }

    m_data.resize(total_size,initial_value);
}

template <typename T>
int NDArray<T>::LinearIndex(const vector<int>& coords) const {
    int idx = 0;
    for (int i = 0; i < m_strides.size(); ++i){
        idx += m_strides[i] * coords[i]; 
    }
    return idx;
}

