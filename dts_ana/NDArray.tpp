template <typename T>
NDArray<T>::NDArray(vector<int> shape) : m_shape(shape){
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
int NDArray<T>::LinearIndex(vector<int> coords) const {
    int idx = 0;
    for (int i = 0; i < m_strides.size(); ++i){
        idx += m_strides[i] * coords[i]; 
    }
    return idx;
}

