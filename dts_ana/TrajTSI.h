#pragma once
#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include "Config.h"
#include "NDArray.h"

using std::array,std::vector;

struct PosData {
    NDArray3D<array<double,3>> pos; // pos = Linear(N_r,N_f,N_v)
    NDArray2D<array<double,3>> boxsize; // boxsize = Linear(N_r,N_f)
    array<double,3> boxsize_mean;

    PosData(int N_r, int N_f, int N_v) : 
        pos(N_r,N_f,N_v), boxsize(N_r,N_f),boxsize_mean({0,0,0}){}

    int N_r() const {return pos.Shape(0);}
    int N_f() const {return pos.Shape(1);}
    int N_v() const {return pos.Shape(2);}
};

class TrajTSI
{       
    private:
        const Config& m_Config;
        PosData m_PosData;

    public:
        TrajTSI(const Config& Config) : m_Config(Config), m_PosData(LoadConfig()){}

        PosData LoadConfig();
        void LoadTSI();
        void AlignTilt(); // Per-frame rotation removing global membrane tilt (mean normal -> z)

        const PosData& getPosData() const {return m_PosData;}
};