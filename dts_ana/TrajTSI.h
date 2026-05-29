#pragma once
#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include "Config.h"
#include "NDArray.h"

using std::array,std::vector;

struct PosData {
    array<double,3> boxsize;
    NDArray<array<double,3>> pos; // pos = Linear(N_r,N_f,N_v)

    PosData(int N_r, int N_f, int N_v) : 
        pos(vector<int>{N_r,N_f,N_v}), boxsize({0,0,0}){}

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

        const PosData& getPosData() const {return m_PosData;}
};