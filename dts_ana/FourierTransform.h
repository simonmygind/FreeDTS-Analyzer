#pragma once
#include "Config.h"
#include "TrajTSI.h"
#include "NDArray.h"
#include <iostream>
#include <vector>
#include <array>
#include <complex>
#include <algorithm>


using std::cout, std::vector, std::array, std::complex;

struct FTData{
    NDArray3D<complex<double>> hq; // hq = (N_r,N_f,N_q)
    vector<array<double,2>> q_mean;
    vector<array<int,2>> mn;
    FTData(int N_r, int N_f, int N_q, vector<array<double,2>> q_in, vector<array<int,2>> mn_in) :
        hq(N_r,N_f,N_q), q_mean(q_in), mn(mn_in){}

    int N_r() const {return hq.Shape(0);}
    int N_f() const {return hq.Shape(1);}
    int N_q() const {return hq.Shape(2);}
};

class FourierTransform{
    public:
        enum class FTMethod {DFT, LSFT};

    private:
        const Config& m_Config;
        const PosData& m_PosData;
        FTData m_FTData; // computed once over the full loaded PosData; StaticFluctuation/AutoCorrelation
                          // then window down to m_Config's current (init_f..final_f, init_rep..final_rep)
                          // sub-range of this same hq array on each call, so the DFT/LSFT itself is never redone
        vector<double> m_tau; // Calculates decorrelation time tau(i) for each q_i
                        // m_tau = (N_q)'
        //bool m_IsotropicFrame;
        FTMethod m_FTMethod;

        
    public:
        FourierTransform(const Config& Config, const PosData& PosData) 
            : m_Config(Config), m_PosData(PosData), m_FTData(LoadFTData()),
            m_FTMethod(Config.ftmethod == "LSFT" ? FTMethod::LSFT : FTMethod::DFT){}

        FTData LoadFTData();
        
        void RunFourierTransform();
        void StaticFluctuation();
        void AnisotropicStaticFluctuation();
        void AutoCorrelation();
        
        inline void DFT(int r,int f); 
        inline void LSFT(int r, int f);

        const FTData& getFTData() const {return m_FTData;}
};