#pragma once
#include "Config.h"
#include "TrajTSI.h"
#include "NDArray.h"
#include <iostream>
#include <vector>
#include <array>
#include <complex.h>


using std::cout, std::vector, std::array, std::complex;

struct FTData{
    std::string FT_method = "DFT";
    vector<array<double,2>> q;

    NDArray<complex<double>> hq; // hq = (N_q,N_r,N_f)
    FTData(int N_q, int N_r, int N_f, vector<array<double,2>> q_in) :
        hq(vector<int>{N_q,N_r,N_f}), q(q_in){}

    int N_q() const {return hq.Shape(0);}
    int N_r() const {return hq.Shape(1);}
    int N_f() const {return hq.Shape(2);}
};

class FourierTransform{
    public:
        enum class FTMethod {DFT, LFST, None};

    private:
        const PosData& m_PosData;
        FTData m_FTData;
        FTMethod m_FTMethod;
        
    public:
        FourierTransform(const PosData& PosData, int N_k, std::string ftmethod = "DFT") 
            : m_PosData(PosData), m_FTData(LoadFTData(N_k)), 
            m_FTMethod(ftmethod == "DFT" ? FTMethod::DFT :
                                           FTMethod::None){}

        FTData LoadFTData(int N_k);
        
        void RunFourierTransform();
        void StaticFluctuation();
        void AutoCorrelation(int N_t);
        
        inline void DFT(int r,int f); 
        inline void LSFT(int r, int f);

        const FTData& getFTData() const {return m_FTData;}  
};