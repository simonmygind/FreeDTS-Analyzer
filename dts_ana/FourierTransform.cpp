#include "FourierTransform.h"
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iostream>

using std::cerr;

FTData FourierTransform::LoadFTData(int N_k){
    int N_q = (N_k + 1) * 2 * N_k;    
    vector<array<double,2>> q(N_q);
    int idx = 0;
    for (int m = -N_k; m <= N_k; ++m) {
        for (int n = 0; n <= N_k; ++n) {
            if (m == 0 && n == 0) continue;
            if (m <= 0 && n == 0) continue;
            else{
                q[idx][0] = 2.0 * M_PI / m_PosData.boxsize[0] * (double)(m);
                q[idx][1] = 2.0 * M_PI / m_PosData.boxsize[1] * (double)(n);
                idx++;
            }
        }
    }
    std::sort(q.begin(),q.end(),[](array<double,2> a, array<double,2> b){
        return sqrt(a[0] * a[0] + a[1]*a[1]) < sqrt(b[0] * b[0] + b[1] * b[1]);
    });

    int N_r = m_PosData.N_r();
    int N_f = m_PosData.N_f();
    FTData FTData(N_q,N_r,N_f,q);
    return FTData;
};

inline void FourierTransform::DFT(int r, int f) {
    // precompute h_mean
    double h_mean = 0;
    for (int v = 0; v < m_PosData.N_v(); ++v)
        h_mean += m_PosData.pos.At({r,f,v})[2];
    h_mean /= (double)(m_PosData.N_v());

    // accumulate per vertex into local array
    vector<complex<double>> hq_local(m_FTData.N_q(),0);
    for (int v = 0; v < m_PosData.N_v(); ++v){
        double px = m_PosData.pos.At({r,f,v})[0];
        double py = m_PosData.pos.At({r,f,v})[1];
        double dh = m_PosData.pos.At({r,f,v})[2] - h_mean;
        for (int i = 0; i < m_FTData.N_q(); ++i){
            double phase = m_FTData.q[i][0]*px + m_FTData.q[i][1]*py;
            hq_local[i] += complex<double>(dh*std::cos(phase), dh*std::sin(phase));
        }
    }

    // write back contiguously
    for (int i = 0; i < m_FTData.N_q(); ++i){
        m_FTData.hq.At({i,r,f}) = hq_local[i] / std::sqrt(m_PosData.N_v());
    }
}

void FourierTransform::RunFourierTransform(){
    for (int r = 0; r < m_FTData.N_r(); r++){
        for(int f = 0; f< m_FTData.N_f(); f++){
            switch(m_FTMethod){
                case FTMethod::DFT: DFT(r,f); break;
                //case FTMethod::LSFT: LSFT(r,f); break;
            }
        }
        cout << "Finished fouriertransforming rep " << r +1 << "\n";
    }  
}

void FourierTransform::StaticFluctuation(){
    //hq = (N_q,N_r,N_f)
    NDArray<double> hq_fmean({m_FTData.N_q(),m_FTData.N_r()}); //hq_fmean = (N_q,N_r)
    for (int i = 0; i < m_FTData.N_q(); ++i)
        for (int r = 0; r < m_FTData.N_r(); ++r)
            for (int f = 0; f < m_FTData.N_f(); ++f)
                hq_fmean.At({i,r}) += std::norm(m_FTData.hq.At({i,r,f})) / (double)(m_FTData.N_f());

    vector<double> s_mean(m_FTData.N_q());
    for (int i = 0; i < m_FTData.N_q(); ++i)
        for (int r = 0; r < m_FTData.N_r(); ++r)
            s_mean[i] += hq_fmean.At({i,r}) / (double)(m_FTData.N_r());

    vector<double> s_var(m_FTData.N_q());
    for (int i = 0; i < m_FTData.N_q(); ++i)
        for (int r = 0; r < m_FTData.N_r(); ++r)
            s_var[i] += (hq_fmean.At({i,r})-s_mean[i])*(hq_fmean.At({i,r})-s_mean[i]) / (double)(m_FTData.N_r()-1);

    std::filesystem::create_directories("StaticUndulationSpectrum");
    std::ofstream outputFile("StaticUndulationSpectrum/StaticUndulationSpectrum.txt");
    outputFile << "Mean Boxsize \n";
    outputFile << "L_x L_y L_z \n";
    outputFile << m_PosData.boxsize[0] << " " << m_PosData.boxsize[1] << " " << m_PosData.boxsize[2] << "\n";
    outputFile << "q |h_q|^2 |h_q|^2_var\n";
    for (int i = 0; i < m_FTData.N_q(); ++i) {
        double q_mag = std::sqrt(m_FTData.q[i][0]*m_FTData.q[i][0] + m_FTData.q[i][1]*m_FTData.q[i][1]);
        outputFile << q_mag << " " << s_mean[i] << " " << s_var[i] << "\n";
    }
}

void FourierTransform::AutoCorrelation(int N_t) {
    if (N_t == 0){
        cerr << "Specify -N_t \n";
        exit(1);
    }
    //hq= (N_q,N_r,N_f)
    //X = (N_q,N_r,N_t) 
    NDArray<complex<double>> X({m_FTData.N_q(),m_FTData.N_r(),N_t});
    for (int i = 0; i < m_FTData.N_q(); ++i)
        for (int r=0; r < m_FTData.N_r(); ++r)    
            for (int t = 0; t < N_t; ++t)                
                for (int f = 0; f < m_FTData.N_f() - t; ++f)
                    X.At({i,r,t}) +=  m_FTData.hq.At({i,r,f+t})*std::conj(m_FTData.hq.At({i,r,f})) / (double)(m_FTData.N_f()-t);
    
    for (int i = 0; i < m_FTData.N_q(); ++i)
        for (int r = 0; r < m_FTData.N_r(); ++r)
            for (int t = 1; t < N_t; ++t)            
                X.At({i,r,t}) /= X.At({i,r,0});

    // X = a(dt)/a(0)
    for (int i = 0; i < m_FTData.N_q(); ++i)
        for (int r = 0; r < m_FTData.N_r(); ++r)
            X.At({i,r,0}) /= X.At({i,r,0});

    //Should be mean = (N_q,N_t)
    NDArray<complex<double>> a_mean({m_FTData.N_q(),N_t});
    for (int i = 0; i < m_FTData.N_q(); ++i)
        for (int r = 0; r < m_FTData.N_r(); ++r)
            for (int t = 0; t < N_t; ++t)            
                a_mean.At({i,t}) += X.At({i,r,t}) / (double)(m_FTData.N_r());

    NDArray<complex<double>> a_var({m_FTData.N_q(),N_t});
    for (int i = 0; i < m_FTData.N_q(); ++i)
        for (int r = 0; r < m_FTData.N_r(); ++r)
            for (int t = 0; t < N_t; ++t)  
                a_var.At({i,t}) += std::norm(X.At({i,r,t}) - a_mean.At({i,t}))/(double)(m_FTData.N_r()-1);
    
    std::filesystem::create_directories("AutoCorrelationFunction");
    for (int i = 0; i < m_FTData.N_q(); i++){
        std::ofstream outputFile("AutoCorrelationFunction/q" + std::to_string(i) +".txt");
        outputFile << "q"+std::to_string(i) +" \n";
        outputFile << sqrt(m_FTData.q[i][0]*m_FTData.q[i][0]+m_FTData.q[i][1]*m_FTData.q[i][1]) << " "<< 0 << " " << 0 <<"\n";
        outputFile << "t mean var \n";
        for (int t = 0; t < N_t; ++t)
            outputFile << t<< " " << std::abs(a_mean.At({i,t})) << " " << a_var.At({i,t}) << "\n";
    }
}

