#include "FourierTransform.h"
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include "omp.h" 
#include "VecMath.h"
#include <chrono>

using std::cerr;

FTData FourierTransform::LoadFTData(){
    int N_r = m_PosData.N_r();
    int N_f = m_PosData.N_f();
    vector<array<int,2>> mn_in;
    //q_mean.reserve(N_q);

    /*
    if (std::abs(m_PosData.boxsize_mean[0] - m_PosData.boxsize_mean[1]) < 1e-6){
        cout << "Isotropic Fourier Transform\n";
        for (int m = -m_Config.N_k; m <= m_Config.N_k; ++m) {
            for (int n = 0; n <= m_Config.N_k; ++n) {
                if (m == 0 && n == 0) continue;
                if (m <= 0 && n == 0) continue;
                if (m*m + n*n > m_Config.N_k*m_Config.N_k) continue;
                else{
                    mn_in.push_back({m,n});
                }
            }
        }
        int N_q = mn_in.size();
        std::sort(mn_in.begin(),mn_in.end(),[](array<int,2> a, array<int,2> b){
            return sqrt(a[0] * a[0] + a[1]*a[1]) < sqrt(b[0] * b[0] + b[1] * b[1]);
        });

        vector<array<double,2>> q_in(N_q,{0.0,0.0});
        for (int r = 0; r < N_r; r++)
            for (int f = 0; f < N_f; f++)
                for(int i = 0; i < N_q; i++){
                    q_in[i][0] += 2*M_PI*(double)(mn_in[i][0])/(double)(m_PosData.boxsize_mean[0]*N_f*N_r);
                    q_in[i][1] += 2*M_PI*(double)(mn_in[i][1])/(double)(m_PosData.boxsize_mean[1]*N_f*N_r);
                }
        FTData FTData(N_r,N_f,N_q,q_in,mn_in);
        return FTData;
    }
    */
    if (std::abs(m_PosData.boxsize_mean[0] - m_PosData.boxsize_mean[1]) < 1e-6){
        cout << "Isotropic Fourier Transform\n";
        double q_max = 2*M_PI*double(m_Config.N_k) / std::min(m_PosData.boxsize_mean[0], m_PosData.boxsize_mean[1]);
        double q_max2 = q_max * q_max;
        struct QCandidate { double q2; double qx, qy; int m, n; };     // Collect candidates as (q2, qx, qy, m, n) so sorting keeps everything in sync
        vector<QCandidate> candidates;

        for (int m = -m_Config.N_k; m <= m_Config.N_k; ++m) {
            for (int n = 0; n <= m_Config.N_k; ++n) {
            //for (int n = -m_Config.N_k; n <= m_Config.N_k; ++n) {
                double qx = 2*M_PI*(double)(m) / m_PosData.boxsize_mean[0];
                double qy = 2*M_PI*(double)(n) / m_PosData.boxsize_mean[1];
                double q2 = qx*qx + qy*qy;
                if (m == 0 && n == 0) continue;
                if (m <= 0 && n == 0) continue; // keep same half-plane convention as isotropic case
                if (q2 > q_max2) continue;
                candidates.push_back({q2, qx, qy, m, n});
            }
        }

        std::sort(candidates.begin(), candidates.end(),
                [](const QCandidate& a, const QCandidate& b){ return a.q2 < b.q2; });

        int N_q = candidates.size();
        vector<array<double,2>> q_in(N_q);
        vector<array<int,2>> mn_in(N_q);
        for (int i = 0; i < N_q; ++i) {
            q_in[i]  = {candidates[i].qx, candidates[i].qy};
            mn_in[i] = {candidates[i].m,  candidates[i].n};
        }

        FTData FTData(N_r, N_f, N_q, q_in, mn_in);
        return FTData;
    }

    else{
        cout << "Anisotropic Fourier Transform \n";
        double q_max = 2*M_PI*double(m_Config.N_k) / std::min(m_PosData.boxsize_mean[0], m_PosData.boxsize_mean[1]);
        double q_max2 = q_max * q_max;
        struct QCandidate { double q2; double qx, qy; int m, n; };     // Collect candidates as (q2, qx, qy, m, n) so sorting keeps everything in sync
        vector<QCandidate> candidates;

        for (int m = -m_Config.N_k; m <= m_Config.N_k; ++m) {
            for (int n = 0; n <= m_Config.N_k; ++n) {
            //for (int n = -m_Config.N_k; n <= m_Config.N_k; ++n) {
                double qx = 2*M_PI*(double)(m) / m_PosData.boxsize_mean[0];
                double qy = 2*M_PI*(double)(n) / m_PosData.boxsize_mean[1];
                double q2 = qx*qx + qy*qy;
                if (m == 0 && n == 0) continue;
                if (m <= 0 && n == 0) continue; // keep same half-plane convention as isotropic case
                if (q2 > q_max2) continue;
                candidates.push_back({q2, qx, qy, m, n});
            }
        }

        std::sort(candidates.begin(), candidates.end(),
                [](const QCandidate& a, const QCandidate& b){ return a.q2 < b.q2; });

        int N_q = candidates.size();
        vector<array<double,2>> q_in(N_q);
        vector<array<int,2>> mn_in(N_q);
        for (int i = 0; i < N_q; ++i) {
            q_in[i]  = {candidates[i].qx, candidates[i].qy};
            mn_in[i] = {candidates[i].m,  candidates[i].n};
        }

        FTData FTData(N_r, N_f, N_q, q_in, mn_in);
        return FTData;
    }
};


inline void FourierTransform::DFT(int r, int f) {
    vector<array<double,2>> q(m_FTData.mn.size());
    for (int i = 0; i < m_FTData.N_q(); ++i){
        q[i][0] = 2*M_PI/m_PosData.boxsize(r,f)[0]*(double)(m_FTData.mn[i][0]);
        q[i][1] = 2*M_PI/m_PosData.boxsize(r,f)[1]*(double)(m_FTData.mn[i][1]); 
    }

    // precompute h_mean
    double h_mean = 0;
    for (int v = 0; v < m_PosData.N_v(); ++v)
        h_mean += m_PosData.pos(r,f,v)[2];
    h_mean /= (double)(m_PosData.N_v());

    // accumulate per vertex into local array
    vector<complex<double>> hq_local(m_FTData.N_q(),0);
    for (int v = 0; v < m_PosData.N_v(); ++v){
        double px = m_PosData.pos(r,f,v)[0];
        double py = m_PosData.pos(r,f,v)[1];
        double dh = m_PosData.pos(r,f,v)[2] - h_mean;
        for (int i = 0; i < m_FTData.N_q(); ++i){
            double phase = q[i][0]*px + q[i][1]*py;
            hq_local[i] += complex<double>(dh*std::cos(phase), dh*std::sin(phase));
        }
    }

    // write back contiguously
    for (int i = 0; i < m_FTData.N_q(); ++i){
        m_FTData.hq(r,f,i) = hq_local[i] / std::sqrt(m_PosData.N_v());
    }
}

inline void FourierTransform::LSFT(int r, int f){
    Eigen::MatrixXd A(m_PosData.N_v(), 2*m_FTData.N_q());
    Eigen::VectorXd dh(m_PosData.N_v());

    vector<array<double,2>> q(m_FTData.mn.size());
    for (int i = 0; i < m_FTData.N_q(); ++i){
        q[i][0] = 2*M_PI/m_PosData.boxsize(r,f)[0]*(double)(m_FTData.mn[i][0]);
        q[i][1] = 2*M_PI/m_PosData.boxsize(r,f)[1]*(double)(m_FTData.mn[i][1]); 
    }

    // precompute h_mean
    double h_mean = 0;
    for (int v = 0; v < m_PosData.N_v(); ++v)
        h_mean += m_PosData.pos(r,f,v)[2];
    h_mean /= (double)(m_PosData.N_v());

    // dh vector
    for (int v = 0; v < m_PosData.N_v(); ++v) {
        dh(v) = m_PosData.pos(r,f,v)[2] - h_mean;
    }

    // accumulate per vertex into local array
    for (int v = 0; v < m_PosData.N_v(); ++v){
        double px = m_PosData.pos(r,f,v)[0];
        double py = m_PosData.pos(r,f,v)[1];
        dh(v) = m_PosData.pos(r,f,v)[2] - h_mean;
        for (int i = 0; i < m_FTData.N_q(); ++i){
            double phase = q[i][0]*px + q[i][1]*py;
            A(v, 2*i)   = std::cos(phase)/std::sqrt(m_PosData.N_v());
            A(v, 2*i+1) = -std::sin(phase)/std::sqrt(m_PosData.N_v());
        }
    }
    //cout << A.jacobiSvd().singularValues() << "\n";
    //cout << "rep" << r+1 << " A0 dot A1 = " << A.col(0).dot(A.col(1)) << "\n";
    //cout << "rep" << r+1 << " A0 dot A5 = " << A.col(0).dot(A.col(5)) << "\n";

    // Solve min||A x - dh||_2   or min||A h_q - dh_v||_2
    Eigen::VectorXd x = A.colPivHouseholderQr().solve(dh);
    //Eigen::VectorXd x = A.householderQr().solve(dh);

    // write back contiguously
    for (int i = 0; i < m_FTData.N_q(); ++i){
        m_FTData.hq(r,f,i) = complex<double>(x(2*i), x(2*i+1))/(double)(2);
    }
}

void FourierTransform::RunFourierTransform(){
    //auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for (int r = 0; r < m_FTData.N_r(); r++){
        //auto start = std::chrono::high_resolution_clock::now();
        for(int f = 0; f< m_FTData.N_f(); f++){
            switch(m_FTMethod){
                case FTMethod::DFT: DFT(r,f); break;
                case FTMethod::LSFT:LSFT(r,f); break;
            }
        }
        //auto end = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed = end - start;
        //#pragma omp critical
        //cout << "Finished FourierTransform rep " << r+1 << " in " << elapsed.count() << " s\n"
        //     << "By thread " << omp_get_thread_num()+1 << " \n"; 
    }  
    //auto end = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double> elapsed = end - start;
    //cout << "Finished FourierTransform in " << elapsed.count() << " s\n";
}

#include <numeric> // for std::iota

namespace {
    // Groups indices [0, N) by a scalar key function, using tolerance-based
    // equality. Does NOT assume the data is pre-sorted by that key.
    template <typename KeyFunc>
    std::vector<std::vector<int>> groupByKey(int N, KeyFunc keyFunc, double tol = 1e-8) {
        std::vector<int> idx(N);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(),
                  [&](int a, int b) { return keyFunc(a) < keyFunc(b); });

        std::vector<std::vector<int>> groups;
        int i = 0;
        while (i < N) {
            double key_i = keyFunc(idx[i]);
            int j = i;
            std::vector<int> group;
            while (j < N && std::abs(keyFunc(idx[j]) - key_i) < tol) {
                group.push_back(idx[j]);
                ++j;
            }
            groups.push_back(std::move(group));
            i = j;
        }
        return groups;
    }
}

void FourierTransform::StaticFluctuation(){
    int N = m_FTData.N_q();
    // Window down to the currently-requested (init_f..final_f, init_rep..final_rep) sub-range of the
    // already-computed hq array (which spans the full loaded PosData) -- no DFT/LSFT recomputation.
    int N_r_active = m_Config.final_rep - m_Config.init_rep + 1;
    int N_f_active = m_Config.final_f - m_Config.init_f + 1;

    // Independent groupings for each quantity
    auto groups_q  = groupByKey(N, [&](int i){ return norm(m_FTData.q_mean[i]); });
    auto groups_qx = groupByKey(N, [&](int i){ return abs(m_FTData.q_mean[i][0]); });
    auto groups_qy = groupByKey(N, [&](int i){ return abs(m_FTData.q_mean[i][1]); });

    int N_unique_q  = (int)groups_q.size();
    int N_unique_qx = (int)groups_qx.size();
    int N_unique_qy = (int)groups_qy.size();

    NDArray2D<double> hquu_norm(N_r_active, N_unique_q,  0.0);
    NDArray2D<double> hqaa_norm(N_r_active, N_unique_q, 0.0);
    NDArray2D<double> hqbb_norm(N_r_active, N_unique_q, 0.0);
    NDArray2D<double> hqaa_qx_norm(N_r_active, N_unique_qx, 0.0);
    NDArray2D<double> hqbb_qy_norm(N_r_active, N_unique_qy, 0.0);

    // Per-replica q, q_x, q_y magnitudes (boxsize varies per replica/frame,
    // so these are frame-averaged within each replica, same as the observables below)
    NDArray2D<double> q_val(N_r_active, N_unique_q, 0.0);
    NDArray2D<double> qx_val(N_r_active, N_unique_qx, 0.0);
    NDArray2D<double> qy_val(N_r_active, N_unique_qy, 0.0);

    for (int r = 0; r < N_r_active; ++r) {
        for (int f = 0; f < N_f_active; ++f) {
            const array<double,3>& boxsize_rf = m_PosData.boxsize(r,f);

            // q grouped by |q|
            for (int k = 0; k < N_unique_q; ++k) {
                const auto& group = groups_q[k];
                double norm_fac = (double)(group.size() * N_f_active);
                for (int j : group) {
                    double qxj = 2*M_PI*(double)(m_FTData.mn[j][0]) / boxsize_rf[0];
                    double qyj = 2*M_PI*(double)(m_FTData.mn[j][1]) / boxsize_rf[1];
                    q_val(r,k) += std::sqrt(qxj*qxj + qyj*qyj) / norm_fac;
                }
            }

            // q_x grouped by |q_x|
            for (int k = 0; k < N_unique_qx; ++k) {
                const auto& group = groups_qx[k];
                double norm_fac = (double)(group.size() * N_f_active);
                for (int j : group) {
                    double qxj = 2*M_PI*(double)(m_FTData.mn[j][0]) / boxsize_rf[0];
                    qx_val(r,k) += std::abs(qxj) / norm_fac;
                }
            }

            // q_y grouped by |q_y|
            for (int k = 0; k < N_unique_qy; ++k) {
                const auto& group = groups_qy[k];
                double norm_fac = (double)(group.size() * N_f_active);
                for (int j : group) {
                    double qyj = 2*M_PI*(double)(m_FTData.mn[j][1]) / boxsize_rf[1];
                    qy_val(r,k) += std::abs(qyj) / norm_fac;
                }
            }

            // <u u> grouped by |q|
            for (int k = 0; k < N_unique_q; ++k) {
                const auto& group = groups_q[k];
                double norm_fac = (double)(group.size() * N_f_active);
                for (int j : group) {
                    complex hq = m_FTData.hq(r,f,j);
                    hquu_norm(r,k) += std::norm(hq) / norm_fac;
                }
            }

            // <a a> grouped by |q|
            for (int k = 0; k < N_unique_q; ++k) {
                const auto& group = groups_q[k];
                double norm_fac = (double)(group.size() * N_f_active);
                for (int j : group) {
                    complex hq = m_FTData.hq(r,f,j);
                    hqaa_norm(r,k) += std::real(hq)*std::real(hq) / norm_fac;
                }
            }

            // <b b> grouped by |q|
            for (int k = 0; k < N_unique_q; ++k) {
                const auto& group = groups_q[k];
                double norm_fac = (double)(group.size() * N_f_active);
                for (int j : group) {
                    complex hq = m_FTData.hq(r,f,j);
                    hqbb_norm(r,k) += std::imag(hq)*std::imag(hq) / norm_fac;
                }
            }

            // <a a>_qx grouped by q_x
            for (int k = 0; k < N_unique_qx; ++k) {
                const auto& group = groups_qx[k];
                double norm_fac = (double)(group.size() * N_f_active);
                for (int j : group) {
                    complex hq = m_FTData.hq(r,f,j);
                    hqaa_qx_norm(r,k) += std::real(hq)*std::real(hq) / norm_fac;
                }
            }

            // <b b>_qy grouped by q_y
            for (int k = 0; k < N_unique_qy; ++k) {
                const auto& group = groups_qy[k];
                double norm_fac = (double)(group.size() * N_f_active);
                for (int j : group) {
                    complex hq = m_FTData.hq(r,f,j);
                    hqbb_qy_norm(r,k) += std::imag(hq)*std::imag(hq) / norm_fac;
                }
            }
        }
    }

    std::string methodDir = m_Config.ftmethod + std::to_string(m_Config.N_k);
    std::string windowName = "init_f" + std::to_string(m_Config.init_f) + "final_f" + std::to_string(m_Config.final_f)
                            + "_init_rep" + std::to_string(m_Config.init_rep) + "final_rep" + std::to_string(m_Config.final_rep);
    if (!m_Config.align) windowName += "raw";
    std::filesystem::create_directories("StaticUndulationSpectrum/"+methodDir);
    std::ofstream outputFile("StaticUndulationSpectrum/"+methodDir+"/"+windowName+".txt");
    outputFile << "Mean Boxsize: \n";
    outputFile << "L_x L_y L_z \n";
    outputFile << m_PosData.boxsize_mean << "\n";
    outputFile << "Fourier Settings:\n";
    outputFile << "N_q_unique N_qx_unique N_qy_unique N_k FTMethod init_f final_f init_rep final_rep\n";
    outputFile << N_unique_q << " " << N_unique_qx << " " << N_unique_qy << " "
               << m_Config.N_k << " " << m_Config.ftmethod << " "
               << m_Config.init_f << " " << m_Config.final_f << " "
               << m_Config.init_rep << " " << m_Config.final_rep << "\n";

    auto writeBlock = [&](const std::string& name, int N_unique, const NDArray2D<double>& vals) {
        outputFile << name << "\n";
        outputFile << "rep:";
        for (int r = 0; r < N_r_active; r++)
            outputFile << " " << r+1;
        outputFile << "\n";
        for (int k = 0; k < N_unique; ++k) {
            for (int r = 0; r < N_r_active; ++r)
                outputFile << vals(r,k) << " ";
            outputFile << "\n";
        }
    };

    writeBlock("q", N_unique_q, q_val);
    writeBlock("q_x", N_unique_qx, qx_val);
    writeBlock("q_y", N_unique_qy, qy_val);
    writeBlock("<u u>", N_unique_q, hquu_norm);
    writeBlock("<u_x u_x>", N_unique_q, hqaa_norm);
    writeBlock("<u_y u_y>", N_unique_q, hqbb_norm);
    writeBlock("<u_x u_x>_qx", N_unique_qx, hqaa_qx_norm);
    writeBlock("<u_y u_y>_qy", N_unique_qy, hqbb_qy_norm);
}
/*
void FourierTransform::StaticFluctuation(){
    vector<int> idx_unique;
    idx_unique.push_back(0);
    int i = 0;
    int N_unique = 0;
    while (i < m_FTData.N_q()) {
        // Find end of degenerate group
        double mag_i = norm(m_FTData.q_mean[i]);
        int j = i;
        while (j < m_FTData.N_q() && std::abs(norm(m_FTData.q_mean[j]) - mag_i) < 1e-8)
            ++j;
        i=j;
        idx_unique.push_back(j);
        N_unique++;
    }
    
    NDArray2D<double> hquu_norm(m_FTData.N_r(),N_unique, 0.0);
    NDArray2D<double> hqaa_norm(m_FTData.N_r(),N_unique, 0.0);
    NDArray2D<double> hqbb_norm(m_FTData.N_r(),N_unique, 0.0);
    NDArray2D<double> hqab_norm(m_FTData.N_r(),N_unique, 0.0);
    for (int r = 0; r < m_FTData.N_r(); ++r)
        for (int f = 0; f < m_FTData.N_f(); ++f)
            for (int k = 0; k < N_unique; ++k)
                for (int j = idx_unique[k]; j < idx_unique[k+1]; ++j){
                    complex hq = m_FTData.hq(r,f,j);
                    hquu_norm(r,k) += std::norm(hq) / (double)((idx_unique[k+1]-idx_unique[k])*m_FTData.N_f());
                    hqaa_norm(r,k) += std::real(hq)*std::real(hq) / (double)((idx_unique[k+1]-idx_unique[k])*m_FTData.N_f());
                    hqbb_norm(r,k) += std::imag(hq)*std::imag(hq) / (double)((idx_unique[k+1]-idx_unique[k])*m_FTData.N_f());
                    hqab_norm(r,k) += std::imag(hq)*std::real(hq) / (double)((idx_unique[k+1]-idx_unique[k])*m_FTData.N_f());
                }
                    
    std::filesystem::create_directories("StaticUndulationSpectrum");
    std::ofstream outputFile("StaticUndulationSpectrum/"+m_Config.ftmethod+std::to_string(m_Config.N_k)+".txt");
    outputFile << "Mean Boxsize: \n";
    outputFile << "L_x L_y L_z \n";
    //outputFile << m_PosData.boxsize_mean << " " << m_PosData.boxsize_mean[1] << " " << boxsize_mean[2] << "\n";
    outputFile << m_PosData.boxsize_mean << "\n";
    outputFile << "Fourier Settings:\n";
    outputFile << "N_q_unique N_k FTMethod\n";
    outputFile << N_unique << " " << m_Config.N_k << " " << m_Config.ftmethod << "\n";
    
    outputFile << "q <u u> rep:";
    for (int r = 0; r < m_FTData.N_r();r++)
        outputFile << " " << r+1;
    outputFile << "\n";
    for (int k = 0; k < N_unique; ++k){
        int idx = idx_unique[k];
        double q_mag = norm(m_FTData.q_mean[idx]);
        outputFile << q_mag << " ";
        for (int r = 0; r < m_FTData.N_r(); ++r)
            outputFile << hquu_norm(r,k) << " ";
        outputFile << "\n";
    }

    outputFile << "q_m <u_x u_x> rep:";
    for (int r = 0; r < m_FTData.N_r();r++)
        outputFile << " " << r+1;
    outputFile << "\n";
    for (int k = 0; k < N_unique; ++k){
        int idx = idx_unique[k];
        double qx = m_FTData.q_mean[idx][0];
        outputFile << qx << " ";
        for (int r = 0; r < m_FTData.N_r(); ++r)
            outputFile << hqaa_norm(r,k) << " ";
        outputFile << "\n";
    }
    
    outputFile << "q_n <u_y u_y> rep:";
    for (int r = 0; r < m_FTData.N_r();r++)
        outputFile << " " << r+1;
    outputFile << "\n";
    for (int k = 0; k < N_unique; ++k){
        int idx = idx_unique[k];
        double qy = m_FTData.q_mean[idx][1];
        outputFile << qy << " ";
        for (int r = 0; r < m_FTData.N_r(); ++r)
            outputFile << hqbb_norm(r,k) << " ";
        outputFile << "\n";
    }

    outputFile << "q <a(k)b(k)> rep:";
    for (int r = 0; r < m_FTData.N_r();r++)
        outputFile << " " << r+1;
    outputFile << "\n";
    for (int k = 0; k < N_unique; ++k){
        int idx = idx_unique[k];
        double qx = m_FTData.q_mean[idx][0];
        double qy = m_FTData.q_mean[idx][1];
        double q = sqrt(qx*qx+qy*qy);
        outputFile << q << " ";
        for (int r = 0; r < m_FTData.N_r(); ++r)
            outputFile << hqab_norm(r,k) << " ";
        outputFile << "\n";
    }
}
*/
/*
void FourierTransform::StaticFluctuation(){
    vector<int> idx_unique;
    idx_unique.push_back(0);
    int i = 0;
    int N_unique = 0;
    while (i < m_FTData.N_q()) {
        // Find end of degenerate group
        double mag_i = norm(m_FTData.q_mean[i]);
        int j = i;
        while (j < m_FTData.N_q() && std::abs(norm(m_FTData.q_mean[j]) - mag_i) < 1e-8)
            ++j;
        i=j;
        idx_unique.push_back(j);
        N_unique++;
    }
    
    NDArray2D<double> hq_norm(m_FTData.N_r(),N_unique, 0.0);
    for (int r = 0; r < m_FTData.N_r(); ++r)
        for (int f = 0; f < m_FTData.N_f(); ++f)
            for (int k = 0; k < N_unique; ++k)
                for (int j = idx_unique[k]; j < idx_unique[k+1]; ++j)
                    hq_norm(r,k) += std::norm(m_FTData.hq(r,f,j)) / (double)((idx_unique[k+1]-idx_unique[k])*m_FTData.N_f());
         
    
    std::filesystem::create_directories("StaticUndulationSpectrum");
    std::ofstream outputFile("StaticUndulationSpectrum/"+m_Config.ftmethod+std::to_string(m_Config.N_k)+".txt");
    outputFile << "Mean Boxsize: \n";
    outputFile << "L_x L_y L_z \n";
    //outputFile << m_PosData.boxsize_mean << " " << m_PosData.boxsize_mean[1] << " " << boxsize_mean[2] << "\n";
    outputFile << m_PosData.boxsize_mean << "\n";
    outputFile << "Fourier Settings:\n";
    outputFile << "N_q_unique N_k FTMethod\n";
    outputFile << N_unique << " " << m_Config.N_k << " " << m_Config.ftmethod << "\n";
    
    outputFile << "q |h_q|^2 rep:";
    for (int r = 0; r < m_FTData.N_r();r++)
        outputFile << " " << r+1;
    outputFile << "\n";

    for (int k = 0; k < N_unique; ++k){
        int idx = idx_unique[k];
        double q_mag = norm(m_FTData.q_mean[idx]);
        outputFile << q_mag << " ";
        for (int r = 0; r < m_FTData.N_r(); ++r)
            outputFile << hq_norm(r,k) << " ";
        outputFile << "\n";
    }
}
*/

/*
void FourierTransform::StaticFluctuation(){
    vector<int> idx_unique;
    idx_unique.push_back(0);
    int i = 0;
    int N_unique = 0;
    while (i < m_FTData.N_q()) {
        // Find end of degenerate group
        double mag_i = norm(m_FTData.q_mean[i]);
        int j = i;
        while (j < m_FTData.N_q() && std::abs(norm(m_FTData.q_mean[j]) - mag_i) < 1e-8)
            ++j;
        i=j;
        idx_unique.push_back(j);
        N_unique++;
    }
    
    NDArray2D<double> hq_norm(m_FTData.N_r(),N_unique, 0.0);
    for (int r = 0; r < m_FTData.N_r(); ++r)
        for (int f = 0; f < m_FTData.N_f(); ++f)
            for (int k = 0; k < N_unique; ++k)
                for (int j = idx_unique[k]; j < idx_unique[k+1]; ++j)
                    hq_norm(r,k) += std::norm(m_FTData.hq(r,f,j)) / (double)((idx_unique[k+1]-idx_unique[k])*m_FTData.N_f());
         
    
    std::filesystem::create_directories("StaticUndulationSpectrum");
    std::ofstream outputFile("StaticUndulationSpectrum/"+m_Config.ftmethod+std::to_string(m_Config.N_k)+".txt");
    outputFile << "Mean Boxsize: \n";
    outputFile << "L_x L_y L_z \n";
    //outputFile << m_PosData.boxsize_mean << " " << m_PosData.boxsize_mean[1] << " " << boxsize_mean[2] << "\n";
    outputFile << m_PosData.boxsize_mean << "\n";
    outputFile << "Fourier Settings:\n";
    outputFile << "N_q_unique N_k FTMethod\n";
    outputFile << N_unique << " " << m_Config.N_k << " " << m_Config.ftmethod << "\n";
    
    outputFile << "q |h_q|^2 rep:";
    for (int r = 0; r < m_FTData.N_r();r++)
        outputFile << " " << r+1;
    outputFile << "\n";

    for (int k = 0; k < N_unique; ++k){
        int idx = idx_unique[k];
        double q_mag = norm(m_FTData.q_mean[idx]);
        outputFile << q_mag << " ";
        for (int r = 0; r < m_FTData.N_r(); ++r)
            outputFile << hq_norm(r,k) << " ";
        outputFile << "\n";
    }


    // Per-replica average over degenerate partners
    NDArray3D<double> hq_norm(m_FTData.N_r(),m_FTData.N_f(),N_unique, 0.0);
    for (int r = 0; r < m_FTData.N_r(); ++r)
        for (int f = 0; f < m_FTData.N_f(); ++f)
            for (int k = 0; k < N_unique; ++k)
                for (int j = idx_unique[k]; j < idx_unique[k+1]; ++j)
                    hq_norm(r,f,k) += std::norm(m_FTData.hq(r,f,j)) / (double)(idx_unique[k+1]-idx_unique[k]);
    
    

    
    NDArray2D<double> hq_fmean(m_FTData.N_r(),N_unique,0.0);
    for (int r = 0; r < m_FTData.N_r(); ++r)
        for (int f = 0; f < m_FTData.N_f(); ++f)
            for (int k = 0; k < N_unique; ++k)
                hq_fmean(r,k) += hq_norm(r,f,k) / (double)(m_FTData.N_f());

    vector<double> s_mean(N_unique,0.0);
    for (int r = 0; r < m_FTData.N_r(); ++r)
        for (int k = 0; k < N_unique; ++k)
            s_mean[k] += hq_fmean(r,k)/(double)(m_FTData.N_r());
    
    vector<double> var_replica(N_unique,0.0);
    for (int r = 0; r < m_FTData.N_r(); ++r)
        for (int k = 0; k < N_unique; ++k)
            var_replica[k] += (hq_fmean(r,k)-s_mean[k])*(hq_fmean(r,k)-s_mean[k]) / (double)(m_FTData.N_r()-1);
    for (int k = 0; k < N_unique; ++k)
        var_replica[k] /= (double)(m_FTData.N_r());
    
    NDArray2D<double> ACF0(m_FTData.N_r(),N_unique,0.0);
    for (int r = 0; r < m_FTData.N_r(); ++r)
        for (int f = 0; f < m_FTData.N_f(); ++f)
            for (int k = 0; k < N_unique; ++k)
                ACF0(r,k) += ((hq_norm(r,f,k) - hq_fmean(r,k)) * (hq_norm(r,f,k) - hq_fmean(r,k)));

    NDArray3D<double> ACF(m_FTData.N_r(),m_N_dt_max,N_unique,0.0);
    for (int r = 0; r < m_FTData.N_r(); ++r)
        for (int dt = 0; dt < m_N_dt_max; ++dt)
            for (int f = 0; f < m_FTData.N_f() - dt; ++f)
                for (int k = 0; k < N_unique; ++k)
                    ACF(r,dt,k) += (hq_norm(r,f+dt,k) - hq_fmean(r,k)) * (hq_norm(r,f,k) - hq_fmean(r,k)) / ACF0(r,k);
    
    NDArray2D<double> tau(m_FTData.N_r(),N_unique,0.0);
    for (int r = 0; r < m_FTData.N_r(); ++r){
        for (int k = 0; k < N_unique; ++k){
            tau(r,k) = 0.5;
            for (int dt = 1; dt < m_N_dt_max; ++dt){
                if (ACF(r,dt,k) < 0.05) break;
                tau(r,k) += 2*ACF(r,dt,k);
            }
        }
    }

    vector<double> var_frames (N_unique,0.0);
    for (int r = 0; r < m_FTData.N_r(); ++r){
        vector<double> var_r(N_unique, 0.0);
        for (int f = 0; f < m_FTData.N_f(); ++f)
            for (int k = 0; k < N_unique; ++k)
                var_r[k] += (hq_norm(r,f,k)-hq_fmean(r,k)) *
                            (hq_norm(r,f,k)-hq_fmean(r,k)) / (double)(m_FTData.N_f()-1);
        for (int k = 0; k < N_unique; ++k){
            double N_eff = (double)(m_FTData.N_f()) / tau(r,k);
            var_frames[k] += var_r[k] / (N_eff * (double)(m_FTData.N_r()));
        }
    }
    
    std::filesystem::create_directories("StaticUndulationSpectrum");
    std::ofstream outputFile("StaticUndulationSpectrum/"+m_Config.ftmethod+std::to_string(m_Config.N_k)+".txt");
    outputFile << "Mean Boxsize: \n";
    outputFile << "L_x L_y L_z \n";
    outputFile << m_PosData.boxsize[0] << " " << m_PosData.boxsize[1] << " " << m_PosData.boxsize[2] << "\n";
    outputFile << "Fourier Settings:\n";
    outputFile << "N_q_unique N_k FTMethod\n";
    outputFile << N_unique << " " << m_Config.N_k << " " << m_Config.ftmethod << "\n";
    outputFile << "q |h_q|^2 var: replica frames total\n";
    outputFile << "q |h_q|^2 var: replica frames total\n";
    for (int k = 0; k < N_unique; ++k){
        cout << "var_replica = " << var_replica[k] << "\n var_frames = " << var_frames[k] << "\n";
        double var_total = var_replica[k] + var_frames[k];
        outputFile << q_magnitude(idx_unique[k]) << " " << s_mean[k] << " " << var_replica[k] << " " << var_frames[k] << " " << var_total << "\n";
    }
}
*/


void FourierTransform::AutoCorrelation() {
    // Window down to the currently-requested (init_f..final_f, init_rep..final_rep) sub-range of the
    // already-computed hq array (which spans the full loaded PosData) -- no DFT/LSFT recomputation.
    // tau_max still caps the max lag exactly as before, just measured against the active N_f.
    int N_r_active = m_Config.final_rep - m_Config.init_rep + 1;
    int N_f_active = m_Config.final_f - m_Config.init_f + 1;
    int N_dt_max = (m_Config.tau_max > 0) ? std::min(m_Config.tau_max, N_f_active) : N_f_active;

    std::string methodDir = m_Config.ftmethod + std::to_string(m_Config.N_k);
    std::string windowName = "init_f" + std::to_string(m_Config.init_f) + "final_f" + std::to_string(m_Config.final_f)
                            + "_init_rep" + std::to_string(m_Config.init_rep) + "final_rep" + std::to_string(m_Config.final_rep);
    if (!m_Config.align) windowName += "raw";
    std::string outDir = "AutoCorrelationFunction/"+methodDir+"/"+windowName;
    std::filesystem::create_directories(outDir);
    //hq = (N_r,N_f,N_q)
    //ACF = (N_q,N_r,N_t) : <|h_q|^2(t)|h_q|^2(0)> averaged over time origins

    m_tau.clear(); // this call's tau(i) per q_i; cleared so repeated calls (one per final_f/final_rep) don't accumulate stale entries

    // time-averaged <|h_q|^2> per (mode,replica), needed to subtract the disconnected
    // part so the correlation function actually decays to 0 for non-OU processes
    NDArray2D<double> X_mean(m_FTData.N_q(), N_r_active, 0.0);
    for (int i = 0; i < m_FTData.N_q(); ++i)
        for (int r = 0; r < N_r_active; ++r)
            for (int f = 0; f < N_f_active; ++f)
                X_mean(i,r) += std::norm(m_FTData.hq(r,f,i)) / (double)(N_f_active);

    NDArray3D<double> ACF(m_FTData.N_q(),N_r_active,N_dt_max,0);
    for (int i = 0; i < m_FTData.N_q(); ++i)
        for (int r=0; r < N_r_active; ++r)
            for (int dt = 0; dt < N_dt_max; ++dt)
                for (int f = 0; f < N_f_active - dt; ++f)
                    ACF(i,r,dt) +=  std::norm(m_FTData.hq(r,f+dt,i))*std::norm(m_FTData.hq(r,f,i)) / (double)(N_f_active-dt);

    // Subtract the disconnected part <X>^2 (X(t) and X(0) become independent as dt
    // grows, so <X(t)X(0)> -> <X>^2, not 0) and normalize by the connected variance.
    for (int i = 0; i < m_FTData.N_q(); ++i)
        for (int r = 0; r < N_r_active; ++r) {
            double X2_mean = X_mean(i,r) * X_mean(i,r);
            double var0 = ACF(i,r,0) - X2_mean;
            for (int dt = 1; dt < N_dt_max; ++dt)
                ACF(i,r,dt) = (ACF(i,r,dt) - X2_mean) / var0;
            ACF(i,r,0) = 1;
        }

    int i = 0;
    int q_index = 0; // ascending group index (q0, q1, q2, ...), unlike i which skips over degenerate partners
    while (i < m_FTData.N_q()) {
        // Find end of degenerate group
        double mag_i = norm(m_FTData.q_mean[i]);
        int j = i;
        while (j < m_FTData.N_q() && std::abs(norm(m_FTData.q_mean[j]) - mag_i) < 1e-10)
            ++j;

        // modes i..j-1 are degenerate
        int N_degenerate = j - i;

        // Per-replica average over degenerate partners
        NDArray2D<complex<double>> hq_group(N_r_active,N_dt_max,0);
        for (int r = 0; r < N_r_active; ++r)
            for(int dt = 0; dt < N_dt_max; ++dt)
                for (int k = i; k < j; ++k)
                    hq_group(r,dt) += ACF(k,r,dt) / (double)(N_degenerate);

        vector<complex<double>> a_mean(N_dt_max,0);
        for (int r = 0; r < N_r_active; ++r)
            for (int dt = 0; dt < N_dt_max; ++dt)
                a_mean[dt] += hq_group(r,dt) / (double)(N_r_active);

        double tau = 0.5; // integrated correlation time, only accumulated up to where the ACF decorrelates
        for (int dt = 1; dt < N_dt_max; dt++){
            if (std::abs(a_mean[dt]) <= 0.1) break;
            tau += std::abs(a_mean[dt]);
        }
        m_tau.push_back(tau);

        // full ACF (including the noise-dominated tail past decorrelation) is still
        // written out in full, so downstream plotting can choose its own [dt1,dt2] window
        vector<double> a_var(N_dt_max,0);
        for (int r = 0; r < N_r_active; ++r)
            for (int dt = 0; dt < N_dt_max; ++dt)
                a_var[dt] += std::norm(hq_group(r,dt)-a_mean[dt]) / (double)(N_r_active-1);

        std::ofstream outputFile(outDir+"/q" + std::to_string(q_index) +".txt");
        outputFile << "q"+std::to_string(q_index) +" tau N_r\n";
        outputFile << mag_i << " " << tau << " " << N_r_active << "\n";
        outputFile << "t mean var\n";
        for (int dt = 0; dt < N_dt_max; ++dt){
            outputFile << dt << " " << std::abs(a_mean[dt]) << " " << a_var[dt] << "\n";
        }

        i = j; // advance to next group
        ++q_index;
    }
}


