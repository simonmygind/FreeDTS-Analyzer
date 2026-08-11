#include "TrajTSI.h"
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <algorithm>
#include <Eigen/Dense>
#include "omp.h"
#include "VecMath.h"

using std::cout;

static double next_double(char*& ptr){
    char* end;
    double x = std::strtod(ptr, &end);
    ptr = end;
    return x;
}

PosData TrajTSI::LoadConfig(){
    int N_f, N_r, N_v = -1;
    // Load out to the largest requested final_f/final_rep so every (final_f, final_rep) combo in
    // the sweep can window down into this same in-memory trajectory without re-reading TSI files.
    int max_final_f = *std::max_element(m_Config.final_f_list.begin(), m_Config.final_f_list.end());
    int max_final_rep = *std::max_element(m_Config.final_rep_list.begin(), m_Config.final_rep_list.end());
    N_f = max_final_f-m_Config.init_f+1;
    N_r = max_final_rep-m_Config.init_rep+1;

    std::ifstream inputFile(m_Config.base_dir+"rep"+std::to_string(m_Config.init_rep)+"/TrajTSI/dts"+std::to_string(m_Config.init_f)+".tsi");
    if (!inputFile) {
        std::cerr << "Error: cannot open initial frame " << m_Config.init_f << " in rep " << m_Config.init_rep <<"\n";
        std::cerr << "or file " + m_Config.base_dir+"rep"+std::to_string(m_Config.init_rep)+"/TrajTSI/dts"+std::to_string(m_Config.init_f)+".tsi \n";
        std::cout << "Specify -base_dir _/";
        exit(1);
    }
    std::string token;
    while (inputFile >> token) {
        if (token == "vertex"){
            inputFile >> N_v;
            continue;
        }
    }
    if (N_v == -1){
        std::cerr << "N_v was not read in corretly";
        exit(1);
    }

    cout << "TSI files of N_v "<< N_v << " vertices" << "\n";
    PosData PosData(N_r, N_f, N_v);
    return PosData;
}


void TrajTSI::LoadTSI(){
    // Invariant variables: 
    const int N_r  = m_PosData.N_r();
    const int N_f  = m_PosData.N_f();
    const int N_v  = m_PosData.N_v();
    const double inv_norm = 1/(double)(N_f*N_r);
    

    cout << "Start Loading TSI files...\n";    
    auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for (int r = 0; r < N_r; ++r) {
        //    auto start = std::chrono::high_resolution_clock::now();
        for (int f = 0; f < N_f; ++f){
            std::string path = m_Config.base_dir+"rep"+std::to_string(r+m_Config.init_rep)+"/TrajTSI/dts"+std::to_string(f+m_Config.init_f)+".tsi";
            FILE* fp = std::fopen(path.c_str(), "r");

            if (!fp){
                std::fprintf(stderr,"Error Cannot open %s\n",path.c_str());
                std::exit(1);
            }

            char buf[4096];
            std::fgets(buf,sizeof(buf),fp);
            std::fgets(buf,sizeof(buf),fp);
            char* ptr = buf + 4; // Skips past 'box'
            array<double,3>& b = m_PosData.boxsize(r,f);
            b[0] = next_double(ptr);
            b[1] = next_double(ptr);
            b[2] = next_double(ptr);
            
            std::fgets(buf, sizeof(buf), fp);
            for (int v = 0; v < N_v; ++v) {
                std::fgets(buf, sizeof(buf), fp);
                char* ptr = buf;
                next_double(ptr); //Skips past vertex id 
                array<double,3>& p = m_PosData.pos(r, f, v);
                p[0]    = next_double(ptr);
                p[1]    = next_double(ptr);
                p[2]    = next_double(ptr);
            }
            std::fclose(fp);
        }
        //#pragma omp critical
        //cout << "Finished loading rep " << r+1 << " in " << elapsed.count() <<" s by thread "<< omp_get_thread_num()+1 << " \n";
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << "Finished loading in " << elapsed.count() <<" s \n";

    for (int r = 0; r < m_PosData.N_r(); ++r)
        for (int f = 0; f < m_PosData.N_f(); ++f)
            m_PosData.boxsize_mean += m_PosData.boxsize(r,f)*inv_norm;

}

// Removes static/systematic membrane tilt before Fourier analysis: per replica,
// fits the best-fit plane to the ENTIRE trajectory's vertex point cloud (pooled
// over all frames) and applies the single minimal rotation that brings that
// time-averaged normal onto z, identically to every frame of that replica.
//
// This must be one fixed rotation per replica, not one per frame: an
// instantaneous per-frame best-fit normal is dominated by whatever long-
// wavelength bending fluctuation that frame happens to be in the middle of, so
// re-fitting and re-rotating every frame independently cancels out exactly the
// low-q undulation signal the spectrum is meant to measure, and scrambles which
// physical direction a given (m,n) mode index corresponds to from frame to
// frame -- averaging |h_q|^2 over frames that no longer share a common in-plane
// basis. A single time-averaged rotation per replica corrects only the static
// bias (e.g. from box geometry or initial condition) while leaving genuine
// thermal fluctuations relative to that fixed reference plane untouched.
void TrajTSI::AlignTilt(){
    const int N_r = m_PosData.N_r();
    const int N_f = m_PosData.N_f();
    const int N_v = m_PosData.N_v();

    cout << "Aligning membrane tilt (rotating each replica's time-averaged mean normal onto z)...\n";
    auto start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for (int r = 0; r < N_r; ++r){
        // Centroid pooled over every vertex of every frame in this replica
        Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
        for (int f = 0; f < N_f; ++f)
            for (int v = 0; v < N_v; ++v){
                const array<double,3>& p = m_PosData.pos(r,f,v);
                centroid += Eigen::Vector3d(p[0],p[1],p[2]);
            }
        centroid /= (double)(N_f) * (double)(N_v);

        // Time-averaged covariance (gyration) tensor about that centroid
        Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
        for (int f = 0; f < N_f; ++f)
            for (int v = 0; v < N_v; ++v){
                const array<double,3>& p = m_PosData.pos(r,f,v);
                Eigen::Vector3d d(p[0]-centroid[0], p[1]-centroid[1], p[2]-centroid[2]);
                cov += d * d.transpose();
            }
        cov /= (double)(N_f) * (double)(N_v);

        // Mean normal = eigenvector of smallest eigenvalue (least-variance direction
        // of a quasi-planar point cloud). computeDirect() uses Eigen's closed-form
        // trigonometric solver for fixed 3x3 matrices instead of the general iterative
        // tridiagonalize+QL path that the plain constructor/compute() would take.
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver;
        solver.computeDirect(cov);
        Eigen::Vector3d n = solver.eigenvectors().col(0); // ascending eigenvalues
        if (n[2] < 0) n = -n; // fix sign convention: normal points toward +z

        // Minimal rotation taking n onto z_hat (Rodrigues' rotation formula)
        Eigen::Vector3d z_hat(0.0,0.0,1.0);
        Eigen::Vector3d axis = n.cross(z_hat);
        double s = axis.norm();
        double c = n.dot(z_hat);

        Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
        if (s > 1e-12){
            Eigen::Matrix3d K;
            K <<      0, -axis[2],  axis[1],
                 axis[2],       0, -axis[0],
                -axis[1], axis[0],       0;
            R = Eigen::Matrix3d::Identity() + K + K*K*((1.0-c)/(s*s));
        }

        // Apply the SAME replica-level rotation to every frame, about the same
        // replica-level centroid (translation only shifts the Fourier phase,
        // which the downstream power-spectrum analysis doesn't use).
        for (int f = 0; f < N_f; ++f)
            for (int v = 0; v < N_v; ++v){
                array<double,3>& p = m_PosData.pos(r,f,v);
                Eigen::Vector3d d(p[0]-centroid[0], p[1]-centroid[1], p[2]-centroid[2]);
                Eigen::Vector3d rp = R * d;
                p[0] = rp[0] + centroid[0];
                p[1] = rp[1] + centroid[1];
                p[2] = rp[2] + centroid[2];
            }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << "Finished tilt alignment in " << elapsed.count() << " s \n";
}
