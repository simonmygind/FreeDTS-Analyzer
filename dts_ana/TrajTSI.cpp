#include "TrajTSI.h"

#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

using std::cout;

PosData TrajTSI::LoadConfig(){
    int N_f, N_r, N_v;
    N_f = m_Config.final_f-m_Config.init_f+1; 
    N_r = m_Config.final_rep-m_Config.init_rep+1;

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
    cout << "TSI files of N_v "<< N_v << " vertices" << "\n";
    PosData PosData(N_r, N_f, N_v);
    return PosData;
}

void TrajTSI::LoadTSI(){
    // boxsize = curent boxsize
    // m_boxsize = mean boxsize
    array<double,3> boxsize;
    std::string token;
    cout << "Start Loading TSI files...\n";
    double dummy;
    for (int r = 0; r < m_PosData.N_r(); ++r) {
        for (int f = 0; f < m_PosData.N_f(); ++f){
            std::ifstream inputFile(m_Config.base_dir+"rep"+std::to_string(r+m_Config.init_rep)+"/TrajTSI/dts"+std::to_string(f+m_Config.init_f)+".tsi");
            if (!inputFile) {
                std::cerr << "Error: cannot open frame " << f+m_Config.init_f << " in rep " << r+m_Config.init_rep <<"\n";
                std::cerr << "or file " + m_Config.base_dir+"rep"+std::to_string(r+m_Config.init_rep)+"/TrajTSI/dts"+std::to_string(f+m_Config.init_f)+".tsi \n";
                std::cerr << "Specify -base_dir *__*/ \n";
                exit(1);
            }
            
            while (inputFile >> token) {
                if (token == "box"){
                    inputFile >> boxsize[0] >> boxsize[1] >> boxsize[2];
                    m_PosData.boxsize[0] += boxsize[0]/(double)(m_PosData.N_f() *m_PosData.N_r());  
                    m_PosData.boxsize[1] += boxsize[1]/(double)(m_PosData.N_f() *m_PosData.N_r());  
                    m_PosData.boxsize[2] += boxsize[2]/(double)(m_PosData.N_f() *m_PosData.N_r());         
                }
                else if (token == std::to_string(m_PosData.N_v())){
                    for (int v = 0; v < m_PosData.N_v(); ++v) {
                        inputFile   >> dummy 
                                    >> m_PosData.pos.At({r,f,v})[0]
                                    >> m_PosData.pos.At({r,f,v})[1] 
                                    >> m_PosData.pos.At({r,f,v})[2];
                    }
                    break;
                }
                
            }
        }
        cout << "Finished loading rep " << r+1 << "\n";
    }
}