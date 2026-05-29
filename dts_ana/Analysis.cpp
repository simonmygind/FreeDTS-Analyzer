
#include "Analysis.h"
#include <fstream>
#include <chrono>

Config Analysis::ParseArgs(const int argc,const char* argv[])
{
    Config Config;
    vector<std::string> analysis;
    for(int i = 1; i < argc; i++){
        //cout << argv[i] << "\n";
        if (std::string(argv[i]) == "-init_f"){
            Config.init_f = atoi(argv[i+1]);
        }
        else if (std::string(argv[i]) == "-final_f"){
            Config.final_f= atoi(argv[i+1]);
        }
        else if (std::string(argv[i]) == "-init_rep"){
            Config.init_rep= atoi(argv[i+1]);
        }
        else if (std::string(argv[i]) == "-final_rep"){
            Config.final_rep= atoi(argv[i+1]);
        }
        else if (std::string(argv[i]) == "-N_k"){
            Config.N_k = atoi(argv[i+1]);
        }
        else if (std::string(argv[i]) == "-N_t"){
            Config.N_t = atoi(argv[i+1]);
        }
        else if (std::string(argv[i]) == "-base_dir"){
            Config.base_dir = std::string(argv[i+1]);
        }
        
        else if (std::string(argv[i]) == "-analysis"){
            int j = 1;
            while (i+j < argc){
                if      (std::string(argv[i+j]) == "SF" ){
                    Config.analysis.push_back(std::string(argv[i+j]));
                }
                else if (std::string(argv[i+j]) == "AC" ){
                    Config.analysis.push_back(std::string(argv[i+j]));
                }
                else if (std::string(argv[i+j]) == "DF"){
                    Config.analysis.push_back(std::string(argv[i+j]));
                }
                else {
                    break;
                }
                j+= 1;
            }      
        }
        else if (std::string(argv[i]) == "-ftmethod"){
            if (argv[i+1] == "DFT"){
                m_Config.ftmethod = "DFT";
            }
            else{
                break;
            }
        }
    }

    if (Config.final_f == 0){
        cout << "Specify -final_f" << std::endl;
        exit(1); 
    }
    
    if (Config.final_f < Config.init_f){
        cout << "rep_init is bigger than rep_final" << std::endl;
        exit(1);
    }
    if (Config.analysis.size() == 0){
        cout << "Specify -analysis \n";
        exit(1);
    }
    
    return Config;
}

void Analysis::RunAnalysis(){
    cout << "Running analysis:";
    for (std::string& a : m_Config.analysis){
        cout << " " << a;
    }
    cout << "\n";


    cout << "From directory: " << m_Config.base_dir << "\n";

    cout << "frame = " << m_Config.init_f <<" to " << m_Config.final_f << "\n";
    cout << "rep = " << m_Config.init_rep <<" to " << m_Config.final_rep << "\n";
    cout << "----------------------------------------------------------\n";
    
    bool construct_TrajTSI = 0, construct_FourierTransform = 0;
    bool run_StaticFluctuation = 0, run_AutoCorrelation = 0;
    for (auto& a : m_Config.analysis){
        
        if (a == "SF") {
            construct_TrajTSI = 1; 
            construct_FourierTransform=1;
            run_StaticFluctuation = 1;
        }
        if (a == "AC"){
            construct_TrajTSI = 1; 
            construct_FourierTransform=1;
            run_AutoCorrelation = 1;
        }
        if (a == "DF"){
            construct_TrajTSI =1;
        }
    }

    
    // Contional construction of TrajTSI 
    if (construct_TrajTSI) {
        
        m_TrajTSI = std::make_unique<TrajTSI>(m_Config);
        auto start = std::chrono::high_resolution_clock::now();
        m_TrajTSI->LoadTSI();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        cout << "Finished Loading TSI Files in " << elapsed.count() << " s \n";
    }
    
    /*
    cout << "PosData N_r: " << m_TrajTSI->getPosData().N_r() << "\n";
    cout << "PosData N_f: " << m_TrajTSI->getPosData().N_f() << "\n";
    cout << "PosData N_v: " << m_TrajTSI->getPosData().N_v() << "\n";
    cout << "boxsize: " << m_TrajTSI->getPosData().boxsize[0] << "\n";
    */

    // Conditional construction of FourierTransform
    if (construct_FourierTransform){
        cout << "Starting Fourier Transform...\n";
        m_FourierTransform = std::make_unique<FourierTransform>(m_TrajTSI->getPosData(),m_Config.N_k,m_Config.ftmethod);
        auto start = std::chrono::high_resolution_clock::now();
        m_FourierTransform->RunFourierTransform();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        cout << "Finished Fourier Transform in " << elapsed.count() << " s \n";

        if(run_StaticFluctuation == 1){
            auto start = std::chrono::high_resolution_clock::now();
            m_FourierTransform->StaticFluctuation();
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            cout << "Finished Static Fluctuation in " << elapsed.count() << " s \n";
        }
    }
}