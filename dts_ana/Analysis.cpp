
#include "Analysis.h"
#include <fstream>
#include <chrono>
#include <algorithm>

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
            Config.final_f_list.clear();
            int j = 1;
            while (i+j < argc && argv[i+j][0] != '-'){
                Config.final_f_list.push_back(atoi(argv[i+j]));
                j += 1;
            }
        }
        else if (std::string(argv[i]) == "-init_rep"){
            Config.init_rep= atoi(argv[i+1]);
        }
        else if (std::string(argv[i]) == "-final_rep"){
            Config.final_rep_list.clear();
            int j = 1;
            while (i+j < argc && argv[i+j][0] != '-'){
                Config.final_rep_list.push_back(atoi(argv[i+j]));
                j += 1;
            }
        }
        else if (std::string(argv[i]) == "-N_k"){
            Config.N_k_list.clear();
            int j = 1;
            while (i+j < argc && argv[i+j][0] != '-'){
                Config.N_k_list.push_back(atoi(argv[i+j]));
                j += 1;
            }
        }
        else if (std::string(argv[i]) == "-tau_max"){
            Config.tau_max = atoi(argv[i+1]);
        }
        else if (std::string(argv[i]) == "-base_dir"){
            Config.base_dir = std::string(argv[i+1]);
        }
        else if (std::string(argv[i]) == "-raw"){
            Config.align = false;
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
            Config.ftmethod_list.clear();
            int j = 1;
            while (i+j < argc){
                std::string tok = argv[i+j];
                if (tok == "DFT" || tok == "LSFT"){
                    Config.ftmethod_list.push_back(tok);
                }
                else break;
                j += 1;
            }
            if (Config.ftmethod_list.empty()){
                std::cerr << "Specify Fourier Transform Method(s)\n";
                std::cerr << "Use DFT and/or LSFT\n";
                exit(1);
            }
        }
    }

    if (Config.final_f_list.empty()){
        cout << "Specify at least one -final_f value\n";
        exit(1);
    }
    if (Config.final_rep_list.empty()){
        cout << "Specify at least one -final_rep value\n";
        exit(1);
    }
    if (*std::min_element(Config.final_f_list.begin(), Config.final_f_list.end()) < Config.init_f){
        cout << "init_f is bigger than final_f" << std::endl;
        exit(1);
    }
    if (*std::min_element(Config.final_rep_list.begin(), Config.final_rep_list.end()) < Config.init_rep){
        cout << "init_rep is bigger than final_rep" << std::endl;
        exit(1);
    }
    if (Config.analysis.size() == 0){
        cout << "Specify -analysis \n";
        exit(1);
    }
    if (Config.N_k_list.empty()){
        cout << "Specify at least one -N_k value\n";
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

    cout << "frame = " << m_Config.init_f << " to "
         << *std::max_element(m_Config.final_f_list.begin(), m_Config.final_f_list.end())
         << " (final_f sweep:";
    for (int f : m_Config.final_f_list) cout << " " << f;
    cout << ")\n";
    cout << "rep = " << m_Config.init_rep << " to "
         << *std::max_element(m_Config.final_rep_list.begin(), m_Config.final_rep_list.end())
         << " (final_rep sweep:";
    for (int r : m_Config.final_rep_list) cout << " " << r;
    cout << ")\n";
    cout << "----------------------------------------------------------\n";
    
    bool construct_TrajTSI = 0, construct_FourierTransform = 0;
    bool run_StaticFluctuation = 0, run_AutoCorrelation = 0;
    for (auto& a : m_Config.analysis){
        
        if (a == "SF") {
            construct_TrajTSI = 1; 
            construct_FourierTransform=1;
            run_StaticFluctuation = 1;
            run_AutoCorrelation = 1;
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

    
    // Conditional construction of TrajTSI 
    if (construct_TrajTSI) {
        
        m_TrajTSI = std::make_unique<TrajTSI>(m_Config);
        auto start = std::chrono::high_resolution_clock::now();
        m_TrajTSI->LoadTSI();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        cout << "Finished Loading TSI Files in " << elapsed.count() << " s \n";

        if (m_Config.align){
            m_TrajTSI->AlignTilt();
        }
        else{
            cout << "Skipping tilt alignment (-raw): using unrotated lab-frame coordinates\n";
        }
    }

    cout << "----------------------------------------------------------\n";
    // Conditional construction of FourierTransform: sweeps every (ftmethod, N_k) combination against
    // the same already-loaded (and already tilt-aligned) PosData -- loaded out to
    // max(final_rep_list)/max(final_f_list) -- so TSI files are only ever read from disk once, and
    // the DFT/LSFT itself only computed once per (ftmethod, N_k), regardless of how many
    // final_f/final_rep amounts-of-data are requested. Each (final_rep, final_f) in the inner sweep
    // just windows StaticFluctuation/AutoCorrelation down to a smaller rep/frame prefix of that same
    // already-computed hq array (see FourierTransform::StaticFluctuation / ::AutoCorrelation).
    if (construct_FourierTransform){
        for (const std::string& method : m_Config.ftmethod_list){
            for (int N_k : m_Config.N_k_list){
                m_Config.ftmethod = method;
                m_Config.N_k = N_k;

                cout << "\n";
                cout << "Running FTMethod = " << m_Config.ftmethod << ", N_k = " << m_Config.N_k << " ....\n";
                m_FourierTransform = std::make_unique<FourierTransform>(m_Config,m_TrajTSI->getPosData());
                cout << "With N_q = " << m_FourierTransform->getFTData().N_q() << "\n";
                auto start = std::chrono::high_resolution_clock::now();
                m_FourierTransform->RunFourierTransform();
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                cout << "Finished Fourier Transform in " << elapsed.count() << " s \n";

                for (int final_rep : m_Config.final_rep_list){
                    for (int final_f : m_Config.final_f_list){
                        m_Config.final_rep = final_rep;
                        m_Config.final_f = final_f;

                        cout << "  Windowing to final_rep = " << final_rep << ", final_f = " << final_f << " ...\n";

                        if(run_StaticFluctuation == 1){
                            start = std::chrono::high_resolution_clock::now();
                            m_FourierTransform->StaticFluctuation();
                            end = std::chrono::high_resolution_clock::now();
                            elapsed = end - start;
                            cout << "  Finished Static Fluctuation in " << elapsed.count() << " s \n";
                        }
                        if(run_AutoCorrelation == 1){
                            start = std::chrono::high_resolution_clock::now();
                            m_FourierTransform->AutoCorrelation();
                            end = std::chrono::high_resolution_clock::now();
                            elapsed = end - start;
                            cout << "  Finished Autocorrelation Function in " << elapsed.count() << " s \n";
                        }
                    }
                }
                cout << "----------------------------------------------------------\n";
            }
        }
    }
}