#pragma once
#include <vector>
#include <fstream>

using std::vector;
struct Config // Contains metadata that can be changed/ set in the terminal
{   
    int init_f = 1;
    int final_f = 1; // final_f currently being run (populated from final_f_list during the sweep)
    vector<int> final_f_list = {1}; // sweep of final_f values to run (populated by -final_f), reusing the same loaded PosData (loaded out to max(final_f_list))
    int init_rep = 1;
    int final_rep = 1; // final_rep currently being run (populated from final_rep_list during the sweep)
    vector<int> final_rep_list = {1}; // sweep of final_rep values to run (populated by -final_rep), reusing the same loaded PosData (loaded out to max(final_rep_list))
    int N_k = 5; // N_k of the FourierTransform currently being run
    vector<int> N_k_list = {5}; // sweep of N_k values to run (populated by -N_k), reusing the same loaded PosData
    int tau_max = -1; // max lag in frames to autocorrelate out to; <=0 means use the full N_f range (the ACF<0.1 break condition usually stops well before that)
    bool align = true; // Rotate each frame so the mean membrane normal is aligned with z before FT (disable with -raw)

    std::string base_dir = "";
    std::string ftmethod = "DFT"; // FTMethod of the FourierTransform currently being run
    vector<std::string> ftmethod_list = {"DFT"}; // sweep of FTMethods to run (populated by -ftmethod)
    vector<std::string> analysis;

};