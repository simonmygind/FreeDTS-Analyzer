#pragma once
#include <vector>
#include <fstream>

using std::vector;
struct Config
{   
    int init_f = 1;
    int final_f = 1;
    int init_rep = 1; 
    int final_rep = 1;
    int N_k = 5;
    int N_t= 0;
    std::string base_dir = "";
    std::string ftmethod = "DFT";
    vector<std::string> analysis;

};