#include "Analysis.h"
#include <chrono>

using std::cout;
int main(const int argc,const char* argv[]){
    cout.setf(std::ios::unitbuf);
    auto start_totalruntime = std::chrono::high_resolution_clock::now();
    Analysis terminal_analysis(argc, argv);
    terminal_analysis.RunAnalysis();
    auto end_totalruntime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_totalruntime - start_totalruntime;
    cout << "Finished Analysis total runtime " << elapsed.count() << " s \n";
    return 0;
}