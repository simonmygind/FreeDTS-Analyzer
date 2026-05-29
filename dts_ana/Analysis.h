#pragma once
#include "TrajTSI.h"
#include "Config.h"
#include "FourierTransform.h"
#include <vector>
#include <iostream>
#include <array>
#include <memory>
#include <optional>

using std::cout, std::vector, std::array;


class Analysis
{
private:
    Config m_Config;
    std::unique_ptr<TrajTSI> m_TrajTSI;
    std::unique_ptr<FourierTransform> m_FourierTransform;

public:
    Analysis(const int argc,const char* argv[]) : m_Config(ParseArgs(argc,argv)){}
    
    Config ParseArgs(const int argc,const char* argv[]);
    void RunAnalysis();
};



