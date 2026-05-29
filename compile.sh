
version=version_1
destination=/lustre/astro/semygind/Analyzer/version_1/
cd /lustre/astro/semygind/Analyzer/${version}/dts_ana
g++ -std=c++17 -O3 -o ANA *.cpp
mv ANA $destination