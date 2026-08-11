version=version_1
destination=/lustre/astro/semygind/Analyzer/version_1/

cd /lustre/astro/semygind/Analyzer/${version}/dts_ana/Build
cmake -DCMAKE_BUILD_TYPE=Release .. 
make
cp ANA $destination