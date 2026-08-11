#!/bin/bash
set -euo pipefail
rm -f ANA
./compile.sh

base_dir="/lustre/astro/semygind/Data/MembraneInclusion/"


FrameTension=IsotropicFrameTension
: '
#/lustre/astro/semygind/Data/MembraneInclusion/IsotropicFrameTension/FromFlat/VFMCoupling/VFFno/EPS1/
InitialConfiguration=FromFlat
Simulation=VFMCoupling
for VFF in no
do
    for EPS in 1 2 5
    do
        for KP in 1 10
        do
            #VFF="$EPS"
            KL=$(awk "BEGIN {printf \"%g\", $KP/2}")

            path_data="${base_dir}${FrameTension}/${InitialConfiguration}/${Simulation}/VFF${VFF}/EPS${EPS}/KP${KP}KL${KL}"
            path_destination="${FrameTension}/${InitialConfiguration}/${Simulation}/VFF${VFF}/EPS${EPS}/KP${KP}KL${KL}"

            mkdir -p $path_destination
            cp run.sh run_temporary.sh
            sed -i "s/simulation/KP${KP}KL${KL}/g" run_temporary.sh
            sed -i "s|pathdata|${path_data}/|g" run_temporary.sh
            cp run_temporary.sh "$path_destination/run.sh"
            
            rm -f "$path_destination/ANA"
            cp ANA $path_destination
            (
                cd $path_destination
                sbatch run.sh
            )
        done
    done
done
'
: '
InitialConfiguration=FromFlat
Simulation=OnlyKP
for KP in 0 0.1 0.5 1 5 10
do
    path_data="${base_dir}${FrameTension}/${InitialConfiguration}/${Simulation}/KP${KP}"
    path_destination="${FrameTension}/${InitialConfiguration}/${Simulation}/KP${KP}"

    mkdir -p $path_destination
    cp run.sh run_temporary.sh
    sed -i "s/simulation/KP${KP}/g" run_temporary.sh
    sed -i "s|pathdata|${path_data}/|g" run_temporary.sh
    cp run_temporary.sh "$path_destination/run.sh"
    
    rm -f "$path_destination/ANA"
    cp ANA $path_destination
    (
        cd $path_destination
        sbatch run.sh
    )
done
'
: '
InitialConfiguration=FromFlat
Simulation=VFMCoupling
kappa=20
for VFF in no
do
    for EPS in 1
    do
        for KP in 1 2 5 10 20 100
        do
            #VFF="$EPS"
            KL=$(awk "BEGIN {printf \"%g\", $KP/2}")

            path_data="${base_dir}${FrameTension}/${InitialConfiguration}/${Simulation}/kappa${kappa}/VFF${VFF}/EPS${EPS}/KP${KP}KL${KL}"
            path_destination="${FrameTension}/${InitialConfiguration}/${Simulation}/kappa${kappa}/VFF${VFF}/EPS${EPS}/KP${KP}KL${KL}"

            mkdir -p $path_destination
            cp run.sh run_temporary.sh
            sed -i "s/simulation/KP${KP}KL${KL}/g" run_temporary.sh
            sed -i "s|pathdata|${path_data}/|g" run_temporary.sh
            cp run_temporary.sh "$path_destination/run.sh"
            
            rm -f "$path_destination/ANA"
            cp ANA $path_destination
            (
                cd $path_destination
                sbatch run.sh
            )
        done
    done
done
'

: '
#/lustre/astro/semygind/Data/MembraneInclusion/IsotropicFrameTension/FromEquilibrium/Activity/VFFno/EPS1/KP10KL5/XI-1/
InitialConfiguration=FromEquilibrium
Simulation=Activity
for VFF in no
do
    for EPS in 1
    do
        for KP in 1 10
        do
            for XI in -1 1
            do
                KL=$(awk "BEGIN {printf \"%g\", $KP/2}")

                path_data="${base_dir}${FrameTension}/${InitialConfiguration}/${Simulation}/VFF${VFF}/EPS${EPS}/KP${KP}KL${KL}/XI${XI}"
                path_destination="${FrameTension}/${InitialConfiguration}/${Simulation}/VFF${VFF}/EPS${EPS}/KP${KP}KL${KL}/XI${XI}"

                mkdir -p $path_destination
                cp run.sh run_temporary.sh
                sed -i "s/simulation/XI${XI}KP${KP}KL${KL}/g" run_temporary.sh
                sed -i "s|pathdata|${path_data}/|g" run_temporary.sh
                cp run_temporary.sh "$path_destination/run.sh"
                
                rm -f "$path_destination/ANA"
                cp ANA $path_destination
                (
                    cd $path_destination
                    sbatch run.sh
                )
            done
        done
    done
done
'
: '
kappa=20
InitialConfiguration=FromFlat
Simulation=OnlyVectorField
for VFF in no
do
    for EPS in 0 0.1 0.5 0.8 1 2 5 10
    do
        #VFF="$EPS"
        path_data="${base_dir}${FrameTension}/${InitialConfiguration}/${Simulation}/kappa${kappa}/VFF${VFF}/EPS${EPS}"
        path_destination="${FrameTension}/${InitialConfiguration}/${Simulation}/kappa${kappa}/VFF${VFF}/EPS${EPS}"

        mkdir -p $path_destination
        cp run.sh run_temporary.sh
        sed -i "s/simulation/EPS${EPS}VFF${VFF}/g" run_temporary.sh
        sed -i "s|pathdata|${path_data}/|g" run_temporary.sh
        cp run_temporary.sh "$path_destination/run.sh"
        
        rm -f "$path_destination/ANA"
        cp ANA $path_destination
        (
            cd $path_destination
            sbatch run.sh
        )
    done
done
'
kappa=10
InitialConfiguration=FromFlat
Simulation=OnlyVectorField
for VFF in no
do
    for EPS in 0 0.5 0.8 1 2 5 10
    do
        #VFF="$EPS"
        path_data="${base_dir}${FrameTension}/${InitialConfiguration}/${Simulation}/kappa${kappa}/VFF${VFF}/EPS${EPS}"
        path_destination="${FrameTension}/${InitialConfiguration}/${Simulation}/kappa${kappa}/VFF${VFF}/EPS${EPS}"

        mkdir -p $path_destination
        cp run.sh run_temporary.sh
        sed -i "s/simulation/EPS${EPS}VFF${VFF}/g" run_temporary.sh
        sed -i "s|pathdata|${path_data}/|g" run_temporary.sh
        cp run_temporary.sh "$path_destination/run.sh"
        
        rm -f "$path_destination/ANA"
        cp ANA $path_destination
        (
            cd $path_destination
            sbatch run.sh
        )
    done
done

: '
InitialConfiguration=FromFlat
Simulation=OnlyVectorField
for VFF in 1 2 5
do
    for EPS in 1
    do
        #VFF="$EPS"
        path_data="${base_dir}${FrameTension}/${InitialConfiguration}/${Simulation}/VFF${VFF}/EPS${EPS}"
        path_destination="${FrameTension}/${InitialConfiguration}/${Simulation}/VFF${VFF}/EPS${EPS}"

        mkdir -p $path_destination
        cp run.sh run_temporary.sh
        sed -i "s/simulation/EPS${EPS}VFF${VFF}/g" run_temporary.sh
        sed -i "s|pathdata|${path_data}/|g" run_temporary.sh
        cp run_temporary.sh "$path_destination/run.sh"
        
        rm -f "$path_destination/ANA"
        cp ANA $path_destination
        (
            cd $path_destination
            sbatch run.sh
        )
    done
done
'
: '
kappa=20
InitialConfiguration=FromFlat
Simulation=OnlyMembrane
path_data="${base_dir}${FrameTension}/${InitialConfiguration}/${Simulation}/kappa${kappa}"
path_destination="${FrameTension}/${InitialConfiguration}/${Simulation}/kappa${kappa}"

mkdir -p $path_destination
cp run.sh run_temporary.sh
sed -i "s/simulation/OnlyMembrane/g" run_temporary.sh
sed -i "s|pathdata|${path_data}/|g" run_temporary.sh
cp run_temporary.sh "$path_destination/run.sh"

rm -f "$path_destination/ANA"
cp ANA $path_destination
(
    cd $path_destination
    sbatch run.sh
)
'
kappa=10
InitialConfiguration=FromFlat
Simulation=OnlyMembrane
path_data="${base_dir}${FrameTension}/${InitialConfiguration}/${Simulation}/kappa${kappa}"
path_destination="${FrameTension}/${InitialConfiguration}/${Simulation}/kappa${kappa}"

mkdir -p $path_destination
cp run.sh run_temporary.sh
sed -i "s/simulation/OnlyMembrane/g" run_temporary.sh
sed -i "s|pathdata|${path_data}/|g" run_temporary.sh
cp run_temporary.sh "$path_destination/run.sh"

rm -f "$path_destination/ANA"
cp ANA $path_destination
(
    cd $path_destination
    sbatch run.sh
)

: '
FrameTension=AnisotropicFrameTension

InitialConfiguration=FromFlat
Simulation=OnlyVectorField
for EPS in 0 0.1 1 2 5 10
do
    VFF="$EPS"
    path_data="${base_dir}${FrameTension}/${InitialConfiguration}/${Simulation}/VFF${VFF}/EPS${EPS}"
    path_destination="${FrameTension}/${InitialConfiguration}/${Simulation}/VFF${VFF}/EPS${EPS}"

    mkdir -p $path_destination
    cp run.sh run_temporary.sh
    sed -i "s/simulation/EPS${EPS}VFF${VFF}/g" run_temporary.sh
    sed -i "s|pathdata|${path_data}/|g" run_temporary.sh
    cp run_temporary.sh "$path_destination/run.sh"
    
    rm -f "$path_destination/ANA"
    cp ANA $path_destination
    (
        cd $path_destination
        sbatch run.sh
    )
done

InitialConfiguration=FromEquilibrium
Simulation=OnlyMembrane
for EPS in no
do
    VFF="$EPS"
    path_data="${base_dir}${FrameTension}/${InitialConfiguration}/${Simulation}/VFF${VFF}/EPS${EPS}"
    path_destination="${FrameTension}/${InitialConfiguration}/${Simulation}/VFF${VFF}/EPS${EPS}"

    mkdir -p $path_destination
    cp run.sh run_temporary.sh
    sed -i "s/simulation/EPS${EPS}VFF${VFF}/g" run_temporary.sh
    sed -i "s|pathdata|${path_data}/|g" run_temporary.sh
    cp run_temporary.sh "$path_destination/run.sh"
    
    rm -f "$path_destination/ANA"
    cp ANA $path_destination
    (
        cd $path_destination
        sbatch run.sh
    )
done
'