#!/bin/bash
KL=5

for KP in 0 0.1 0.5 1 5 10
do
    mkdir -p KP${KP}KL${KL}
    
    #Modifies run_new file
    cp run.sh run_temporary.sh
    sed -i "s/kp/$KP/g" run_temporary.sh
    sed -i "s/kl/$KL/g" run_temporary.sh
    cp run_temporary.sh "KP${KP}KL${KL}/run.sh"

    cp ANA KP${KP}KL${KL}
    cd KP${KP}KL${KL}
    sbatch run.sh
    cd ..
done