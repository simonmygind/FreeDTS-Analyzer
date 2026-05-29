#!/bin/bash

Upper_Limit_HOURS=3600

while read -r jobid start_time; do
    start_epoch=$(date -d "$start_time" +%s)
    now=$(date +%s)
    elapsed=$((now - start_epoch))
    
    if [ $elapsed -lt $Upper_Limit_HOURS ]; then
        echo "Cancelling job $jobid (elapsed: $((elapsed/60)) minutes)"
        scancel $jobid
    fi
done < <(squeue -u $USER -h -o "%i %S")