#!/bin/bash

THRESHOLD=$((1 * 60 * 60))  # Measured in seconds
NAME=EPS
#1 10 0.1
while read -r jobid jobname start_time; do
    # Skip jobs not starting with "KP"
    [[ "$jobname" != ${NAME}* ]] && continue

    start_epoch=$(date -d "$start_time" +%s)
    now=$(date +%s)
    elapsed=$((now - start_epoch))

    if [ "$elapsed" -lt "$THRESHOLD" ]; then
        echo "Cancelling $jobid ($jobname, running $((elapsed / 3600))h $((elapsed % 3600 / 60))m)"
        scancel "$jobid"
    fi
done < <(squeue -u "$USER" -h -o "%i %j %S")