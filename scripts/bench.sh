#!/bin/bash

total_time=0
iterations=500
threads=32
successful_runs=0

for ((i = 0 ; i < iterations ; i++)); do
    # Capture the output
    output=$(make run-omp INPUT=data/1000-res/3JB9.pdb OUTPUT=results/ THREADS=$threads 2>&1)
    
    # Extract the time - looking for pattern like "64|0.00034(s)"
    # This assumes the number before | is the thread count
    exec_time=$(echo "$output" | grep -oP "${threads}\|\K[0-9]+\.?[0-9]*(?=\(s\))")
    
    # Verify we got a numeric value
    if [[ -z "$exec_time" ]] || [[ ! "$exec_time" =~ ^[0-9]+\.?[0-9]*$ ]]; then
        echo "Warning: Invalid time value extracted on iteration $((i+1)): '$exec_time'"
        echo "Output was: $output"
        continue
    fi
    
    # Add to the running total
    total_time=$(echo "$total_time + $exec_time" | bc -l)
    successful_runs=$((successful_runs + 1))
    echo "$exec_time" 
    # Optional: show progress
    if (( (i + 1) % 100 == 0 )); then
        echo "Completed $((i+1))/$iterations iterations..."
    fi
done

# Calculate the average
if [ $successful_runs -gt 0 ]; then
    average_time=$(echo "scale=6; $total_time / $successful_runs" | bc -l)
    echo "Average execution time over $successful_runs successful runs: $average_time s"
else
    echo "Error: No successful runs to average"
fi
