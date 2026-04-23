#!/bin/bash

# Benchmark script for ProteinContactMapC
# Executes seq_main, omp_main, and mpi_main with all models in data folder
# Creates CSV files with benchmark results

# Configuration
ITERATIONS=50 
THREADS_LIST=(4 8 16 32 64 128 255)
PROCESSORS_LIST=(4 8 16 32 64)
DATA_FOLDERS=("data/10-res" "data/100-res" "data/1000-res")
OUTPUT_DIR="results"
RESULTS_DIR="benchmark_results"

# Parse command line arguments
RUN_SEQ=false
RUN_OMP=false
RUN_MPI=false

if [ $# -eq 0 ]; then
    echo "Usage: $0 [--all|--seq|--omp|--mpi]"
    echo "  --all  : Run all benchmarks (sequential, OpenMP, and MPI)"
    echo "  --seq  : Run only sequential benchmark"
    echo "  --omp  : Run only OpenMP benchmark"
    echo "  --mpi  : Run only MPI benchmark"
    exit 1
fi

for arg in "$@"; do
    case $arg in
        --all)
            RUN_SEQ=true
            RUN_OMP=true
            RUN_MPI=true
            ;;
        --seq)
            RUN_SEQ=true
            ;;
        --omp)
            RUN_OMP=true
            ;;
        --mpi)
            RUN_MPI=true
            ;;
        *)
            echo "Unknown option: $arg"
            echo "Usage: $0 [--all|--seq|--omp|--mpi]"
            exit 1
            ;;
    esac
done

# Create results directory
mkdir -p "$RESULTS_DIR"

# Create log file with timestamp
LOG_FILE="$RESULTS_DIR/benchmark_$(date +%Y%m%d_%H%M%S).log"
touch "$LOG_FILE"

# Function to log messages
log() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $*" >> "$LOG_FILE"
}

# Function to get model size from PDB file
get_model_size() {
    local file=$1
    # Count CA (alpha carbon) atoms in the PDB file
    local size=$(grep "^ATOM.*CA" "$file" | wc -l)
    echo "$size"
}

# Function to extract filename without extension
get_filename() {
    basename "$1" .pdb
}

# Function to parse alg-time from outputs
parse_alg_time() {
    local output=$1
    echo "$output" | grep "alg-time" | head -n 1 | awk -F'|' '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}'
}

log "==================================================================="
log "Starting Protein Contact Map Benchmarks"
log "Log file: $LOG_FILE"
log "Run Sequential: $RUN_SEQ"
log "Run OpenMP: $RUN_OMP"
log "Run MPI: $RUN_MPI"
log "==================================================================="
log ""

# ===================================================================
# SEQUENTIAL VERSION BENCHMARK
# ===================================================================
if [ "$RUN_SEQ" = true ]; then
log "==================================================================="
log "Benchmarking Sequential Version (seq_main)"
log "=================================================================="

for data_folder in "${DATA_FOLDERS[@]}"; do
    if [ ! -d "$data_folder" ]; then
        log "Warning: Directory $data_folder does not exist, skipping..."
        continue
    fi
    
    for pdb_file in "$data_folder"/*.pdb; do
        if [ ! -f "$pdb_file" ]; then
            log "Warning: No PDB files found in $data_folder"
            continue
        fi
        
        model_name=$(get_filename "$pdb_file")
        model_size=$(get_model_size "$pdb_file")
        
        # Create model-specific directory
        MODEL_DIR="$RESULTS_DIR/$model_name"
        mkdir -p "$MODEL_DIR"
        
        SEQ_CSV="$MODEL_DIR/sequential_${model_name}_${model_size}res.csv"
        echo "iteration,alg_time,time_unit" > "$SEQ_CSV"
        
        log "Processing: $model_name (size: $model_size residues)"
        
        for ((i=1; i<=ITERATIONS; i++)); do
            output=$(make run-seq INPUT="$pdb_file" OUTPUT="$OUTPUT_DIR" 2>&1)
            
            if [ $? -ne 0 ]; then
                log "  Error on iteration $i for $model_name - skipping to next model"
                continue 2
            fi
            
            alg_time=$(parse_alg_time "$output")
            
            if [ -z "$alg_time" ]; then
                log "  Error: Could not extract time for iteration $i - skipping to next model"
                continue 2
            fi
            
            echo "$i,$alg_time,s" >> "$SEQ_CSV"
            
            # Log progress at intervals
            if [ $((i % 50)) -eq 0 ] || [ $i -eq $ITERATIONS ]; then
                log "  Progress: $i/$ITERATIONS iterations completed"
            fi
        done
        log "  Results saved to $SEQ_CSV"
        log ""
    done
done

log "Sequential benchmark complete."
log ""
fi

# ===================================================================
# OPENMP VERSION BENCHMARK
# ===================================================================
if [ "$RUN_OMP" = true ]; then
log "==================================================================="
log "Benchmarking OpenMP Version (omp_main)"
log "=================================================================="

for data_folder in "${DATA_FOLDERS[@]}"; do
    if [ ! -d "$data_folder" ]; then
        continue
    fi
    
    for pdb_file in "$data_folder"/*.pdb; do
        if [ ! -f "$pdb_file" ]; then
            continue
        fi
        
        model_name=$(get_filename "$pdb_file")
        model_size=$(get_model_size "$pdb_file")
        
        # Create model-specific directory
        MODEL_DIR="$RESULTS_DIR/$model_name"
        mkdir -p "$MODEL_DIR"
        
        OMP_CSV="$MODEL_DIR/openmp_${model_name}_${model_size}res.csv"
        echo "iteration,alg_time,time_unit,num_threads" > "$OMP_CSV"
        
        log "Processing: $model_name (size: $model_size residues)"
        
        total_iterations=$((${#THREADS_LIST[@]} * ITERATIONS))
        current_iteration=0
        
        for threads in "${THREADS_LIST[@]}"; do
            
            for ((i=1; i<=ITERATIONS; i++)); do
                output=$(make run-omp INPUT="$pdb_file" OUTPUT="$OUTPUT_DIR" THREADS="$threads" 2>&1)
                
                if [ $? -ne 0 ]; then
                    log "    Error on iteration $i with $threads threads - skipping to next model"
                    continue 3
                fi
                
                alg_time=$(parse_alg_time "$output")
                
                if [ -z "$alg_time" ]; then
                    log "    Error: Could not extract time for iteration $i - skipping to next model"
                    continue 3
                fi
                
                echo "$i,$alg_time,s,$threads" >> "$OMP_CSV"
                
                # Log progress at intervals
                current_iteration=$((current_iteration + 1))
                if [ $((current_iteration % 50)) -eq 0 ] || [ $current_iteration -eq $total_iterations ]; then
                    log "  Progress: $current_iteration/$total_iterations iterations completed (threads: $threads)"
                fi
            done
        done
        log "  Results saved to $OMP_CSV"
        log ""
    done
done

log "OpenMP benchmark complete."
log ""
fi

# ===================================================================
# MPI VERSION BENCHMARK
# ===================================================================
if [ "$RUN_MPI" = true ]; then
log "==================================================================="
log "Benchmarking MPI Version (mpi_main)"
log "=================================================================="

for data_folder in "${DATA_FOLDERS[@]}"; do
    if [ ! -d "$data_folder" ]; then
        continue
    fi
    
    for pdb_file in "$data_folder"/*.pdb; do
        if [ ! -f "$pdb_file" ]; then
            continue
        fi
        
        model_name=$(get_filename "$pdb_file")
        model_size=$(get_model_size "$pdb_file")
        
        # Create model-specific directory
        MODEL_DIR="$RESULTS_DIR/$model_name"
        mkdir -p "$MODEL_DIR"
        
        MPI_CSV="$MODEL_DIR/mpi_${model_name}_${model_size}res.csv"
        echo "iteration,alg_time,time_unit,num_processors,num_threads" > "$MPI_CSV"
        
        log "Processing: $model_name (size: $model_size residues)"
        
        total_iterations=$((${#PROCESSORS_LIST[@]} * ${#THREADS_LIST[@]} * ITERATIONS))
        current_iteration=0
        
        for procs in "${PROCESSORS_LIST[@]}"; do
            for threads in "${THREADS_LIST[@]}"; do
                
                for ((i=1; i<=ITERATIONS; i++)); do
                    output=$(make run-mpi INPUT="$pdb_file" OUTPUT="$OUTPUT_DIR" THREADS="$threads" PROCS="$procs" 2>&1 | grep -v "UCX  WARN" | grep -v "Rank:")
                    
                    if [ $? -ne 0 ]; then
                        log "    Error on iteration $i with $procs procs, $threads threads - skipping to next model"
                        continue 4
                    fi
                    
                    alg_time=$(parse_alg_time "$output")
                    
                    if [ -z "$alg_time" ]; then
                        log "    Error: Could not extract time for iteration $i - skipping to next model"
                        continue 4
                    fi
                    
                    echo "$i,$alg_time,s,$procs,$threads" >> "$MPI_CSV"
                    
                    # Log progress at intervals
                    current_iteration=$((current_iteration + 1))
                    if [ $((current_iteration % 50)) -eq 0 ] || [ $current_iteration -eq $total_iterations ]; then
                        log "  Progress: $current_iteration/$total_iterations iterations completed (procs: $procs, threads: $threads)"
                    fi
                done
            done
        done
        log "  Results saved to $MPI_CSV"
        log ""
    done
done

log "MPI benchmark complete."
log ""
fi

# ===================================================================
# SUMMARY
# ===================================================================
log "==================================================================="
log "Benchmark Complete!"
log "==================================================================="
log ""
log "Results have been saved to: $RESULTS_DIR/<model_name>/"
log ""
log "Directory structure:"
log "  $RESULTS_DIR/"
log "    ├── <model_name>/"
log "    │   ├── sequential_<model>_<size>res.csv"
log "    │   ├── openmp_<model>_<size>res.csv"
log "    │   └── mpi_<model>_<size>res.csv"
log ""
log "You can analyze these CSV files using the plotting scripts."
log "==================================================================="
