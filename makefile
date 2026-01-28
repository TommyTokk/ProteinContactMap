# Compiler settings
CXX = g++
MPICXX = /usr/local/openmpi/bin/mpic++
MPIRUN = /usr/local/openmpi/bin/mpirun
CXXFLAGS = -std=c++20 -Wall -O3
OMPFLAGS = -fopenmp
INCLUDES = -I.

# Directories
SRC_DIR = src
LIB_DIR = lib
BIN_DIR = bin
BUILD_DIR = build

# Source files
UTILS_SRC = $(SRC_DIR)/utils.cpp

# Object files
UTILS_OBJ = $(BUILD_DIR)/utils.o
UTILS_OMP_OBJ = $(BUILD_DIR)/utils_omp.o
UTILS_MPI_OBJ = $(BUILD_DIR)/utils_mpi.o

# Executables
SEQ_EXEC = $(BIN_DIR)/seq_main
OMP_EXEC = $(BIN_DIR)/omp_main
MPI_EXEC = $(BIN_DIR)/mpi_main

# Default target
.PHONY: all
all: $(SEQ_EXEC) $(OMP_EXEC) $(MPI_EXEC)

# Create necessary directories
$(BIN_DIR) $(BUILD_DIR):
	mkdir -p $@

# Sequential version
$(UTILS_OBJ): $(UTILS_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(SEQ_EXEC): seq_main.cpp $(UTILS_OBJ) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

# OpenMP version
$(UTILS_OMP_OBJ): $(UTILS_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(INCLUDES) -c $< -o $@

$(OMP_EXEC): omp_main.cpp $(UTILS_OMP_OBJ) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(INCLUDES) $^ -o $@

# MPI version
$(UTILS_MPI_OBJ): $(UTILS_SRC) | $(BUILD_DIR)
	$(MPICXX) $(CXXFLAGS) $(OMPFLAGS) $(INCLUDES) -c $< -o $@

$(MPI_EXEC): mpi_main.cpp $(UTILS_MPI_OBJ) | $(BIN_DIR)
	$(MPICXX) $(CXXFLAGS) $(OMPFLAGS) $(INCLUDES) $^ -o $@

# Individual build targets
.PHONY: seq omp mpi
seq: $(SEQ_EXEC)
omp: $(OMP_EXEC)
mpi: $(MPI_EXEC)

# Run targets
.PHONY: run-seq run-omp run-mpi
run-seq: $(SEQ_EXEC)
	@if [ -z "$(INPUT)" ] || [ -z "$(OUTPUT)" ]; then \
		echo "Error: Please specify INPUT and OUTPUT variables"; \
		exit 1; \
	fi
	$(SEQ_EXEC) $(INPUT) $(OUTPUT)

run-omp: $(OMP_EXEC)
	@if [ -z "$(INPUT)" ] || [ -z "$(OUTPUT)" ] || [ -z "$(THREADS)" ]; then \
		echo "Error: Please specify INPUT, OUTPUT, and THREADS variables"; \
		exit 1; \
	fi
	$(OMP_EXEC) $(INPUT) $(OUTPUT) $(THREADS)

run-mpi: $(MPI_EXEC)
	@if [ -z "$(INPUT)" ] || [ -z "$(OUTPUT)" ] || [ -z "$(THREADS)" ] || [ -z "$(PROCS)" ]; then \
		echo "Error: Please specify INPUT, OUTPUT, THREADS, and PROCS variables"; \
		exit 1; \
	fi
	$(MPIRUN) -np $(PROCS) $(MPI_EXEC) $(INPUT) $(OUTPUT) $(THREADS)

# Clean targets
.PHONY: clean clean-build clean-bin clean-all
clean-build:
	rm -rf $(BUILD_DIR)

clean-bin:
	rm -rf $(BIN_DIR)

clean: clean-build clean-bin

clean-all: clean
	rm -f seq_main omp_main mpi_main