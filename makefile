# Compiler settings
CXX = g++
MPICXX = mpicxx
CXXFLAGS = -Wall -Wextra -std=c++17 -O3 -Iinclude -fopenmp -pthread
MPICXXFLAGS = -Wall -Wextra -std=c++17 -O3 -Iinclude -fPIC

# Directories
LIB_DIR = lib
INCLUDE_DIR = include
BUILD_DIR = build

# Project files
TARGET = tsp
OBJECTS = $(patsubst $(LIB_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(wildcard $(LIB_DIR)/*.cpp)) $(BUILD_DIR)/main.o
OBJECTS_MPI = $(patsubst $(LIB_DIR)/%.cc,$(BUILD_DIR)/%.o,$(wildcard $(LIB_DIR)/*.cc)) $(BUILD_DIR)/main.o
OBJECTS_MPI_C = $(BUILD_DIR)/ACO.o $(BUILD_DIR)/ga_cpu.o $(BUILD_DIR)/reader.o $(BUILD_DIR)/solver.o

# Default input file
testcase ?= tsp_graph/a280.tsp
algo ?= ga	# ga, aco
parallel ?= serial	# serial, omp, thread, mpi

ifeq ($(parallel), mpi)
all: clean $(OBJECTS_MPI) $(OBJECTS_MPI_C)
	$(MPICXX) $(MPICXXFLAGS) $(OBJECTS_MPI) $(OBJECTS_MPI_C) -o $(TARGET)
else
all: clean $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)
endif

$(BUILD_DIR)/main.o: main.cpp $(INCLUDE_DIR)/*.h
ifneq ($(parallel), mpi)
	$(CXX) $(CXXFLAGS) -c $< -o $@
endif

$(BUILD_DIR)/main.o: main.cc $(INCLUDE_DIR)/*.hpp
ifeq ($(parallel), mpi)
	$(MPICXX) $(MPICXXFLAGS) -c $< -o $@
endif

$(BUILD_DIR)/%.o: $(LIB_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
ifneq ($(parallel), mpi)
	$(CXX) $(CXXFLAGS) -c $< -o $@
else
	$(MPICXX) $(MPICXXFLAGS) -c $< -o $@
endif

$(BUILD_DIR)/%.o: $(LIB_DIR)/%.cc
	@mkdir -p $(BUILD_DIR)
	$(MPICXX) $(MPICXXFLAGS) -c $< -o $@

# Run the program with input file
run: clean all
	./$(TARGET) $(testcase) $(algo) $(parallel)

# Clean built files
clean:
	rm -rf $(BUILD_DIR) $(TARGET) *.out

cuda:
	nvcc -std=c++17 -O3 -Xcompiler '-fPIC -use_fast_math' aco.cu -o aco_cuda.out
	./aco_cuda.out $(testcase)
