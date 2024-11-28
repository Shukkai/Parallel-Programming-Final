# Compiler settings
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11 -O3 -Iinclude

# Directories
LIB_DIR = lib
INCLUDE_DIR = include
BUILD_DIR = build

# Project files
TARGET = tsp
OBJECTS = $(patsubst $(LIB_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(wildcard $(LIB_DIR)/*.cpp)) $(BUILD_DIR)/main.o

# Default input file
INPUT_FILE ?= tsp_graph/ch130.tsp
INPUT_TYPE ?= ga

# Main target
all: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

$(BUILD_DIR)/main.o: main.cpp $(INCLUDE_DIR)/*.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(LIB_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Run the program with input file
run: all
	./$(TARGET) $(INPUT_FILE) $(INPUT_TYPE)

# Clean built files
clean:
	rm -rf $(BUILD_DIR) $(TARGET)


# Help target
help:
	@echo "Available targets:"
	@echo "  make        - Build the TSP solver"
	@echo "  make run    - Run the solver with default input (a280.tsp)"
	@echo "  make clean  - Remove built files"
	@echo ""
	@echo "To specify a different input file:"
	@echo "  make run INPUT_FILE=your_file.tsp"

# # Specify dependencies
# main.o: main.cpp header/reader.h header/solver.h
# reader.o: reader.cpp header/reader.h
# solver.o: solver.cpp header/solver.h header/reader.h

.PHONY: all clean run help