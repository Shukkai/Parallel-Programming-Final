# Compiler settings
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11 -O2 

# Project files
TARGET = tsp
SOURCES = main.cpp reader.cpp solver.cpp
HEADERS = reader.h solver.h
OBJECTS = $(SOURCES:.cpp=.o)

# Default input file
INPUT_FILE ?= tsp_graph/a280.tsp

# Main target
all: $(TARGET)

# Compile source files to object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link object files to create executable
$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET)

# Run the program with input file
run: $(TARGET)
	./$(TARGET) $(INPUT_FILE)

# Clean built files
clean:
	rm -f $(TARGET) $(OBJECTS)

# Help target
help:
	@echo "Available targets:"
	@echo "  make        - Build the TSP solver"
	@echo "  make run    - Run the solver with default input (a280.tsp)"
	@echo "  make clean  - Remove built files"
	@echo ""
	@echo "To specify a different input file:"
	@echo "  make run INPUT_FILE=your_file.tsp"

# Specify dependencies
main.o: main.cpp header/reader.h header/solver.h
reader.o: reader.cpp header/reader.h
solver.o: solver.cpp header/solver.h header/reader.h

.PHONY: all clean run help