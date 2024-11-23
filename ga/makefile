# Compiler settings
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11 -O2 

# Project files
TARGET = ga
SOURCES = main.cpp ga_cpu.cpp
HEADERS = ga_cpu.h
OBJECTS = $(SOURCES:.cpp=.o)

# Default input file
INPUT_FILE ?= ../tsp_graph/a280.tsp

# Main target
all: $(TARGET)

# Compile source files to object files
$(OBJ_DIR)/%.o: %.cpp $(HEADERS)
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

# Specify dependencies
main.o: main.cpp header/ga_cpu.h
ga_cpu.o: ga_cpu.cpp header/ga_cpu.h

.PHONY: all clean run help