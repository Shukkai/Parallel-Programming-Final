// main.cpp
#include "header/ga_cpu.h"
#include <iostream>
#include <utility>
#include <chrono>

using namespace std;

int main(int argc, char *argv[])
{

    GeneticTSP ga(100, 200000, 0.5, 0.8);

    // read file
    std::string filename = (argc > 1) ? argv[1] : "a280.tsp";
    if (!ga.readfile(filename))
    {
        std::cerr << "Failed to read TSP file." << std::endl;
        return 1;
    }

    // Solve TSP
    std::cout << "\nSolving TSP..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<int> bestTour = ga.solve();

    // Print result
    std::cout << "\nBest tour found:\n";
    for (int city : bestTour)
    {
        std::cout << city << " -> ";
    }

    std::cout << bestTour[0] << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "\nSolution found in " << duration.count() << " milliseconds" << std::endl;
    ;

    return 0;
}