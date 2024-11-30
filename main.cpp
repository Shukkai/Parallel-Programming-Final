#include "ACO.h"
#include "ga_cpu.h"
#include "reader.h"
#include "solver.h"
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char *argv[])
{
    std::string filename = (argc > 0) ? argv[1] : "a280.tsp";
    std::string type = (argc > 1) ? argv[2] : "aco";

    TSPReader reader;
    if (!reader.readFile(filename)) {
        std::cerr << "Failed to read TSP file." << std::endl;
        return 1;
    }

    std::cout << "\nSolving TSP..." << std::endl;
    if (type == "aco") {
        ACO solver(reader.getPoints());

        auto start = std::chrono::high_resolution_clock::now();
        solver.solve();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        // write results to file
        std::ofstream file;
        file.open("res/aco.txt");
        file << "Best tour found:\n";
        for (int city : solver.getTour()) {
            file << city << "\n";
        }

        std::cout << solver.getDistance() << "\nSolution found in " << duration.count() << " milliseconds" << std::endl;
    }
    else if (type == "ga") {
        GeneticTSP ga(100, 10000, 0.05, 0.8);

        // Solve TSP
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<int> bestTour = ga.solve();
        auto end = std::chrono::high_resolution_clock::now();

        // write results to file
        std::ofstream file;
        file.open("res/ga.txt");
        file << "Best tour found:\n";
        for (int city : bestTour) {
            file << city << "\n";
        }

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "\nSolution found in " << duration.count() << " milliseconds" << std::endl;
    }

    return 0;
}