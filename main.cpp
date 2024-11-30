#include "ACO.h"
#include "reader.h"
#include "solver.h"
#include "ga_cpu.h"
#include <chrono>
#include <iostream>

int main(int argc, char *argv[])
{
    std::string filename = (argc > 0) ? argv[1] : "a280.tsp";
    std::string type = (argc > 1) ? argv[2] : "aco";
    if(type == "aco"){
        TSPReader reader;
        if (!reader.readFile(filename)) {
            std::cerr << "Failed to read TSP file." << std::endl;
            return 1;
        }
        // Print the original data
        // reader.printData();
        ACO solver(reader.getPoints());
        solver.solve();
        std::cout << solver.getDistance() << std::endl;
    }
    else{
        GeneticTSP ga(100, 10000, 0.05, 0.8);
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
    }

    // // Solve TSP
    // std::cout << "\nSolving TSP..." << std::endl;
    // auto start = std::chrono::high_resolution_clock::now();

    // TSPAlgorithm solver(reader.getPoints());

    // // First solve using Nearest Neighbor
    // solver.solveNearestNeighbor();
    // std::cout << "\nInitial solution (Nearest Neighbor):" << std::endl;
    // solver.printTour();

    // // Then improve using 2-opt
    // std::cout << "\nImproving solution with 2-opt..." << std::endl;
    // solver.improve2Opt();
    // std::cout << "\nFinal solution (after 2-opt):" << std::endl;
    // solver.printTour();

    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // std::cout << "\nSolution found in " << duration.count() << " milliseconds" << std::endl;

    return 0;
}