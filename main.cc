#include "reader.h"
#include "solver.h"
#include "ACO_mpi.hpp"
#include "ga_mpi.hpp"
#include <chrono>

int main(int argc, char *argv[])
{
    std::string filename = (argc > 1) ? argv[1] : "tsp_graph/a280.tsp";
    std::string type = (argc > 2) ? argv[2] : "aco";
    std::string parallel = (argc > 3) ? argv[3] : "mpi";
    TSPReader reader;
    if (!reader.readFile(filename))
    {
        std::cerr << "Failed to read TSP file." << std::endl;
        return 1;
    }

    std::vector<int> bestTour;
    double bestDistance = 0.0, elapsedTime = 0.0;
    std::cout << "Solving TSP...\n";
    if (type == "aco")
    {
        auto start = std::chrono::high_resolution_clock::now();

        if (parallel == "mpi")
        {
            int rank, size;
            MPI_Init(nullptr, nullptr);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            std::cout << rank << " " << size << std::endl;

            ACOMpi solver(reader.getPoints());
            solver.solve();
            bestTour = solver.getTour();
            bestDistance = solver.getDistance();

            MPI_Barrier(MPI_COMM_WORLD);
            if (rank != 0) {
                MPI_Finalize();
                return 0;
            } else {
                MPI_Finalize();
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        elapsedTime = duration.count();
    }
    else if (type == "ga")
    {
        // Solve TSP
        auto start = std::chrono::high_resolution_clock::now();
        std::pair<std::vector<int>, double> result;
        if (parallel == "mpi") {
            int rank, size;
            MPI_Init(nullptr, nullptr);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);

            ga_mpi gampi(reader.getPoints(), 100, 10000, 0.05, 0.8);
            start = std::chrono::high_resolution_clock::now();
            result = gampi.solve();

            MPI_Barrier(MPI_COMM_WORLD);
            if (rank != 0) {
                MPI_Finalize();
                return 0;
            } else {
                MPI_Finalize();
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        bestTour = result.first;
        bestDistance = result.second;
        elapsedTime = duration.count();
    }

    // write results to file
    std::string res_fname = "res/" + filename + "_" + type + "_mpi.txt";

    std::ofstream file;
    file.open(res_fname);
    file << "Testcase: " << filename << "\n";
    file << "Time taken: " << elapsedTime << " milliseconds\n";
    file << "Total distance: " << bestDistance << "\n";
    file << "Best tour found:\n";
    for (int city : bestTour)
    {
        file << city << "\n";
    }

    // print results to console
    std::cout << "Distance: " << bestDistance << "\nSolution found in " << elapsedTime << " milliseconds\n";
    
    return 0;
}
