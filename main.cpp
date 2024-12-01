#include "ACO.h"
#include "ACO_omp.h"
#include "ACO_thread.h"
#include "ga_cpu.h"
#include "reader.h"
#include "solver.h"
#include <chrono>
#include "ga_omp.h"

enum parallel_type
{
    SERIAL,
    OMP,
    THREAD,
    CUDA
};

inline int hashit(std::string const &inString)
{
    if (inString == "serial")
        return SERIAL;
    if (inString == "omp")
        return OMP;
    if (inString == "thread")
        return THREAD;
    if (inString == "cuda")
        return CUDA;
    return SERIAL;
}

int main(int argc, char *argv[])
{
    std::string filename = (argc > 0) ? argv[1] : "a280.tsp";
    std::string type = (argc > 1) ? argv[2] : "aco";
    std::string parallel = (argc > 2) ? argv[3] : "serial";
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

        if (parallel == "omp")
        {
            ACOOmp solver(reader.getPoints());
            solver.solve();
            bestTour = solver.getTour();
            bestDistance = solver.getDistance();
        }
        else if (parallel == "thread")
        {
            ACOThread solver(reader.getPoints(), 4);
            solver.solve();
            bestTour = solver.getTour();
            bestDistance = solver.getDistance();
        }
        else
        {
            ACO solver(reader.getPoints());
            solver.solve();
            bestTour = solver.getTour();
            bestDistance = solver.getDistance();
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
        if (parallel == "omp")
        {
            ga_omp gaomp(reader.getPoints(), 100, 10000, 0.05, 0.8);
            start = std::chrono::high_resolution_clock::now();
            result = gaomp.solve();
        }
        else
        {
            GeneticTSP ga(reader.getPoints(), 100, 10000, 0.05, 0.8);
            start = std::chrono::high_resolution_clock::now();
            result = ga.solve();
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        bestTour = result.first;
        bestDistance = result.second;
        elapsedTime = duration.count();
    }

    // write results to file
    std::string res_fname = "res/" + type;
    switch (hashit(parallel))
    {
    case 1:
        res_fname += "_omp.txt";
        break;
    case 2:
        res_fname += "_thread.txt";
        break;
    case 3:
        res_fname += "_cuda.txt";
        break;
    default:
        res_fname += ".txt";
        break;
    }

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