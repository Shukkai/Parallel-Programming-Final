// ga_mpi.cpp
#include "ga_mpi.h"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>
#include <random>
#include <string>

ga_mpi::ga_mpi(const std::vector<Point> &pts, int popSize, int gens, double mutRate, double crossRate)
    : cities(pts), populationSize(popSize), numGenerations(gens), mutationRate(mutRate), crossoverRate(crossRate)
{
    // MPI_Init(nullptr, nullptr);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);   

    // srand(rank + 1);
}

ga_mpi::~ga_mpi() {
    // MPI_Finalize();
}

double ga_mpi::distance(const Point &c1, const Point &c2)
{
    double dx = c1.x - c2.x;
    double dy = c1.y - c2.y;
    return sqrt(dx * dx + dy * dy);
}

double ga_mpi::calculateFitness(const std::vector<int> &path)
{
    double total = 0;
    for (size_t i = 0; i < path.size() - 1; i++)
    {
        total += distance(cities[path[i]], cities[path[i + 1]]);
    }
    total += distance(cities[path.back()], cities[path[0]]); // Return to start
    return total;
}

std::vector<int> ga_mpi::selectParent(const std::vector<std::vector<int>> &currentPopulation, std::mt19937 &rng)
{
    int tournamentSize = 5;
    std::uniform_int_distribution<> dist(0, currentPopulation.size() - 1);
    double bestFitness = std::numeric_limits<double>::infinity();
    std::vector<int> bestPath;

    for (int i = 0; i < tournamentSize; i++)
    {
        int idx = dist(rng);
        double fitness = calculateFitness(currentPopulation[idx]);
        if (fitness < bestFitness)
        {
            bestFitness = fitness;
            bestPath = currentPopulation[idx];
        }
    }
    return bestPath;
}

std::vector<int> ga_mpi::crossover(const std::vector<int> &parent1, const std::vector<int> &parent2, double crossoverRate, std::mt19937 &rng)
{
    std::uniform_real_distribution<> dis(0.0, 1.0);
    if (dis(rng) > crossoverRate)
    {
        return parent1;
    }

    int n = parent1.size();
    std::vector<int> child(n, -1);
    std::uniform_int_distribution<> dist(0, n - 1);
    int start = dist(rng);
    int end = dist(rng);
    if (start > end)
        std::swap(start, end);

    for (int i = start; i <= end; i++)
    {
        child[i] = parent1[i];
    }

    int j = 0;
    for (int i = 0; i < n; i++)
    {
        if (i >= start && i <= end)
            continue;
        while (std::find(child.begin() + start, child.begin() + end + 1, parent2[j]) != child.begin() + end + 1)
        {
            j++;
        }
        child[i] = parent2[j++];
    }
    return child;
}

void ga_mpi::mutate(std::vector<int> &path, double mutationRate, std::mt19937 &rng)
{
    std::uniform_real_distribution<> dis(0.0, 1.0);
    if (dis(rng) < mutationRate)
    {
        std::uniform_int_distribution<> dist(0, path.size() - 1);
        int i = dist(rng);
        int j = dist(rng);
        std::swap(path[i], path[j]);
    }
}

void ga_mpi::improve2Opt(std::vector<int> &tour)
{
    int n = tour.size();
    bool improved = false;
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            double beforeDistance =
                distance(cities[tour[i]], cities[tour[i + 1]]) + distance(cities[tour[j]], cities[tour[(j + 1) % n]]);

            double afterDistance =
                distance(cities[tour[i]], cities[tour[j]]) + distance(cities[tour[i + 1]], cities[tour[(j + 1) % n]]);

            if (afterDistance < beforeDistance)
            {
                std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                improved = true;
                break;
            }
        }
        if (improved)
            break;
    }
}

std::pair<std::vector<int>, double> ga_mpi::solve()
{
    int individualsPerProcess = populationSize / size;
    int remainder = populationSize % size;
    if (rank < remainder) {
        individualsPerProcess += 1;
    }

    std::mt19937 rng(rank + 1);

    std::vector<std::vector<int>> localPopulation;
    for (int i = 0; i < individualsPerProcess; i++) {
        std::vector<int> tour(cities.size());
        for (size_t j = 0; j < cities.size(); j++)
            tour[j] = j;
        std::shuffle(tour.begin(), tour.end(), rng);
        localPopulation.push_back(tour);
    }

    std::vector<int> localBestTour;
    double localBestDistance = std::numeric_limits<double>::infinity();
    std::vector<int> globalBestTour(cities.size());
    double globalBestDistance = std::numeric_limits<double>::infinity();

    // Main GA loop
    for (int gen = 0; gen < numGenerations; gen++)
    {
        // Evaluate fitness for the local population
        for (auto &tour : localPopulation)
        {
            double fitness = calculateFitness(tour);
            if (fitness < localBestDistance)
            {
                localBestDistance = fitness;
                localBestTour = tour;
            }
        }

        // Gather the best distances from all processes
        double allBestDistances[size];
        MPI_Allgather(&localBestDistance, 1, MPI_DOUBLE, allBestDistances, 1, MPI_DOUBLE, MPI_COMM_WORLD);

        // Gather the best tours from all processes
        int tourSize = cities.size();
        std::vector<int> allBestTours;
        if (rank == 0) {
            allBestTours.resize(size * tourSize);
        }
        MPI_Gather(localBestTour.data(), tourSize, MPI_INT,
                   allBestTours.data(), tourSize, MPI_INT, 0, MPI_COMM_WORLD);

        // Determine the global best tour
        if (rank == 0) {
            globalBestDistance = allBestDistances[0];
            globalBestTour.assign(allBestTours.begin(), allBestTours.begin() + tourSize);
            for (int i = 1; i < size; i++) {
                if (allBestDistances[i] < globalBestDistance) {
                    globalBestDistance = allBestDistances[i];
                    globalBestTour.assign(allBestTours.begin() + i * tourSize, allBestTours.begin() + (i + 1) * tourSize);
                }
            }
        }

        // Broadcast the global best distance and tour to all processes
        MPI_Bcast(&globalBestDistance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(globalBestTour.data(), tourSize, MPI_INT, 0, MPI_COMM_WORLD);

        // Print the best distance (optional)
        // if (rank == 0) {
        //     std::cout << "Generation " << gen << " Best Distance: " << globalBestDistance << std::endl;
        // }

        // Create a new population
        std::vector<std::vector<int>> newLocalPopulation;
        std::uniform_real_distribution<> dis(0.0, 1.0);
        for (int i = 0; i < individualsPerProcess; i++) {
            // Selection
            std::vector<int> parent1 = selectParent(localPopulation, rng);
            std::vector<int> parent2 = selectParent(localPopulation, rng);

            // Crossover
            std::vector<int> offspring = crossover(parent1, parent2, crossoverRate, rng);

            // Mutation
            mutate(offspring, mutationRate, rng);

            // 2-Opt Improvement
            improve2Opt(offspring);

            newLocalPopulation.push_back(offspring);

            // Update local best
            double fitness = calculateFitness(offspring);
            if (fitness < localBestDistance) {
                localBestDistance = fitness;
                localBestTour = offspring;
            }
        }

        // Replace the old population with the new population
        localPopulation = newLocalPopulation;

        // Reset local best for the next generation
        localBestDistance = std::numeric_limits<double>::infinity();
    }

    std::cout << "2\n";

    if (rank == 0) {
        std::cout << "Final Best Distance: " << globalBestDistance << std::endl;
        std::cout << "Final Best Tour: ";
        for (const auto &city : globalBestTour) {
            std::cout << city << " ";
        }
        std::cout << std::endl;
    }

    // Return the best tour and its distance from the root process
    if (rank == 0) {
        return {globalBestTour, globalBestDistance};
    } else {
        return {std::vector<int>(), 0.0};
    }
}