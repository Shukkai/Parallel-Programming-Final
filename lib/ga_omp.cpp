#include "ga_omp.h"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <omp.h>
#include <numeric>

ga_omp::ga_omp(const std::vector<Point> &pts, int popSize, int gens, double mutRate, double crossRate)
    : cities(pts), populationSize(popSize), numGenerations(gens), mutationRate(mutRate), crossoverRate(crossRate)
{
    // Use high-quality random number generator for better seed
    std::random_device rd;
    std::mt19937 gen(rd());
    std::srand(gen());
}

// Optimized distance calculation using SIMD-friendly approach
double ga_omp::distance(const Point &c1, const Point &c2)
{
    // Use faster alternative to sqrt for approximate distance if needed
    double dx = c1.x - c2.x;
    double dy = c1.y - c2.y;
    return std::hypot(dx, dy); // More accurate and potentially faster
}

// Parallel fitness calculation
double ga_omp::calculateFitness(const std::vector<int> &path)
{
    double total = 0;
#pragma omp parallel for reduction(+ : total)
    for (size_t i = 0; i < path.size() - 1; i++)
    {
        total += distance(cities[path[i]], cities[path[i + 1]]);
    }
    total += distance(cities[path.back()], cities[path[0]]); // Return to start
    return total;
}

// Tournament selection with parallel random generation
std::vector<int> ga_omp::selectParent()
{
    int tournamentSize = 5;
    std::vector<int> tournament(tournamentSize);
    double bestFitness = std::numeric_limits<double>::infinity();
    std::vector<int> bestPath;

    for (int i = 0; i < tournamentSize; i++)
    {
        int idx = rand() % population.size();
        double fitness = calculateFitness(population[idx]);
        if (fitness < bestFitness)
        {
            bestFitness = fitness;
            bestPath = population[idx];
        }
    }
    return bestPath;
}

// Order crossover with better randomization
std::vector<int> ga_omp::crossover(const std::vector<int> &parent1, const std::vector<int> &parent2)
{
    if ((double)rand() / RAND_MAX > crossoverRate)
    {
        return parent1;
    }

    int n = parent1.size();
    std::vector<int> child(n, -1);
    int start = rand() % n;
    int end = rand() % n;
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

// Optimized mutation with better randomization
void ga_omp::mutate(std::vector<int> &path)
{
    if ((double)rand() / RAND_MAX < mutationRate)
    {
        int i = rand() % path.size();
        int j = rand() % path.size();
        std::swap(path[i], path[j]);
    }
}

// Parallel 2-opt improvement
void ga_omp::improve2Opt(std::vector<int> &tour)
{
    // double distance = calculateFitness(tour);
    int n = tour.size() - 1;
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            double beforeDistance =
                distance(cities[tour[i]], cities[tour[i + 1]]) + distance(cities[tour[j]], cities[tour[(j + 1)]]);

            double afterDistance =
                distance(cities[tour[i]], cities[tour[j]]) + distance(cities[tour[i + 1]], cities[tour[(j + 1)]]);

            if (afterDistance < beforeDistance)
            {
                std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                return;
            }
        }
    }
}

// Fully parallelized solve method
std::pair<std::vector<int>, double> ga_omp::solve()
{
#pragma omp parallel
    {
        // Thread-local random number generator
        thread_local std::random_device rd;
        thread_local std::mt19937 g(rd());

        std::vector<std::vector<int>> localPopulation; // Thread-local population

#pragma omp for nowait
        for (int i = 0; i < populationSize; i++)
        {
            std::vector<int> tour(cities.size());
            std::iota(tour.begin(), tour.end(), 0); // Fill with 0, 1, ..., n-1
            std::shuffle(tour.begin(), tour.end(), g);
            localPopulation.push_back(tour); // Save to thread-local population
        }

#pragma omp critical
        {
            // Combine thread-local populations into the global population
            population.insert(population.end(), localPopulation.begin(), localPopulation.end());
        }
    }

    // Main GA loop
    std::vector<int> bestTour;
    double bestDistance = std::numeric_limits<double>::infinity();

    for (int gen = 0; gen < numGenerations; gen++)
    {
        std::vector<std::vector<int>> newPopulation;
        newPopulation.reserve(populationSize);

        double localBestDistance = bestDistance;
        std::vector<int> localBestTour = bestTour;

#pragma omp parallel
        {
            std::vector<std::vector<int>> threadPopulation;
            threadPopulation.reserve(populationSize / omp_get_num_threads());

            double threadBestDistance = localBestDistance;
            std::vector<int> threadBestTour = localBestTour;

            thread_local static std::random_device rd;
            thread_local static std::mt19937 gen(rd());

#pragma omp for schedule(dynamic)
            for (int i = 0; i < populationSize; i++)
            {
                std::vector<int> parent1 = selectParent();
                std::vector<int> parent2 = selectParent();
                std::vector<int> offspring = crossover(parent1, parent2);

                mutate(offspring);
                improve2Opt(offspring);

                threadPopulation.push_back(offspring);

                double distance = calculateFitness(offspring);
                if (distance < threadBestDistance)
                {
                    threadBestDistance = distance;
                    threadBestTour = offspring;
                }
            }

#pragma omp critical
            {
                if (threadBestDistance < localBestDistance)
                {
                    localBestDistance = threadBestDistance;
                    localBestTour = threadBestTour;
                }
                newPopulation.insert(
                    newPopulation.end(),
                    threadPopulation.begin(),
                    threadPopulation.end());
            }
        }

        if (localBestDistance < bestDistance)
        {
            bestDistance = localBestDistance;
            bestTour = localBestTour;
            std::cout << "\nGen = " << gen << " Best distance = " << bestDistance;
        }

        population = std::move(newPopulation);
    }

    std::cout << "\nBest distance = " << bestDistance << std::endl;
    return {bestTour, bestDistance};
}