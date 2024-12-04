// ga.cpp
#include "ga_thread.h"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <mutex>
#include <thread>

ga_thread::ga_thread(const std::vector<Point> &pts, int popSize, int gens, double mutRate, double crossRate, int numofthread)
    : cities(pts), populationSize(popSize), numGenerations(gens), mutationRate(mutRate), crossoverRate(crossRate), numofthread(numofthread)
{
    // srand(time(nullptr));
}

double ga_thread::distance(const Point &c1, const Point &c2)
{
    double dx = c1.x - c2.x;
    double dy = c1.y - c2.y;
    return sqrt(dx * dx + dy * dy);
}

double ga_thread::calculateFitness(const std::vector<int> &path)
{
    double total = 0;
    for (size_t i = 0; i < path.size() - 1; i++)
    {
        total += distance(cities[path[i]], cities[path[i + 1]]);
    }
    total += distance(cities[path.back()], cities[path[0]]); // Return to start
    return total;
}

std::vector<int> ga_thread::selectParent(std::mt19937 &rng)
{
    int tournamentSize = 5;
    std::uniform_int_distribution<> dist(0, population.size() - 1);

    std::vector<int> tournament(tournamentSize);
    double bestFitness = std::numeric_limits<double>::infinity();
    std::vector<int> bestPath;

    for (int i = 0; i < tournamentSize; i++)
    {
        // int idx = rand() % population.size();
        int idx = dist(rng);
        double fitness = calculateFitness(population[idx]);
        if (fitness < bestFitness)
        {
            bestFitness = fitness;
            bestPath = population[idx];
        }
    }
    return bestPath;
}

std::vector<int> ga_thread::crossover(const std::vector<int> &parent1, const std::vector<int> &parent2, std::mt19937 &rng)
{
    std::uniform_real_distribution<> dis_real(0.0, 1.0);
    if (dis_real(rng) > crossoverRate)
    {
        return parent1;
    }

    int n = parent1.size();
    std::vector<int> child(n, -1);
    std::uniform_int_distribution<> dist_int(0, n - 1);
    int start = dist_int(rng);
    int end = dist_int(rng);

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

void ga_thread::mutate(std::vector<int> &path, std::mt19937 &rng)
{
    std::uniform_real_distribution<> dis_real(0.0, 1.0);
    if (dis_real(rng) < mutationRate)
    {
        std::uniform_int_distribution<> dist_int(0, path.size() - 1);
        int i = dist_int(rng);
        int j = dist_int(rng);
        std::swap(path[i], path[j]);
    }
}

void ga_thread::improve2Opt(std::vector<int> &tour)
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

std::pair<std::vector<int>, double> ga_thread::solve()
{
    // Initialize population with random tours
    for (int i = 0; i < populationSize; i++)
    {
        std::vector<int> tour(cities.size());
        for (size_t j = 0; j < cities.size(); j++)
            tour[j] = j;

        std::random_device rd;
        std::mt19937 g(rd()); // Mersenne Twister pseudo-random generator
        std::shuffle(tour.begin(), tour.end(), g);
        population.push_back(tour);
    }

    // Best solution tracking
    std::vector<int> bestTour;
    double bestDistance = std::numeric_limits<double>::infinity();

    // Function to process a chunk of the population in a thread
    auto processChunk = [&](int threadId, int startIdx, int endIdx, std::vector<std::vector<int>> &threadPopulation)
    {
        std::random_device rd;
        std::mt19937 rng(rd());
        std::vector<int> localBestTour;
        double localBestDistance = std::numeric_limits<double>::infinity();

        for (int i = startIdx; i < endIdx; i++)
        {
            std::vector<int> parent1 = selectParent(rng);
            std::vector<int> parent2 = selectParent(rng);
            std::vector<int> offspring = crossover(parent1, parent2, rng);
            mutate(offspring, rng);
            improve2Opt(offspring);

            double distance = calculateFitness(offspring);
            if (distance < localBestDistance)
            {
                localBestDistance = distance;
                localBestTour = offspring;
            }

            threadPopulation.push_back(offspring);
        }

        // Update global best solution
        std::lock_guard<std::mutex> lock(bestMutex);
        if (localBestDistance < bestDistance)
        {
            bestDistance = localBestDistance;
            bestTour = localBestTour;
            std::cout << "\nThread " << threadId << " updated best distance: " << bestDistance;
        }
    };

    // Main GA loop
    for (int gen = 0; gen < numGenerations; gen++)
    {
        std::vector<std::vector<int>> newPopulation;

        // Multithreading setup
        int chunkSize = populationSize / numofthread;
        std::vector<std::thread> threads;
        std::vector<std::vector<std::vector<int>>> threadPopulations(numofthread);

        // Create threads
        for (int t = 0; t < numofthread; t++)
        {
            int startIdx = t * chunkSize;
            int endIdx = (t == numofthread - 1) ? populationSize : startIdx + chunkSize;
            threads.emplace_back(processChunk, t, startIdx, endIdx, std::ref(threadPopulations[t]));
        }

        // Join threads
        for (auto &th : threads)
        {
            th.join();
        }

        // Merge thread populations
        for (const auto &threadPop : threadPopulations)
        {
            newPopulation.insert(newPopulation.end(), threadPop.begin(), threadPop.end());
        }

        population = newPopulation;
    }

    std::cout << "\nFinal best distance = " << bestDistance << std::endl;
    return {bestTour, bestDistance};
}
