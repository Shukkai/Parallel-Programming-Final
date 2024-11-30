// ga.cpp
#include "ga_cpu.h"
#include <random>
#include <algorithm>
#include <ctime>
#include <limits>
#include <iostream>
#include <string>
#include <cmath>

GeneticTSP::GeneticTSP(int popSize, int gens, double mutRate, double crossRate)
    : populationSize(popSize), numGenerations(gens), mutationRate(mutRate), crossoverRate(crossRate)
{
    srand(time(nullptr));
}

bool GeneticTSP::readfile(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    std::string line;
    bool reading_coords = false;

    while (std::getline(file, line))
    {
        if (reading_coords)
        {
            std::istringstream iss(line);
            City c;
            if (iss >> c.id >> c.x >> c.y)
            {
                cities.push_back(c);
            }
        }
        else
        {
            // if (line.find("NAME") != std::string::npos)
            // {
            //     name = line.substr(line.find(":") + 1);
            //     name = name.substr(name.find_first_not_of(" \t"));
            // }
            // else if (line.find("COMMENT") != std::string::npos)
            // {
            //     comment = line.substr(line.find(":") + 1);
            // }
            // else if (line.find("TYPE") != std::string::npos)
            // {
            //     type = line.substr(line.find(":") + 1);
            // }
            // else if (line.find("DIMENSION") != std::string::npos)
            // {
            //     dimension = std::stoi(line.substr(line.find(":") + 1));
            // }
            // else if (line.find("EDGE_WEIGHT_TYPE") != std::string::npos)
            // {
            //     edge_weight_type = line.substr(line.find(":") + 1);
            // }
            if (line.find("NODE_COORD_SECTION") != std::string::npos)
            {
                reading_coords = true;
            }
        }
    }

    file.close();
    return true;
}

double GeneticTSP::distance(const City &c1, const City &c2)
{
    double dx = c1.x - c2.x;
    double dy = c1.y - c2.y;
    return sqrt(dx * dx + dy * dy);
}

double GeneticTSP::calculateFitness(const std::vector<int> &path)
{
    double total = 0;
    for (size_t i = 0; i < path.size() - 1; i++)
    {
        total += distance(cities[path[i]], cities[path[i + 1]]);
    }
    total += distance(cities[path.back()], cities[path[0]]); // Return to start
    return total;
}

std::vector<int> GeneticTSP::selectParent()
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

std::vector<int> GeneticTSP::crossover(const std::vector<int> &parent1, const std::vector<int> &parent2)
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

void GeneticTSP::mutate(std::vector<int> &path)
{
    if ((double)rand() / RAND_MAX < mutationRate)
    {
        int i = rand() % path.size();
        int j = rand() % path.size();
        std::swap(path[i], path[j]);
    }
}

void GeneticTSP::improve2Opt(std::vector<int> &tour)
{
    // double distance = calculateFitness(tour);
    int n = tour.size() - 1;
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            double beforeDistance =
                distance(cities[tour[i]], cities[tour[i + 1]]) +
                distance(cities[tour[j]], cities[tour[(j + 1)]]);

            double afterDistance =
                distance(cities[tour[i]], cities[tour[j]]) +
                distance(cities[tour[i + 1]], cities[tour[(j + 1)]]);

            if (afterDistance < beforeDistance)
            {
                std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                return;
            }
        }
    }
}

std::vector<int> GeneticTSP::solve()
{
    // Initialize population with random tours
    for (int i = 0; i < populationSize; i++)
    {
        std::vector<int> tour(cities.size());
        for (size_t j = 0; j < cities.size(); j++)
            tour[j] = j;
        std::random_device rd;
        std::mt19937 g(rd()); // Mersenne Twister pseudo-random generator

        // Shuffle the array using std::shuffle
        std::shuffle(tour.begin(), tour.end(), g);
        // std::random_shuffle(tour.begin(), tour.end());
        population.push_back(tour);
    }

    // Main GA loop
    std::vector<int> bestTour;
    double bestDistance = std::numeric_limits<double>::infinity();

    for (int gen = 0; gen < numGenerations; gen++)
    {
        std::vector<std::vector<int>> newPopulation;

        while ((int)newPopulation.size() < populationSize)
        {
            std::vector<int> parent1 = selectParent();
            std::vector<int> parent2 = selectParent();
            std::vector<int> offspring = crossover(parent1, parent2);
            mutate(offspring);

            improve2Opt(offspring);

            newPopulation.push_back(offspring);

            double distance = calculateFitness(offspring);
            if (distance < bestDistance)
            {
                bestDistance = distance;
                bestTour = offspring;
                std::cout << "\nGen = "<<gen<< " Best distance = " << bestDistance;
            }
        }

        population = newPopulation;
    }
    std::cout << "\nBest distance = " << bestDistance<<std::endl;
    return bestTour;
}