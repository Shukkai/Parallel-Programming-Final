#include "ACO.h"

int ACO::selectNextCity(int current, const std::vector<bool> &visited)
{
    std::vector<double> probabilities(numCities, 0.0);
    double totalProb = 0.0;

    for (int nextCity = 0; nextCity < numCities; nextCity++) {
        if (!visited[nextCity]) {
            probabilities[nextCity] = std::pow(pheromones[current][nextCity], ALPHA) *
                                      std::pow(1.0 / calculateDistance(points[current], points[nextCity]), BETA);
            totalProb += probabilities[nextCity];
        }
    }

    double r = uniform_dist(gen) * totalProb;
    for (int i = 0; i < numCities; i++) {
        if (!visited[i]) {
            r -= probabilities[i];
            if (r <= 0) {
                return i;
            }
        }
    }

    for (int i = 0; i < numCities; i++) {
        if (!visited[i]) {
            return i;
        }
    }

    return -1;
}

void ACO::updatePheromones(const std::vector<int> &allTours)
{
    // evaporation
    for (int i = 0; i < numCities; i++) {
        for (int j = 0; j < numCities; j++) {
            pheromones[i][j] *= (1.0 - RHO);
        }
    }

    int city1, city2;
    double deposit = 1.0 / globalBestDistance;
    for (int i = 0; i < numCities - 1; i++) {
        city1 = allTours[i];
        city2 = allTours[i + 1];
        pheromones[city1][city2] += deposit;
        pheromones[city2][city1] += deposit;
        pheromones[city1][city2] = std::clamp(pheromones[city1][city2], tau_min, tau_max);
        pheromones[city2][city1] = std::clamp(pheromones[city2][city1], tau_min, tau_max);
    }

    // Add return to start
    city1 = allTours[numCities - 1];
    city2 = allTours[0];
    pheromones[city1][city2] += deposit;
    pheromones[city2][city1] += deposit;
    pheromones[city1][city2] = std::clamp(pheromones[city1][city2], tau_min, tau_max);
    pheromones[city2][city1] = std::clamp(pheromones[city2][city1], tau_min, tau_max);
}

std::vector<int> ACO::contructSolution()
{
    std::vector<bool> visited(numCities, false);
    std::vector<int> tour;
    tour.reserve(numCities);

    int current = std::uniform_int_distribution<int>(0, numCities - 1)(gen);
    tour.push_back(current);
    visited[current] = true;

    for (int i = 1; i < numCities; i++) {
        int next = selectNextCity(current, visited);
        tour.push_back(next);
        visited[next] = true;
        current = next;
    }

    // 2-opt
    _2Opt(tour);

    return tour;
}

void ACO::_2Opt(std::vector<int> &tour)
{
    while (improve2Opt(tour)) {
        // Continue until no more improvements can be made
    }
}

void ACO::reverse(std::vector<int> &tour, int start, int end)
{
    while (start < end) {
        std::swap(tour[start], tour[end]);
        start++;
        end--;
    }
}

bool ACO::improve2Opt(std::vector<int> &tour)
{
    int n = tour.size() - 1; // Don't include last city
    bool improved = false;

    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            double beforeDistance = calculateDistance(points[tour[i]], points[tour[i + 1]]) +
                                    calculateDistance(points[tour[j]], points[tour[(j + 1)]]);

            double afterDistance = calculateDistance(points[tour[i]], points[tour[j]]) +
                                   calculateDistance(points[tour[i + 1]], points[tour[(j + 1)]]);

            if (afterDistance < beforeDistance) {
                reverse(tour, i + 1, j);
                improved = true;
            }
        }
    }
    return improved;
}

void ACO::solve()
{
    for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
        std::vector<int> tour;
        double tourLength;

        for (int ant = 0; ant < NUM_ANTS; ant++) {
            tour = contructSolution();
            tourLength = calculateTourDistance(tour);

            if (tourLength < globalBestDistance) {
                bestTour = tour;
                globalBestDistance = tourLength;
                // std::cout << "Iteration " << iter << ": Best tour length = " << globalBestDistance << std::endl;
            }
        }

        tau_max = 1.0 / (RHO * globalBestDistance);
        tau_min =
            tau_max * (1.0 - std::pow(0.05, 1.0 / numCities)) / ((numCities / 2 - 1) * std::pow(0.05, 1.0 / numCities));

        updatePheromones(bestTour);
    }

    totalDistance = globalBestDistance;
}
