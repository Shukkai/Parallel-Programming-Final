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

void ACO::updatePheromones(const std::vector<std::vector<int>> &allTours, const std::vector<double> &tourLengths)
{
    // evaporation
    for (int i = 0; i < numCities; i++) {
        for (int j = 0; j < numCities; j++) {
            pheromones[i][j] *= (1.0 - RHO);
            pheromones[i][j] = std::clamp(pheromones[i][j], tau_min, tau_max);
        }
    }

    // deposit
    for (int i = 0; i < NUM_ANTS; i++) {
        double contribution = Q / tourLengths[i];
        int city1, city2;
        for (int j = 0; j < numCities - 1; j++) {
            city1 = allTours[i][j];
            city2 = allTours[i][j + 1];
            pheromones[city1][city2] += contribution;
            pheromones[city2][city1] += contribution;
        }

        // Add return to start
        city1 = allTours[i][numCities - 1];
        city2 = allTours[i][0];
        pheromones[city1][city2] += contribution;
        pheromones[city2][city1] += contribution;
    }
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
                // distance = calculateTourDistance(tour);
                improved = true;
            }
        }
    }
    return improved;
}

void ACO::solve()
{
    for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
        std::vector<std::vector<int>> allTours(NUM_ANTS);
        std::vector<double> tourLengths(NUM_ANTS, 0.0);

        for (int ant = 0; ant < NUM_ANTS; ant++) {
            allTours[ant] = contructSolution();
            tourLengths[ant] = calculateTourDistance(allTours[ant]);

            if (totalDistance == 0 || tourLengths[ant] < totalDistance) {
                bestTour = allTours[ant];
                totalDistance = tourLengths[ant];
                // std::cout << "Iteration " << iter << ": Best tour length = " << totalDistance << std::endl;
            }
        }

        tau_max = 1.0 / (RHO * totalDistance);
        tau_min =
            tau_max * (1.0 - std::pow(0.05, 1.0 / numCities)) / ((numCities / 2 - 1) * std::pow(0.05, 1.0 / numCities));

        updatePheromones(allTours, tourLengths);
    }
}
