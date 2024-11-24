#include "ACO.h"

int ACO::selectNextCity(int current, const std::vector<bool> &visited)
{
    std::vector<double> probabilities(numCities, 0.0);
    double sum = 0.0;

    for (int i = 0; i < numCities; i++) {
        if (!visited[i]) {
            probabilities[i] = std::pow(pheromones[current][i], ALPHA) *
                               std::pow(1.0 / calculateDistance(points[current], points[i]), BETA);
            sum += probabilities[i];
        }
    }

    double r = static_cast<double>(rand()) / RAND_MAX * sum;
    double total = 0.0;

    for (int i = 0; i < numCities; i++) {
        if (!visited[i]) {
            total += probabilities[i];
            if (total >= r) {
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

void ACO::solve()
{
    srand(time(0));

    for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
        std::vector<std::vector<int>> allTours(NUM_ANTS);
        std::vector<double> tourLengths(NUM_ANTS, 0.0);

        for (int ant = 0; ant < NUM_ANTS; ant++) {
            std::vector<bool> visited(numCities, false);
            int currentCity = rand() % numCities;
            allTours[ant].push_back(currentCity);
            visited[currentCity] = true;

            for (int step = 1; step < numCities; step++) {
                int nextCity = selectNextCity(currentCity, visited);
                allTours[ant].push_back(nextCity);
                visited[nextCity] = true;
                tourLengths[ant] += calculateDistance(points[currentCity], points[nextCity]);
                currentCity = nextCity;
            }

            // Add return to start
            tourLengths[ant] += calculateDistance(points[allTours[ant][numCities - 1]], points[allTours[ant][0]]);

            if (totalDistance == 0 || tourLengths[ant] < totalDistance) {
                bestTour = allTours[ant];
                totalDistance = tourLengths[ant];
                std::cout << "Iteration " << iter << ": Best tour length = " << totalDistance << std::endl;
            }
        }

        updatePheromones(allTours, tourLengths);
    }
}