#include "ACO_omp.h"

void ACOOmp::solve()
{
    for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
        std::vector<int> tour;
        double tourLength;

#pragma omp parallel for private(tour, tourLength) schedule(dynamic, 4)
        for (int ant = 0; ant < NUM_ANTS; ant++) {
            tour = contructSolution();
            tourLength = calculateTourDistance(tour);

#pragma omp critical
            {
                if (tourLength < globalBestDistance) {
                    bestTour = tour;
                    globalBestDistance = tourLength;
                    // std::cout << "Iteration " << iter << ": Best tour length = " << globalBestDistance << std::endl;
                }
            }
        }

        tau_max = 1.0 / (RHO * globalBestDistance);
        tau_min =
            tau_max * (1.0 - std::pow(0.05, 1.0 / numCities)) / ((numCities / 2 - 1) * std::pow(0.05, 1.0 / numCities));

        updatePheromones(bestTour);
    }

    totalDistance = globalBestDistance;
}
