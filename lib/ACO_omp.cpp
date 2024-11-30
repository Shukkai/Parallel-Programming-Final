#include "ACO_omp.h"

void ACOOmp::solve()
{
    for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
        std::vector<std::vector<int>> allTours(NUM_ANTS);
        std::vector<double> tourLengths(NUM_ANTS, 0.0);

#pragma omp parallel for schedule(dynamic, 4)
        for (int ant = 0; ant < NUM_ANTS; ant++) {
            allTours[ant] = contructSolution();
            tourLengths[ant] = calculateTourDistance(allTours[ant]);

#pragma omp critical
            {
                if (totalDistance == 0 || tourLengths[ant] < totalDistance) {
                    bestTour = allTours[ant];
                    totalDistance = tourLengths[ant];
                    // std::cout << "Iteration " << iter << ": Best tour length = " << totalDistance << std::endl;
                }
            }
        }

        tau_max = 1.0 / (RHO * totalDistance);
        tau_min =
            tau_max * (1.0 - std::pow(0.05, 1.0 / numCities)) / ((numCities / 2 - 1) * std::pow(0.05, 1.0 / numCities));
        updatePheromones(allTours, tourLengths);
    }
}
