#include "ACO_thread.h"

void ACOThread::worker(int ants_per_thread)
{
    std::vector<int> localBestTour;
    double localBestDistance = std::numeric_limits<double>::infinity();

    for (int ant = 0; ant < ants_per_thread; ant++) {
        std::vector<int> tour = contructSolution();
        double tourLength = calculateTourDistance(tour);

        if (tourLength < localBestDistance) {
            localBestTour = tour;
            localBestDistance = tourLength;
        }
    }

    // mutex
    std::lock_guard<std::mutex> lock(mtx);
    if (localBestDistance < globalBestDistance) {
        bestTour = localBestTour;
        globalBestDistance = localBestDistance;
    }
}

void ACOThread::solve()
{
    for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
        int ants_per_thread = NUM_ANTS / threads_num;
        std::vector<std::thread> threads;

        for (int t = 0; t < threads_num; t++) {
            threads.emplace_back(&ACOThread::worker, this, ants_per_thread);
        }

        for (std::thread &t : threads) {
            t.join();
        }

        tau_max = 1.0 / (RHO * globalBestDistance);
        tau_min =
            tau_max * (1.0 - std::pow(0.05, 1.0 / numCities)) / ((numCities / 2 - 1) * std::pow(0.05, 1.0 / numCities));

        updatePheromones(bestTour);
    }

    totalDistance = globalBestDistance;
}
