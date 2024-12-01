#include "ACO_thread.h"

void ACOThread::worker(int ants_per_thread)
{
    for (int ant = 0; ant < ants_per_thread; ant++) {
        std::vector<int> tour;
        double tourLength;

        tour = contructSolution();
        tourLength = calculateTourDistance(tour);

        // mutex
        std::lock_guard<std::mutex> lock(mtx);
        if (tourLength < globalBestDistance) {
            bestTour = tour;
            globalBestDistance = tourLength;
        }
    }
}

void ACOThread::solve()
{
    for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
        int ants_per_thread = NUM_ANTS / threads_num;
        std::vector<std::thread> threads;

        for (int t = 0; t < threads_num; t++) {
            threads.push_back(std::thread(&ACOThread::worker, this, ants_per_thread));
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
