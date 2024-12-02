#include "ACO_mpi.h"

int ACOMpi::selectNextCity(int current, const std::vector<bool> &visited)
{
    std::vector<double> probabilities(numCities, 0.0);
    double totalProb = 0.0;

    for (int nextCity = 0; nextCity < numCities; nextCity++) {
        if (!visited[nextCity]) {
            probabilities[nextCity] = std::pow(pheromones[current * numCities + nextCity], ALPHA) *
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

void ACOMpi::updatePheromones(const std::vector<int> &allTours)
{
    // evaporation
    for (int i = 0; i < numCities; i++) {
        for (int j = 0; j < numCities; j++) {
            pheromones[i * numCities + j] *= (1.0 - RHO);
        }
    }

    int city1, city2;
    double deposit = 1.0 / globalBestDistance;
    for (int i = 0; i < numCities - 1; i++) {
        city1 = allTours[i];
        city2 = allTours[i + 1];
        pheromones[city1 * numCities + city2] += deposit;
        pheromones[city2 * numCities + city1] += deposit;
        pheromones[city1 * numCities + city2] = std::clamp(pheromones[city1 * numCities + city2], tau_min, tau_max);
        pheromones[city2 * numCities + city1] = std::clamp(pheromones[city2 * numCities + city1], tau_min, tau_max);
    }

    // Add return to start
    city1 = allTours[numCities - 1];
    city2 = allTours[0];
    pheromones[city1 * numCities + city2] += deposit;
    pheromones[city2 * numCities + city1] += deposit;
    pheromones[city1 * numCities + city2] = std::clamp(pheromones[city1 * numCities + city2], tau_min, tau_max);
    pheromones[city2 * numCities + city1] = std::clamp(pheromones[city2 * numCities + city1], tau_min, tau_max);
}

void ACOMpi::solve() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    gen = std::mt19937(42 * rank);

    std::cout << "Process " << rank << " of " << size << " is starting." << std::endl;

    int antsPerProcess = NUM_ANTS / size;
    std::vector<int> localBestTour(numCities);
    double localBestDistance = std::numeric_limits<double>::infinity();

    // Initialize global best variables
    std::vector<int> globalBestTour(numCities);
    double globalBestDistance = std::numeric_limits<double>::infinity();

    for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
        // Reset local best for each iteration
        localBestDistance = std::numeric_limits<double>::infinity();
        std::vector<int> localBestTourIter(numCities);

        // Each process constructs tours for its assigned ants
        for (int ant = 0; ant < antsPerProcess; ant++) {
            std::vector<int> tour = contructSolution();
            double tourLength = calculateTourDistance(tour);

            if (tourLength < localBestDistance) {
                localBestDistance = tourLength;
                localBestTourIter = tour;
            }
        }

        // Update the local best tour and distance if a better one is found
        // if (localBestDistance < localBestDistance) { // This condition seems redundant; consider revising
        //     localBestDistance = localBestDistance;
        //     localBestTour = localBestTourIter;
        // }

        // Gather all local best distances at the root process
        std::vector<double> allDistances(size);
        MPI_Gather(&localBestDistance, 1, MPI_DOUBLE, allDistances.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Gather all local best tours at the root process
        std::vector<int> allBestTours;
        if (rank == 0) {
            allBestTours.resize(size * numCities);
        }
        MPI_Gather(localBestTourIter.data(), numCities, MPI_INT,
                  allBestTours.data(), numCities, MPI_INT, 0, MPI_COMM_WORLD);

        // Root process identifies the global best distance and corresponding tour
        if (rank == 0) {
            // globalBestDistance = allDistances[0];
            // globalBestTour = std::vector<int>(allBestTours.begin(), allBestTours.begin() + numCities);
            for (int i = 0; i < size; i++) {
                if (allDistances[i] < globalBestDistance) {
                    globalBestDistance = allDistances[i];
                    globalBestTour.assign(
                        allBestTours.begin() + i * numCities,
                        allBestTours.begin() + (i + 1) * numCities
                    );
                }
            }

            // Update pheromone parameters based on the global best distance
            tau_max = 1.0 / (RHO * globalBestDistance);
            tau_min = tau_max * (1.0 - std::pow(0.05, 1.0 / numCities)) /
                      ((numCities / 2 - 1) * std::pow(0.05, 1.0 / numCities));

            std::cout << "Iteration " << iter << ": Global Best Distance = " << globalBestDistance << std::endl;
        }

        // Broadcast the global best distance, tau_max, and tau_min to all processes
        MPI_Bcast(&globalBestDistance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&tau_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&tau_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Broadcast the global best tour to all processes
        MPI_Bcast(globalBestTour.data(), numCities, MPI_INT, 0, MPI_COMM_WORLD);

        std::cout << "Process " << rank << " received global best distance: " << globalBestDistance << std::endl;

        // Root process updates pheromones based on the global best tour
        if (rank == 0) {
            updatePheromones(globalBestTour);
        }
        
        // std::cout << "QAQ\n";
        // Broadcast the updated pheromones to all processes
        MPI_Bcast(pheromones.data(), numCities * numCities, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // std::cout << "QAQ2\n";

        // Optional: Synchronize processes
        MPI_Barrier(MPI_COMM_WORLD);

        std::cout << "Iteration " << iter << " completed by process " << rank << "." << std::endl;
    }

    // Finalize the results
    if (rank == 0) {
        totalDistance = globalBestDistance;
        std::cout << "Final Best Distance: " << totalDistance << std::endl;
        // Optionally, print the best tour
        std::cout << "Best Tour: ";
        for (const auto& city : globalBestTour) {
            std::cout << city << " ";
        }
        std::cout << std::endl;
    }

    // Finalize MPI
    // MPI_Finalize();

    std::cout << "Process " << rank << " has completed execution." << std::endl;
}