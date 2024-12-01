#include <omp.h>

void ga_omp::solve()
{
    // Initialize population with random tours in parallel
    #pragma omp parallel
    {
        std::random_device rd;
        std::mt19937 g(rd()); // Seed with thread ID to avoid same seeds

        #pragma omp for
        for (int i = 0; i < populationSize; i++) {
            std::vector<int> tour(cities.size());
            for (size_t j = 0; j < cities.size(); j++)
                tour[j] = j;

            // Shuffle the array using std::shuffle
            std::shuffle(tour.begin(), tour.end(), g);
            
            #pragma omp critical
            {
                population.push_back(tour);
            }
        }
    }

    // Main GA loop with parallel sections
    std::vector<int> bestTour;
    double bestDistance = std::numeric_limits<double>::infinity();

    for (int gen = 0; gen < numGenerations; gen++) {
        std::vector<std::vector<int>> newPopulation;
        newPopulation.reserve(populationSize);

        double localBestDistance = bestDistance;
        std::vector<int> localBestTour = bestTour;

        // Parallel generation of offspring
        #pragma omp parallel reduction(min:localBestDistance)
        {
            std::vector<std::vector<int>> localNewPopulation;
            localNewPopulation.reserve(populationSize / omp_get_num_threads() + 1);

            #pragma omp for nowait
            for (int i = 0; i < populationSize; i++) {
                std::vector<int> parent1 = selectParent();
                std::vector<int> parent2 = selectParent();
                std::vector<int> offspring = crossover(parent1, parent2);
                mutate(offspring);

                improve2Opt(offspring);

                double distance = calculateFitness(offspring);

                #pragma omp critical
                {
                    localNewPopulation.push_back(offspring);
                    
                    // Update best tour if found a better solution
                    if (distance < localBestDistance) {
                        localBestDistance = distance;
                        localBestTour = offspring;
                        std::cout << "\nGen = " << gen << " Best distance = " << localBestDistance;
                    }
                }
            }

            // Combine local new populations
            #pragma omp critical
            {
                newPopulation.insert(
                    newPopulation.end(), 
                    localNewPopulation.begin(), 
                    localNewPopulation.end()
                );
            }
        }

        // Update global best tour
        if (localBestDistance < bestDistance) {
            bestDistance = localBestDistance;
            bestTour = localBestTour;
        }

        population = newPopulation;
    }

    std::cout << "\nBest distance = " << bestDistance << std::endl;
    return {bestTour, bestDistance};
}