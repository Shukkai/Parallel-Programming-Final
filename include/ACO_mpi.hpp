#ifndef TSP_ACO_MPI_H
#define TSP_ACO_MPI_H

#include "ACO.h"
#include "reader.h"
#include "solver.h"
#include <mpi.h>

class ACOMpi : public ACO {
  private:
    // parameters
    const int NUM_ANTS = 100;
    const int NUM_ITERATIONS = 100;
    const double ALPHA = 1.0;
    const double BETA = 2.0;
    const double RHO = 0.2;
    const double Q = 1.0;

  protected:
    std::vector<double> pheromones;
  public:
    ACOMpi(const std::vector<Point> &pts) : ACO(pts) 
    {
        numCities = points.size();
        pheromones.resize(numCities * numCities, 1.0);
        // gen = std::mt19937(42);
    }
    void solve();
    void updatePheromones(const std::vector<int> &allTours);
    int selectNextCity(int current, const std::vector<bool> &visited);
};

#endif