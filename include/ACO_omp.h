#ifndef TSP_ACO_OMP_H
#define TSP_ACO_OMP_H

#include "ACO.h"
#include "reader.h"
#include "solver.h"
#include <omp.h>

class ACOOmp : public ACO {
  private:
    // parameters
    const int NUM_ANTS = 100;
    const int NUM_ITERATIONS = 100;
    const double ALPHA = 1.0;
    const double BETA = 2.0;
    const double RHO = 0.2;
    const double Q = 1.0;

  public:
    ACOOmp(const std::vector<Point> &pts) : ACO(pts) {}
    void solve();
};

#endif