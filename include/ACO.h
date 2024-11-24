#ifndef TSP_ACO_H
#define TSP_ACO_H

#include "reader.h"
#include "solver.h"
#include <cmath>
#include <iostream>
#include <vector>

class ACO : public TSPSolver {
  private:
    // parameters
    const int NUM_ANTS = 50;
    const int NUM_ITERATIONS = 500;
    const double ALPHA = 1.0;
    const double BETA = 2.0;
    const double RHO = 0.2;
    const double Q = 1.0;

    int numCities;

    std::vector<std::vector<double>> pheromones;
    std::vector<int> bestTour;

    // functions
    int selectNextCity(int current, const std::vector<bool> &visited);
    void updatePheromones(const std::vector<std::vector<int>> &allTours, const std::vector<double> &tourLengths);

  public:
    ACO(const std::vector<Point> &pts) : TSPSolver(pts)
    {
        numCities = points.size();
        pheromones.resize(numCities, std::vector<double>(numCities, 1.0));
        std::cout << "ACO constructor" << std::endl;
    }

    // Main solving method
    void solve();

    double getDistance() { return totalDistance; }
};

#endif