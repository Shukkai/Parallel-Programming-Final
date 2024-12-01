#ifndef TSP_ACO_H
#define TSP_ACO_H

#include "reader.h"
#include "solver.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

class ACO : public TSPSolver {
  private:
    // parameters
    const int NUM_ANTS = 100;
    const int NUM_ITERATIONS = 100;
    const double ALPHA = 1.0;
    const double BETA = 2.0;
    const double RHO = 0.2;
    const double Q = 1.0;

  protected:
    // data
    double tau_max, tau_min;
    int numCities;
    std::vector<std::vector<double>> pheromones;
    std::vector<int> bestTour;
    double globalBestDistance = std::numeric_limits<double>::infinity();

    // rng
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> uniform_dist;

    // functions
    int selectNextCity(int current, const std::vector<bool> &visited);
    void updatePheromones(const std::vector<int> &allTours);
    std::vector<int> contructSolution();
    void _2Opt(std::vector<int> &tour);
    void reverse(std::vector<int> &tour, int start, int end);
    bool improve2Opt(std::vector<int> &tour);

  public:
    ACO(const std::vector<Point> &pts) : TSPSolver(pts)
    {
        numCities = points.size();
        pheromones.resize(numCities, std::vector<double>(numCities, 1.0));
        gen = std::mt19937(42);
        // uniform_dist = std::uniform_real_distribution<double>(0.0, numCities - 1);
    }

    // Main solving method
    void solve();
    const std::vector<int> &getTour() const { return bestTour; }
};

#endif