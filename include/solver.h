#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

#include "reader.h"
#include <cmath>
#include <vector>

class TSPSolver {
  protected:
    // data
    std::vector<Point> points;
    double totalDistance;

    // functions
    double calculateDistance(const Point &p1, const Point &p2) const;
    double calculateTourDistance(const std::vector<int> &tour) const;

  public:
    TSPSolver(const std::vector<Point> &pts) : points(pts), totalDistance(0) {}
    double getDistance() { return totalDistance; }
    void printTour() const;
};

class TSPAlgorithm {
  private:
    std::vector<Point> points;
    std::vector<int> currentTour;
    double totalDistance;

    // Helper functions
    double calculateDistance(const Point &p1, const Point &p2) const;
    int findNearestNeighbor(int currentCity, const std::vector<bool> &visited) const;
    bool improve2Opt(std::vector<int> &tour, double &distance);
    double calculateTourDistance(const std::vector<int> &tour) const;
    void reverse(std::vector<int> &tour, int start, int end);

  public:
    TSPAlgorithm(const std::vector<Point> &pts) : points(pts), totalDistance(0) {}

    // Main solving methods
    void solveNearestNeighbor();
    void improve2Opt();

    // Getters
    const std::vector<int> &getTour() const { return currentTour; }
    double getTotalDistance() const { return totalDistance; }
    void printTour() const;
};

#endif // TSP_SOLVER_H