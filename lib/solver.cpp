#include "solver.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

double TSPSolver::calculateDistance(const Point &p1, const Point &p2) const
{
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return std::sqrt(dx * dx + dy * dy);
}

int TSPAlgorithm::findNearestNeighbor(int currentCity, const std::vector<bool> &visited) const
{
    double minDistance = std::numeric_limits<double>::max();
    int nearestCity = -1;

    for (size_t i = 0; i < points.size(); ++i)
    {
        if (!visited[i] && i != static_cast<size_t>(currentCity))
        {
            double distance = calculateDistance(points[currentCity], points[i]);
            if (distance < minDistance)
            {
                minDistance = distance;
                nearestCity = i;
            }
        }
    }
    return nearestCity;
}

void TSPAlgorithm::solveNearestNeighbor()
{
    int n = points.size();
    std::vector<bool> visited(n, false);
    currentTour.clear();
    currentTour.reserve(n);

    // Start from the first city
    int currentCity = 0;
    currentTour.push_back(currentCity);
    visited[currentCity] = true;

    // Find nearest neighbor for each remaining city
    for (int i = 1; i < n; ++i)
    {
        currentCity = findNearestNeighbor(currentCity, visited);
        if (currentCity == -1)
            break;
        currentTour.push_back(currentCity);
        visited[currentCity] = true;
    }

    // Add return to start
    currentTour.push_back(currentTour[0]);

    // Calculate total distance
    totalDistance = calculateTourDistance(currentTour);
}

void TSPAlgorithm::reverse(std::vector<int> &tour, int start, int end)
{
    while (start < end)
    {
        std::swap(tour[start], tour[end]);
        start++;
        end--;
    }
}

double TSPAlgorithm::calculateTourDistance(const std::vector<int> &tour) const
{
    double distance = 0;
    for (size_t i = 0; i < tour.size() - 1; ++i)
    {
        distance += calculateDistance(points[tour[i]], points[tour[i + 1]]);
    }
    return distance;
}

bool TSPAlgorithm::improve2Opt(std::vector<int> &tour, double &distance)
{
    int n = tour.size() - 1; // Don't include last city (same as first)
    bool improved = false;

    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            double beforeDistance = calculateDistance(points[tour[i]], points[tour[i + 1]]) +
                                    calculateDistance(points[tour[j]], points[tour[(j + 1)]]);

            double afterDistance = calculateDistance(points[tour[i]], points[tour[j]]) +
                                   calculateDistance(points[tour[i + 1]], points[tour[(j + 1)]]);

            if (afterDistance < beforeDistance)
            {
                reverse(tour, i + 1, j);
                distance = calculateTourDistance(tour);
                improved = true;
            }
        }
    }
    return improved;
}

void TSPAlgorithm::improve2Opt()
{
    while (improve2Opt(currentTour, totalDistance))
    {
        // Continue until no more improvements can be made
    }
}

void TSPAlgorithm::printTour() const
{
    std::cout << "\nTour path:" << std::endl;
    for (size_t i = 0; i < currentTour.size(); ++i)
    {
        std::cout << points[currentTour[i]].id;
        if (i < currentTour.size() - 1)
        {
            std::cout << " -> ";
        }
    }
    std::cout << "\nTotal distance: " << totalDistance << std::endl;
}