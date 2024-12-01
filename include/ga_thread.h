#ifndef GA_THREAD_H
#define GA_THREAD_H

#include "reader.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <mutex>
#include <limits>

class ga_thread
{
protected:
    std::vector<Point> cities;
    std::vector<std::vector<int>> population;

    int populationSize;
    int numGenerations;
    double mutationRate;
    double crossoverRate;

    double distance(const Point &c1, const Point &c2);
    double calculateFitness(const std::vector<int> &path);
    std::vector<int> selectParent();
    std::vector<int> crossover(const std::vector<int> &parent1, const std::vector<int> &parent2);
    void mutate(std::vector<int> &path);
    void improve2Opt(std::vector<int> &tour);

    // threading data
    std::mutex bestMutex;
    int numofthread;

public:
    ga_thread(const std::vector<Point> &pts, int popSize = 100, int gens = 1000, double mutRate = 0.01,
              double crossRate = 0.8, int numofthread = 4);

    std::pair<std::vector<int>, double> solve();
};

#endif