#ifndef GA_H
#define GA_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
struct City
{
    int id, x, y;
};

class GeneticTSP
{
private:
    std::vector<City> cities;
    std::vector<std::vector<int>> population;
    int populationSize;
    int numGenerations;
    double mutationRate;
    double crossoverRate;

    double distance(const City &c1, const City &c2);
    double calculateFitness(const std::vector<int> &path);
    std::vector<int> selectParent();
    std::vector<int> crossover(const std::vector<int> &parent1, const std::vector<int> &parent2);
    void mutate(std::vector<int> &path);

public:
    GeneticTSP(int popSize = 100, int gens = 1000, double mutRate = 0.01, double crossRate = 0.8);
    bool readfile(const std::string &filename);
    std::vector<int> solve();
};

#endif