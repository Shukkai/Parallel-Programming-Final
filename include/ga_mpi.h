#ifndef GA_MPI_H
#define GA_MPI_H

#include "reader.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <mpi.h>
#include <random>

class ga_mpi
{
protected:
  std::vector<Point> cities;
  std::vector<std::vector<int>> population;
  int populationSize;
  int numGenerations;
  double mutationRate;
  double crossoverRate;

  int rank;
  int size;

  double distance(const Point &c1, const Point &c2);
  double calculateFitness(const std::vector<int> &path);
  std::vector<int> selectParent(const std::vector<std::vector<int>> &currentPopulation, std::mt19937 &rng);
  std::vector<int> crossover(const std::vector<int> &parent1, const std::vector<int> &parent2, double crossoverRate, std::mt19937 &rng);
  void mutate(std::vector<int> &path, double mutationRate, std::mt19937 &rng);
  void improve2Opt(std::vector<int> &tour);

public:
  ga_mpi(const std::vector<Point> &pts, int popSize = 100, int gens = 1000, double mutRate = 0.01,
             double crossRate = 0.8);
  ~ga_mpi();

  std::pair<std::vector<int>, double> solve();
};

#endif