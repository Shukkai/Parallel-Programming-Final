#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <omp.h>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#define NUM_ANTS 10
#define NUM_ITERATIONS 200
#define ALPHA 1.0
#define BETA 2.0
#define RHO 0.2
#define Q 1.0
#define GRID 2
#define BLOCK 5

struct Point {
    int id;
    float x;
    float y;
};

class TSPReader {
  private:
    std::string name;
    std::string comment;
    std::string type;
    int dimension;
    std::string edge_weight_type;
    std::vector<Point> points;

  public:
    bool readFile(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }

        std::string line;
        bool reading_coords = false;

        while (std::getline(file, line)) {
            if (reading_coords) {
                std::istringstream iss(line);
                Point p;
                if (iss >> p.id >> p.x >> p.y) {
                    points.push_back(p);
                }
            }
            else {
                if (line.find("NAME") != std::string::npos) {
                    name = line.substr(line.find(":") + 1);
                    name = name.substr(name.find_first_not_of(" \t"));
                }
                else if (line.find("COMMENT") != std::string::npos) {
                    comment = line.substr(line.find(":") + 1);
                }
                else if (line.find("TYPE") != std::string::npos) {
                    type = line.substr(line.find(":") + 1);
                }
                else if (line.find("DIMENSION") != std::string::npos) {
                    dimension = std::stoi(line.substr(line.find(":") + 1));
                }
                else if (line.find("EDGE_WEIGHT_TYPE") != std::string::npos) {
                    edge_weight_type = line.substr(line.find(":") + 1);
                }
                else if (line.find("NODE_COORD_SECTION") != std::string::npos) {
                    reading_coords = true;
                }
            }
        }

        file.close();
        return true;
    }
    void printData() const
    {
        std::cout << "Name: " << name << std::endl;
        std::cout << "Comment: " << comment << std::endl;
        std::cout << "Type: " << type << std::endl;
        std::cout << "Dimension: " << dimension << std::endl;
        std::cout << "Edge Weight Type: " << edge_weight_type << std::endl;
        std::cout << "\nCoordinates:" << std::endl;
        for (const auto &point : points) {
            std::cout << "Point " << point.id << ": (" << point.x << ", " << point.y << ")" << std::endl;
        }
    }
    const std::vector<Point> &getPoints() const { return points; }
};

class TSPSolver {
  protected:
    // data
    std::vector<Point> points;
    float totalDistance;

    // functions
    float calculateDistance(const Point &p1, const Point &p2) const
    {
        float dx = p1.x - p2.x;
        float dy = p1.y - p2.y;
        return std::sqrt(dx * dx + dy * dy);
    }
    float calculateTourDistance(const std::vector<int> &tour) const
    {
        float distance = 0;
        for (size_t i = 0; i < tour.size() - 1; ++i) {
            distance += calculateDistance(points[tour[i]], points[tour[i + 1]]);
        }
        distance += calculateDistance(points[tour.back()], points[tour.front()]);
        return distance;
    }

  public:
    TSPSolver(const std::vector<Point> &pts) : points(pts), totalDistance(std::numeric_limits<float>::infinity()) {}
    float getDistance() { return totalDistance; }
    // void printTour() const;
};

struct Ant {
    int tour[300] = {-1};
    bool visited[300] = {false};
    int numCities;

    Ant(int numCities) : numCities(numCities) {}
    Ant() {}
};

__global__ void setup_kernel(curandState *state, unsigned long t)
{
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    curand_init(t, idx, 0, &state[idx]);
}

__device__ float calculateDistanceCuda(const Point &p1, const Point &p2)
{
    float dx = p1.x - p2.x;
    float dy = p1.y - p2.y;
    return std::sqrt(dx * dx + dy * dy);
}

__device__ int selectNextCityCuda(curandState *my_curandstate, int idx, Ant *d_ant, float *d_pheromones,
                                  const Point *points, const int numCities, int current)
{
    float *prob = new float[numCities], total_prob = 0.0;
    memset(prob, 0, numCities * sizeof(float));

    for (int nextCity = 0; nextCity < numCities; nextCity++) {
        if (!d_ant[idx].visited[nextCity]) {
            prob[nextCity] = pow(d_pheromones[current * numCities + nextCity], ALPHA) *
                             pow(1.0 / calculateDistanceCuda(points[current], points[nextCity]), BETA);
            total_prob += prob[nextCity];
        }
    }

    float r = curand_uniform(&(my_curandstate[idx])) * total_prob;
    for (int i = 0; i < numCities; i++) {
        if (!d_ant[idx].visited[i]) {
            r -= prob[i];
            if (r <= 0) {
                return i;
            }
        }
    }

    for (int i = 0; i < numCities; i++) {
        if (!d_ant[idx].visited[i]) {
            return i;
        }
    }

    return -1;
}

__global__ void contructSolutionCuda(curandState *my_curandstate, Ant *d_ant, float *d_pheromones, const Point *points,
                                     const int numCities)
{
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    d_ant[idx].numCities = numCities;

    float myrandf = curand_uniform(&(my_curandstate[idx]));
    myrandf *= ((float)(d_ant[idx].numCities - 1) + 0.999999);
    int current = (int)truncf(myrandf);

    d_ant[idx].tour[0] = current;
    d_ant[idx].visited[current] = true;

    // if (idx == 0)
    // printf("%lf\n", d_pheromones[1]);

    for (int i = 1; i < numCities; i++) {
        int next = selectNextCityCuda(my_curandstate, idx, d_ant, d_pheromones, points, numCities, current);
        d_ant[idx].tour[i] = next;
        d_ant[idx].visited[next] = true;
        current = next;
    }

    // printf("Thread %d, Done\n", idx);
    // printf("%d\n", d_ant[idx].tour[0]);
}

__global__ void restartAnts(Ant *d_ant)
{
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    for (int i = 0; i < d_ant[idx].numCities; i++) {
        d_ant[idx].tour[i] = -1;
        d_ant[idx].visited[i] = false;
    }
}

class ACOCUDA : public TSPSolver {
  protected:
    // data
    float tau_max, tau_min;
    int numCities;
    // std::vector<std::vector<float>> pheromones;
    std::vector<float> pheromones;
    std::vector<int> bestTour;
    float globalBestDistance = std::numeric_limits<float>::infinity();

    // rng
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<float> uniform_dist;

    // // functions
    // int selectNextCity(int current, const std::vector<bool> &visited)
    // {
    //     std::vector<float> probabilities(numCities, 0.0);
    //     float totalProb = 0.0;

    //     for (int nextCity = 0; nextCity < numCities; nextCity++) {
    //         if (!visited[nextCity]) {
    //             probabilities[nextCity] = std::pow(pheromones[current][nextCity], ALPHA) *
    //                                       std::pow(1.0 / calculateDistance(points[current], points[nextCity]),
    //                                       BETA);
    //             totalProb += probabilities[nextCity];
    //         }
    //     }

    //     float r = uniform_dist(gen) * totalProb;
    //     for (int i = 0; i < numCities; i++) {
    //         if (!visited[i]) {
    //             r -= probabilities[i];
    //             if (r <= 0) {
    //                 return i;
    //             }
    //         }
    //     }

    //     for (int i = 0; i < numCities; i++) {
    //         if (!visited[i]) {
    //             return i;
    //         }
    //     }

    //     return -1;
    // }
    void updatePheromones(const std::vector<int> &allTours)
    {
        // evaporation
        for (int i = 0; i < numCities * numCities; i++) {
            pheromones[i] *= (1.0 - RHO);
        }

        int city1, city2;
        float deposit = 1.0 / globalBestDistance;
        for (int i = 0; i < numCities - 1; i++) {
            city1 = allTours[i];
            city2 = allTours[i + 1];
            pheromones[city1 * numCities + city2] += deposit;
            pheromones[city2 * numCities + city1] += deposit;
            pheromones[city1 * numCities + city2] = std::clamp(pheromones[city1 * numCities + city2], tau_min, tau_max);
            pheromones[city2 * numCities + city1] = std::clamp(pheromones[city2 * numCities + city1], tau_min, tau_max);
        }

        // Add return to start
        city1 = allTours[numCities - 1];
        city2 = allTours[0];
        pheromones[city1 * numCities + city2] += deposit;
        pheromones[city2 * numCities + city1] += deposit;
        pheromones[city1 * numCities + city2] = std::clamp(pheromones[city1 * numCities + city2], tau_min, tau_max);
        pheromones[city2 * numCities + city1] = std::clamp(pheromones[city2 * numCities + city1], tau_min, tau_max);
    }
    // std::vector<int> contructSolution()
    // {
    //     std::vector<bool> visited(numCities, false);
    //     std::vector<int> tour;
    //     tour.reserve(numCities);

    //     int current = std::uniform_int_distribution<int>(0, numCities - 1)(gen);
    //     tour.push_back(current);
    //     visited[current] = true;

    //     for (int i = 1; i < numCities; i++) {
    //         int next = selectNextCity(current, visited);
    //         tour.push_back(next);
    //         visited[next] = true;
    //         current = next;
    //     }

    //     // 2-opt
    //     _2Opt(tour);

    //     return tour;
    // }
    void _2Opt(std::vector<int> &tour)
    {
        while (improve2Opt(tour)) {
            // Continue until no more improvements can be made
        }
    }
    void reverse(std::vector<int> &tour, int start, int end)
    {
        while (start < end) {
            std::swap(tour[start], tour[end]);
            start++;
            end--;
        }
    }
    bool improve2Opt(std::vector<int> &tour)
    {
        int n = tour.size() - 1; // Don't include last city
        bool improved = false;

        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                float beforeDistance = calculateDistance(points[tour[i]], points[tour[i + 1]]) +
                                       calculateDistance(points[tour[j]], points[tour[(j + 1)]]);

                float afterDistance = calculateDistance(points[tour[i]], points[tour[j]]) +
                                      calculateDistance(points[tour[i + 1]], points[tour[(j + 1)]]);

                if (afterDistance < beforeDistance) {
                    reverse(tour, i + 1, j);
                    improved = true;
                }
            }
        }
        return improved;
    }

  public:
    ACOCUDA(const std::vector<Point> &pts) : TSPSolver(pts)
    {
        numCities = points.size();
        pheromones.resize(numCities * numCities, 1.0);
        // gen = std::mt19937(42);
        // uniform_dist = std::uniform_real_distribution<float>(0.0, numCities - 1);
    }

    // Main solving method
    void solve()
    {
        curandState *d_state;
        time_t t;
        time(&t);
        cudaMalloc(&d_state, NUM_ANTS * sizeof(curandState));
        setup_kernel<<<GRID, BLOCK>>>(d_state, t);
        cudaDeviceSynchronize();

        Ant *d_ants;
        cudaMalloc(&d_ants, NUM_ANTS * sizeof(Ant));

        float *d_pheromones;
        cudaMalloc(&d_pheromones, numCities * numCities * sizeof(float));
        // cudaMemcpy(d_pheromones, pheromones.data(), numCities * numCities * sizeof(float), cudaMemcpyHostToDevice);

        Point *d_points;
        cudaMalloc(&d_points, numCities * sizeof(Point));
        cudaMemcpy(d_points, points.data(), numCities * sizeof(Point), cudaMemcpyHostToDevice);

        std::vector<int> tour(numCities);
        for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
            // printf("Iter: %d, %lf\n", iter, pheromones[1]);
            cudaMemcpy(d_pheromones, pheromones.data(), numCities * numCities * sizeof(float), cudaMemcpyHostToDevice);

            contructSolutionCuda<<<GRID, BLOCK>>>(d_state, d_ants, d_pheromones, d_points, numCities);
            cudaDeviceSynchronize();

            for (int ant = 0; ant < NUM_ANTS; ant++) {
                cudaMemcpy(tour.data(), d_ants[ant].tour, numCities * sizeof(int), cudaMemcpyDeviceToHost);
                cudaMemcpy(pheromones.data(), d_pheromones, numCities * numCities * sizeof(float),
                           cudaMemcpyDeviceToHost);

                _2Opt(tour);
                float tourLength = calculateTourDistance(tour);

                if (tourLength < globalBestDistance) {
                    bestTour = tour;
                    globalBestDistance = tourLength;
                    // std::cout << "Iteration " << iter << ": Best tour length = " << globalBestDistance <<
                    // std::endl;
                }
            }

            // printf("Iter: %d, Tour Length: %f\n", iter, globalBestDistance);

            tau_max = 1.0 / (RHO * globalBestDistance);
            tau_min = tau_max * (1.0 - std::pow(0.05, 1.0 / numCities)) /
                      ((numCities / 2 - 1) * std::pow(0.05, 1.0 / numCities));
            updatePheromones(bestTour);

            restartAnts<<<GRID, BLOCK>>>(d_ants);
        }

        totalDistance = globalBestDistance;
    }
    const std::vector<int> &getTour() const { return bestTour; }
};

int main(int argc, char *argv[])
{
    std::string filename = (argc > 0) ? argv[1] : "a280.tsp";
    TSPReader reader;
    if (!reader.readFile(filename)) {
        std::cerr << "Failed to read TSP file." << std::endl;
        return 1;
    }
    // reader.printData();

    std::vector<int> bestTour;
    float bestDistance = 0.0, elapsedTime = 0.0;
    std::cout << "Solving TSP...\n";

    auto start = std::chrono::high_resolution_clock::now();

    ACOCUDA solver(reader.getPoints());
    solver.solve();
    bestTour = solver.getTour();
    bestDistance = solver.getDistance();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    elapsedTime = duration.count();

    // write results to file
    std::string res_fname = "res/aco_cuda.txt";
    std::ofstream file;
    file.open(res_fname);
    file << "Testcase: " << filename << "\n";
    file << "Time taken: " << elapsedTime << " milliseconds\n";
    file << "Total distance: " << bestDistance << "\n";
    file << "Best tour found:\n";
    for (int city : bestTour) {
        file << city << "\n";
    }

    // print results to console
    std::cout << "Distance: " << bestDistance << "\nSolution found in " << elapsedTime << " milliseconds\n";
    return 0;
}