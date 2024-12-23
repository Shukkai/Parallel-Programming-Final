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

#define NUM_ANTS 50
#define NUM_ITERATIONS 100
#define ALPHA 1.0
#define BETA 2.0
#define RHO 0.2
#define Q 1.0
#define GRID 10
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

__global__ void setup_kernel(curandState *state, unsigned long t)
{
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    curand_init(t, idx, 0, &state[idx]);
}

__device__ float calculateDistanceCuda(const Point &p1, const Point &p2)
{
    float dx = p1.x - p2.x;
    float dy = p1.y - p2.y;
    return sqrtf(dx * dx + dy * dy);
}

__device__ int selectNextCityCuda(curandState *my_curandstate, int idx, float *d_pheromones, const Point *points,
                                  const int numCities, int current, bool *visited)
{
    float *prob = new float[numCities], total_prob = 0.0;
    memset(prob, 0, numCities * sizeof(float));

    for (int nextCity = 0; nextCity < numCities; nextCity++) {
        if (!visited[nextCity]) {
            prob[nextCity] = powf(d_pheromones[current * numCities + nextCity], ALPHA) *
                             powf(1.0 / calculateDistanceCuda(points[current], points[nextCity]), BETA);
            total_prob += prob[nextCity];
        }
    }

    float r = curand_uniform(&(my_curandstate[idx])) * total_prob;
    for (int i = 0; i < numCities; i++) {
        if (!visited[i]) {
            r -= prob[i];
            if (r <= 0) {
                delete[] prob;
                return i;
            }
        }
    }

    for (int i = 0; i < numCities; i++) {
        if (!visited[i]) {
            delete[] prob;
            return i;
        }
    }

    return -1;
}

__device__ void _2OptCuda(int *tour, const Point *points, const int numCities)
{
    bool improved = true;
    while (improved) {
        improved = false;
        for (int i = 0; i < numCities - 1; i++) {
            for (int j = i + 1; j < numCities; j++) {
                float beforeDistance = calculateDistanceCuda(points[tour[i]], points[tour[i + 1]]) +
                                       calculateDistanceCuda(points[tour[j]], points[tour[(j + 1)]]);

                float afterDistance = calculateDistanceCuda(points[tour[i]], points[tour[j]]) +
                                      calculateDistanceCuda(points[tour[i + 1]], points[tour[(j + 1)]]);

                if (afterDistance < beforeDistance) {
                    int start = i + 1, end = j;
                    while (start < end) {
                        int temp = tour[start];
                        tour[start] = tour[end];
                        tour[end] = temp;
                        start++;
                        end--;
                    }
                    improved = true;
                }
            }
        }
    }
}

__device__ float calculateTourDistanceCuda(int *tour, const Point *points, const int numCities)
{
    float distance = 0;
    for (int i = 0; i < numCities - 1; i++) {
        distance += calculateDistanceCuda(points[tour[i]], points[tour[i + 1]]);
    }
    distance += calculateDistanceCuda(points[tour[numCities - 1]], points[tour[0]]);
    return distance;
}

__global__ void contructSolutionCuda(curandState *my_curandstate, int *ants_tour, float *d_pheromones,
                                     const Point *points, const int numCities, float *bestDist, int *lock)
{
    int idx = threadIdx.x + blockDim.x * blockIdx.x;

    int *d_ant = new int[numCities];
    for (int i = 0; i < numCities; i++) {
        d_ant[i] = -1;
    }

    float myrandf = curand_uniform(&(my_curandstate[idx]));
    myrandf *= ((float)(numCities - 1) + 0.999999);
    int current = (int)truncf(myrandf);

    d_ant[0] = current;
    bool *visited = new bool[numCities];
    memset(visited, false, numCities * sizeof(bool));
    visited[current] = true;

    for (int i = 1; i < numCities; i++) {
        int next = selectNextCityCuda(my_curandstate, idx, d_pheromones, points, numCities, current, visited);
        d_ant[i] = next;
        visited[next] = true;
        current = next;
    }

    _2OptCuda(d_ant, points, numCities);
    float tour_len = calculateTourDistanceCuda(d_ant, points, numCities);

    bool leave = true;
    while (leave) {
        if (atomicCAS(lock, 0, 1) == 0) {
            if (tour_len < *bestDist) {
                *bestDist = tour_len;

                // copy ant tour to global memory
                for (int i = 0; i < numCities; i++) {
                    ants_tour[i] = d_ant[i];
                }
            }

            leave = false;
            atomicExch(lock, 0);
        }
        // break;
    }

    // free
    delete[] visited;
    delete[] d_ant;
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

    // functions
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

  public:
    ACOCUDA(const std::vector<Point> &pts) : TSPSolver(pts)
    {
        numCities = points.size();
        pheromones.resize(numCities * numCities, 1.0);
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

        int *ants_tour;
        cudaMalloc(&ants_tour, numCities * sizeof(int));
        cudaMemset(ants_tour, -1, numCities * sizeof(int));

        float *d_pheromones;
        cudaMalloc(&d_pheromones, numCities * numCities * sizeof(float));
        cudaMemcpy(d_pheromones, pheromones.data(), numCities * numCities * sizeof(float), cudaMemcpyHostToDevice);

        Point *d_points;
        cudaMalloc(&d_points, numCities * sizeof(Point));
        cudaMemcpy(d_points, points.data(), numCities * sizeof(Point), cudaMemcpyHostToDevice);

        std::vector<int> tour(numCities);
        float tourLength = globalBestDistance;
        for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
            float *cudaBestDist;
            cudaMalloc(&cudaBestDist, sizeof(float));
            cudaMemcpy(cudaBestDist, &globalBestDistance, sizeof(float), cudaMemcpyHostToDevice);

            int *lock;
            cudaMalloc(&lock, sizeof(int));
            cudaMemset(lock, 0, sizeof(int));

            contructSolutionCuda<<<GRID, BLOCK>>>(d_state, ants_tour, d_pheromones, d_points, numCities, cudaBestDist,
                                                  lock);
            cudaDeviceSynchronize();

            cudaMemcpy(tour.data(), ants_tour, numCities * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(pheromones.data(), d_pheromones, numCities * numCities * sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(&tourLength, cudaBestDist, sizeof(float), cudaMemcpyDeviceToHost);

            if (tourLength < globalBestDistance) {
                bestTour = tour;
                globalBestDistance = tourLength;
            }

            tau_max = 1.0 / (RHO * globalBestDistance);
            tau_min = tau_max * (1.0 - std::pow(0.05, 1.0 / numCities)) /
                      ((numCities / 2 - 1) * std::pow(0.05, 1.0 / numCities));
            updatePheromones(bestTour);
        }

        totalDistance = globalBestDistance;
    }
    const std::vector<int> &getTour() const { return bestTour; }
};

int main(int argc, char *argv[])
{
    std::string filename = (argc > 1) ? argv[1] : "a280.tsp";
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
    std::string res_fname = "res/" + filename + "_aco_cuda.txt";
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
