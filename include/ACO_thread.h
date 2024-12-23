#ifndef TSP_ACO_THREAD_H
#define TSP_ACO_THREAD_H

#include "ACO.h"
#include "reader.h"
#include "solver.h"
#include <mutex>
#include <thread>

class ACOThread : public ACO {
  private:
    // parameters
    const int NUM_ANTS = 100;
    const int NUM_ITERATIONS = 100;
    const double ALPHA = 1.0;
    const double BETA = 2.0;
    const double RHO = 0.2;
    const double Q = 1.0;
    int threads_num;

    // thread functions
    std::mutex mtx;
    void worker(int ants_per_thread);

  public:
    ACOThread(const std::vector<Point> &pts, int threads_num = 4) : ACO(pts), threads_num(threads_num) {}
    void solve();
};

#endif