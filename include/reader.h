#ifndef TSP_READER_H
#define TSP_READER_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Point {
    int id;
    double x;
    double y;
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
    bool readFile(const std::string &filename);
    void printData() const;
    const std::vector<Point> &getPoints() const;
};

#endif // TSP_READER_H