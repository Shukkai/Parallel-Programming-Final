#include "reader.h"

bool TSPReader::readFile(const std::string &filename)
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

void TSPReader::printData() const
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

const std::vector<Point> &TSPReader::getPoints() const { return points; }