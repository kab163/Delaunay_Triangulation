#include "Triangle.hpp"
#include <vector>

using std::vector;

#ifndef TRIANGULATION_H
#define TRIANGULATION_H

class Triangulation {
    public:
        std::vector<Triangle> triangles;
        int num_points;

        void WriteOutTriangle(char *filename);
};
    
#endif
