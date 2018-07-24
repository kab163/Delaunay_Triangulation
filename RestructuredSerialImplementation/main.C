#include "Triangulation.hpp"

#include <stdio.h>

#define NUM_POINTS 20000

double *
PointsGenerator(int numPoints, int dim = 2) {
    double *array = new double[numPoints*dim];
    //srand(time(NULL));
    
    for(int i = 0; i < numPoints; i++) {
        for (int j = 0 ; j < dim; j++) {
            double rand_value = rand() % 1000000 / 1000000.0;
            array[dim * i + j] = rand_value;
        }
    }

    return array;
}

int main(int argc, char **argv) {
    Triangulation tri;
    double *pts = PointsGenerator(NUM_POINTS);

    tri.Initialize(NUM_POINTS, pts);
    tri.AddBoundingTri();
    tri.AddPoints();
    tri.DelBoundingTri();
    tri.WriteOutTriangleVTK((char *) "lonzo.vtk");

    free(pts);
}
