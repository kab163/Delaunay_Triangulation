#include "Triangle.hpp"
#include "visit_writer.c"

#include <iostream>
#include <vector>
#include <cmath>

using std::vector;
using std::cerr;
using std::endl;

#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#define EPSILON 0.00001

//MARK: Definition
class Triangulation {
    public:
        //MARK: Properties
        std::vector<Triangle>   triangles;
        int                     num_pts;
        double                  bounding_box[4];
        double                  bounding_tri[6];
        double *                pts;
        
        //MARK: Methods
        void                    WriteOutTriangleVTK(char *filename);
        void                    WriteOutTriangleTGL(char *filename);
        void                    ReadFromTGL(char *filename); 
        void                    Initialize(int num_pts, double *pts);
        void                    AddPoint(double, double);
        void                    DelBoundingTri();
        void                    FindBoundingBox();
        void                    FindBoundingTri();
        void                    AddBoundingTri();
};

//MARK: Declarations

void 
Triangulation::WriteOutTriangleVTK(char *filename) {
    int ncells = triangles.size();
    cerr << "NUMBER OF TRIANGLES is " << ncells << endl;

    int *celltypes = new int[ncells];
    for (int i = 0 ; i < ncells ; i++)
        celltypes[i] = VISIT_TRIANGLE;

    int dimensions = 3;
    int vertices_per_cell = 3;
    int npts = ncells * vertices_per_cell * dimensions;
    float *pts = new float[npts];
    int *conn = new int[ncells * vertices_per_cell];
    int offset = 0;

    for (int i = 0 ; i < ncells ; i++) {
        pts[offset + 0] = triangles[i].p1[0];
        pts[offset + 1] = triangles[i].p1[1];
        pts[offset + 2] = 0;
        offset += 3;
        pts[offset + 0] = triangles[i].p2[0];
        pts[offset + 1] = triangles[i].p2[1];
        pts[offset + 2] = 0;
        offset += 3;
        pts[offset + 0] = triangles[i].p3[0];
        pts[offset + 1] = triangles[i].p3[1];
        pts[offset + 2] = 0;
        offset += 3;
    }

    for (int i = 0 ; i < 3 * ncells ; i++) {
        conn[i] = i;
    }

    write_unstructured_mesh(filename, 0, npts/3, pts, 
                            ncells, celltypes, conn, 0,
                            NULL, NULL, NULL, NULL);
}

void 
Triangulation::Initialize(int num_pts, double *pts) {
    this -> pts = pts;
    this -> num_pts = num_pts;
     
}

void
Triangulation::DelBoundingTri() {
    int ncells = triangles.size();
    int edge;
    int j;

    for (j = ncells - 1; j >= 0; j--) {
        if ((fabs(triangles[j].p1[0] - bounding_tri[2]) < EPSILON) ||
            (fabs(triangles[j].p2[0] - bounding_tri[2]) < EPSILON) ||
            (fabs(triangles[j].p3[0] - bounding_tri[2]) < EPSILON) ||
            (fabs(triangles[j].p1[0] - bounding_tri[4]) < EPSILON) ||
            (fabs(triangles[j].p2[0] - bounding_tri[4]) < EPSILON) ||
            (fabs(triangles[j].p3[0] - bounding_tri[4]) < EPSILON) ||
            (fabs(triangles[j].p1[1] - bounding_tri[1]) < EPSILON) ||
            (fabs(triangles[j].p2[1] - bounding_tri[1]) < EPSILON) ||
            (fabs(triangles[j].p3[1] - bounding_tri[1]) < EPSILON)) {

            if (triangles[j].triangle_across_e1) {
            
            }
        }
    }
}
        
void
Triangulation::FindBoundingBox() {
    int i;

    double x_min = 0.0;
    double x_max = 0.0;
    double y_min = 0.0;
    double y_max = 0.0;

    for (i = 0; i < num_pts; i++) {
        if (pts[2 * i] < x_min)
            x_min = pts[2 * i];
        if (pts[2 * i] > x_max)
            x_max = pts[2 * i];
        if (pts[2 * i + 1] < y_min)
            y_min = pts[2 * i + 1];
        if (pts[2 * i + 1] > y_max)
            y_max = pts[2 * i + 1];
    } 

    bounding_box[0] = x_min - 1.0;
    bounding_box[1] = x_max + 1.0;
    bounding_box[2] = y_min - 1.0;
    bounding_box[3] = y_max + 1.0;
}

void
Triangulation::FindBoundingTri() {
   FindBoundingBox();
   
   bounding_tri[0] = (bounding_box[0] + bounding_box[1]) / 2;
   bounding_tri[1] = 2 * bounding_box[3] - bounding_box[2];
   bounding_tri[2] = 2 * bounding_box[0] - bounding_tri[0];
   bounding_tri[3] = bounding_box[2];
   bounding_tri[4] = 2 * bounding_box[1] - bounding_tri[0];
   bounding_tri[5] = bounding_box[3];
}

void
Triangulation::AddBoundingTri() {
    FindBoundingTri();

    Triangle bt;
    bt.p1[0] = bounding_tri[0];
    bt.p1[1] = bounding_tri[1];
    bt.p2[0] = bounding_tri[2];
    bt.p2[1] = bounding_tri[3];
    bt.p3[0] = bounding_tri[4];
    bt.p3[1] = bounding_tri[5];

    triangles.emplace_back(bt);
}

#endif
