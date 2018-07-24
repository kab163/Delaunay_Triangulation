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
        void                    AddPoints();
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

    triangles.reserve(4 * num_pts + 1);
}

void
Triangulation::AddPoint(double x, double y) {
    int i, edge;
    int num_triangles = triangles.size();

    for (i = 0; i < num_triangles; i++) {
        if (triangles[i].ContainsPoint(x, y)) {
            if ((fabs(x - triangles[i].p1[0]) < EPSILON && fabs(y - triangles[i].p1[1]) < EPSILON) ||
                (fabs(x - triangles[i].p2[0]) < EPSILON && fabs(y - triangles[i].p2[1]) < EPSILON) ||
                (fabs(x - triangles[i].p3[0]) < EPSILON && fabs(y - triangles[i].p3[1]) < EPSILON)) {
                cerr << "REDUNDANT POINTS - NOT ADDING" << endl;
                return;
            }
            //else if, collinear
            else {
                Triangle original_triangle = triangles[i];
                Triangle *TA = original_triangle.triangle_across_e1;
                Triangle *TB = original_triangle.triangle_across_e2;
                Triangle *TC = original_triangle.triangle_across_e3;

                triangles[i].p3[0] = x;
                triangles[i].p3[1] = y;
                Triangle *T1 = &(triangles[i]);

                Triangle new_triangle1;
                new_triangle1.p1[0] = x;
                new_triangle1.p1[1] = y;
                new_triangle1.p2[0] = original_triangle.p2[0];
                new_triangle1.p2[1] = original_triangle.p2[1];
                new_triangle1.p3[0] = original_triangle.p3[0];
                new_triangle1.p3[1] = original_triangle.p3[1];
                triangles.emplace_back(new_triangle1);
                Triangle *T3 = &(triangles[num_triangles]);

                Triangle new_triangle2;
                new_triangle2.p1[0] = original_triangle.p1[0];
                new_triangle2.p1[1] = original_triangle.p1[1];
                new_triangle2.p2[0] = x;
                new_triangle2.p2[1] = y;
                new_triangle2.p3[0] = original_triangle.p3[0];
                new_triangle2.p3[1] = original_triangle.p3[1];
                triangles.emplace_back(new_triangle2);
                Triangle *T2 = &(triangles[num_triangles + 1]);

                T1 -> triangle_across_e1 = TA;
                T1 -> triangle_across_e2 = T3;
                T1 -> triangle_across_e3 = T2;
                T2 -> triangle_across_e1 = T1;
                T2 -> triangle_across_e2 = T3;
                T2 -> triangle_across_e3 = TB;
                T3 -> triangle_across_e1 = T1;
                T3 -> triangle_across_e2 = TC;
                T3 -> triangle_across_e3 = T2;

                if (TA != NULL) {
                    edge = original_triangle.what_edge_e1;
                    if (edge == 1)
                        TA -> triangle_across_e1 = T1;
                    else if (edge == 2)
                        TA -> triangle_across_e2 = T1;
                    else if (edge == 3)
                        TA -> triangle_across_e3 = T1;
                }
                if (TB != NULL) {
                    edge = original_triangle.what_edge_e3;
                    if (edge == 1)
                        TB -> triangle_across_e1 = T2;
                    else if (edge == 2)
                        TB -> triangle_across_e2 = T2;
                    else if (edge == 3)
                        TB -> triangle_across_e3 = T2;
                }
                if (TC != NULL) {
                    edge = original_triangle.what_edge_e2;
                    if (edge == 1)
                        TC -> triangle_across_e1 = T3;
                    else if (edge == 2)
                        TC -> triangle_across_e2 = T3;
                    else if (edge == 3)
                        TC -> triangle_across_e3 = T3;
                }
            }
            return;
        }
    }
}

void
Triangulation::AddPoints() {
    int i;
    cerr << "Entered ADD POINTS" << endl;
    for (i = 0; i < num_pts; i++) {
        AddPoint(pts[2 * i], pts[2 * i + 1]);
    }
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

            if (triangles[j].triangle_across_e1 != NULL) {
                edge = triangles[j].what_edge_e1;
                if (edge == 1) {
                    triangles[j].triangle_across_e1 -> triangle_across_e1 = NULL;
                    triangles[j].triangle_across_e1 -> what_edge_e1 = 0;
                }
                else if (edge == 2) {
                    triangles[j].triangle_across_e1 -> triangle_across_e2 = NULL;
                    triangles[j].triangle_across_e1 -> what_edge_e2 = 0;
                }
                else if (edge == 3) {
                    triangles[j].triangle_across_e1 -> triangle_across_e3 = NULL;
                    triangles[j].triangle_across_e1 ->what_edge_e3 = 0;
                }
            }
            if (triangles[j].triangle_across_e2 != NULL) {
                edge = triangles[j].what_edge_e2;
                if (edge == 1) {
                    triangles[j].triangle_across_e2 -> triangle_across_e1 = NULL;
                    triangles[j].triangle_across_e2 -> what_edge_e1 = 0;
                }
                else if (edge == 2) {
                    triangles[j].triangle_across_e2 -> triangle_across_e2 = NULL;
                    triangles[j].triangle_across_e2 -> what_edge_e2 = 0;
                }
                else if (edge == 3) {
                    triangles[j].triangle_across_e2 -> triangle_across_e3 = NULL;
                    triangles[j].triangle_across_e2 -> what_edge_e3 = 0;
                }
            }
            if (triangles[j].triangle_across_e3 != NULL) {
                edge = triangles[j].what_edge_e3;
                if (edge == 1) {
                    triangles[j].triangle_across_e3 -> triangle_across_e1 = NULL;
                    triangles[j].triangle_across_e3 -> what_edge_e1 = 0;
                }
                else if (edge == 2) {
                    triangles[j].triangle_across_e3 -> triangle_across_e2 = NULL;
                    triangles[j].triangle_across_e3 -> what_edge_e2 = 0;
                }
                else if (edge == 3) {
                    triangles[j].triangle_across_e3 -> triangle_across_e3 = NULL;
                    triangles[j].triangle_across_e3 -> what_edge_e3 = 0;
                }
            }

            for (int v = 0; v < triangles.size(); v++) {
                if (triangles[v].triangle_across_e1 > &(triangles[j]))
                    --triangles[v].triangle_across_e1;
                if (triangles[v].triangle_across_e2 > &(triangles[j]))
                    --triangles[v].triangle_across_e2;
                if (triangles[v].triangle_across_e3 > &(triangles[j]))
                    --triangles[v].triangle_across_e3;
            }

            triangles.erase(triangles.begin() + j);
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
   bounding_tri[5] = bounding_box[2];
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
