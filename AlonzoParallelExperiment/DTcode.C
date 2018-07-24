#include <iostream>
#include <cmath>

#include "visit_writer.c"
#include <vector>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <random>

using std::vector;
using std::cerr;
using std::endl;

#define NUM_POINTS 20000


// OUR CONVENTION
// 
//      p1
//     /  \
// e1 /    \ e3
//   /      \
//  p2-------p3
//      e2
//      
// Between p1 and p2 is e1
// Between p2 and p3 is e2
// Between p1 and p3 is e3
// (and the picture could be flipped, rotated, etc.)
//

double * 
PointsGenerator(int numPoints, int dim = 2)
{
    double *array = new double[numPoints*dim];
    srand(time(NULL));   //USE THIS LINE TO SEED RAND
    for (int i = 0 ; i < numPoints ; i++)
    {
        for (int j = 0 ; j < dim ; j++)
        {
            double rand_value = rand() % 1000000 / 1000000.0;
            array[dim*i+j] = rand_value;
        }
    }

    return array;
}

double
sign(double *p1, double *p2, double *p3)
{
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1]);
}

class OneTriangle
{
  public:
    
    double     p1[2]; 
    double     p2[2]; 
    double     p3[2]; 
    OneTriangle  *triangle_across_e1;
    OneTriangle  *triangle_across_e2;
    OneTriangle  *triangle_across_e3;

    bool      ContainsPoint(double x, double y);

    OneTriangle(const OneTriangle &c)  //Copy Constructor
    {
        memcpy(p1, c.p1, sizeof(double) * 2);
        memcpy(p2, c.p2, sizeof(double) * 2);
        memcpy(p3, c.p3, sizeof(double) * 2);
        triangle_across_e1 = c.triangle_across_e1;
        triangle_across_e2 = c.triangle_across_e2;
        triangle_across_e3 = c.triangle_across_e3;
    }
    
    OneTriangle()
    {
        triangle_across_e1 = NULL;
        triangle_across_e2 = NULL;
        triangle_across_e3 = NULL;
    }
};

bool
OneTriangle::ContainsPoint(double x, double y)
{
    double p4[2];
    p4[0] = x;
    p4[1] = y;
    
    bool b1 = sign(p4, p1, p2) < 0.0f;
    bool b2 = sign(p4, p2, p3) < 0.0f;
    bool b3 = sign(p4, p3, p1) < 0.0f;

    return ((b1 == b2) && (b2 == b3));
}

class DelaunayTriangulation
{
  public:
    double *Initialize(double *, int);
    void   AddPoint(double, double);
    bool   CircumcircleCheck(double*, double*, double*, double*);
    void   Verify();
    void   VerifyMeetDC();                                           // verify meet delaunay condition
    bool   SumAngles(double *, double *);                              // method used to verify
    double *FindVectors(int, OneTriangle *);                         // helps verify DC
    void   DelBoundingTri(double *);
    void   WriteOutTriangle(char *filename);
    bool   AltCircumcircleCheck(double*, double*, double*, double*);
    int    WhatEdge(double *, double *, OneTriangle *);
    void   PrintTri(OneTriangle *);
    double *FindBoundingBox(double *, int);
    bool   isCollinear(double, double, double, double, double, double);

  private:
    std::vector<OneTriangle>  triangles;
    double DetHelp(double, double, double, double);
    void EdgeFlip(int, double*, int);
    
};

void DelaunayTriangulation::Verify()
{
    int ncells = triangles.size();
    int iteration = 0;
    int totalFlips = 0;
    int numTrianglesFlipped = 0;
    bool done = false;
    FILE *flipLog = fopen("Logging/flip.log", "w");     //Create file for logging flips per iteration

    while (!done) {
        numTrianglesFlipped = 0;
        for (int j = 0; j < ncells; j++) {   
            if (triangles[j].triangle_across_e1 != NULL) {
                if (AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e1->p2)) {
                    numTrianglesFlipped++; 
                    EdgeFlip(j,triangles[j].triangle_across_e1->p2, 1);
                }
                else if (AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e1->p1)) {
                    numTrianglesFlipped++; 
                    EdgeFlip(j,triangles[j].triangle_across_e1->p1, 1);
                }
                else if (AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e1->p3)) {
                    numTrianglesFlipped++; 
                    EdgeFlip(j,triangles[j].triangle_across_e1->p3, 1);
                }
            }

            if (triangles[j].triangle_across_e2 != NULL) {
                if (AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e2->p1)) { 
                    numTrianglesFlipped++;
                    EdgeFlip(j,triangles[j].triangle_across_e2->p1, 2);
                }
                else if (AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e2->p2)) { 
                    numTrianglesFlipped++;
                    EdgeFlip(j,triangles[j].triangle_across_e2->p2, 2);
                }
                else if (AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e2->p3)) { 
                    numTrianglesFlipped++;
                    EdgeFlip(j,triangles[j].triangle_across_e2->p3, 2);
                }
            } 

            if (triangles[j].triangle_across_e3 != NULL) {
                if(AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e3->p3)) { 
                    numTrianglesFlipped++;
                    EdgeFlip(j, triangles[j].triangle_across_e3->p3, 3);
                }
                else if(AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e3->p2)) { 
                    numTrianglesFlipped++;
                    EdgeFlip(j, triangles[j].triangle_across_e3->p2, 3);
                }
                else if(AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e3->p1)) { 
                    numTrianglesFlipped++;
                    EdgeFlip(j, triangles[j].triangle_across_e3->p1, 3);
                }
            }
        }

        fprintf(flipLog, "%d\t%d\n", iteration, numTrianglesFlipped);

        totalFlips += numTrianglesFlipped;
        done = (numTrianglesFlipped == 0 ? true : false);
        iteration++;
    }

    printf("Iteration count: %d\n", iteration);
    printf("Total flips: %d\n", totalFlips);
    fclose(flipLog);
}

void DelaunayTriangulation::DelBoundingTri(double *bounding_tri) 
{
    /*
      Here is where I delete the first, bounding triangle - update any triangles who have a triangle_across_e*
      that is this bounding triangle. The DT should now be complete.
    */

    int ncells = triangles.size();
    int edge;
    int j;
    double EPSILON = 0.00001f;

    for (j = ncells - 1; j >= 0; j--) { 
        if ((fabs(triangles[j].p1[0] - (bounding_tri[2])) < EPSILON) ||
            (fabs(triangles[j].p2[0] - (bounding_tri[2])) < EPSILON) ||
            (fabs(triangles[j].p3[0] - (bounding_tri[2])) < EPSILON) ||
            (fabs(triangles[j].p1[0] - (bounding_tri[4])) < EPSILON) ||
            (fabs(triangles[j].p2[0] - (bounding_tri[4])) < EPSILON) ||
            (fabs(triangles[j].p3[0] - (bounding_tri[4])) < EPSILON) ||
            (fabs(triangles[j].p1[1] - (bounding_tri[1])) < EPSILON) ||
            (fabs(triangles[j].p2[1] - (bounding_tri[1])) < EPSILON) ||
            (fabs(triangles[j].p3[1] - (bounding_tri[1])) < EPSILON)) 
        {
            if (triangles[j].triangle_across_e1 != NULL) {
                edge = WhatEdge(triangles[j].p1, triangles[j].p2, triangles[j].triangle_across_e1);
                if (edge == 1) 
                    triangles[j].triangle_across_e1->triangle_across_e1 = NULL;
                else if (edge == 2)
                    triangles[j].triangle_across_e1->triangle_across_e2 = NULL;
                else if (edge == 3)
                    triangles[j].triangle_across_e1->triangle_across_e3 = NULL;
            }
            
            if (triangles[j].triangle_across_e2 != NULL){
                edge = WhatEdge(triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e2);
                if (edge == 1)
                    triangles[j].triangle_across_e2->triangle_across_e1 = NULL;
                else if (edge == 2)
                    triangles[j].triangle_across_e2->triangle_across_e2 = NULL;
                else if (edge == 3)
                    triangles[j].triangle_across_e2->triangle_across_e3 = NULL;
            }
            
            if (triangles[j].triangle_across_e3 != NULL) {
                edge = WhatEdge(triangles[j].p3, triangles[j].p1, triangles[j].triangle_across_e3);
                if (edge == 1)
                    triangles[j].triangle_across_e3->triangle_across_e1 = NULL;
                else if (edge == 2)
                    triangles[j].triangle_across_e3->triangle_across_e2 = NULL;
                else if (edge == 3)
                    triangles[j].triangle_across_e3->triangle_across_e3 = NULL;
            }
            
            triangles[j].p1[0] = -100000000.0f;             //Mark triangle for deletion
        }
    }

    for (j = ncells - 1; j >= 0; j--) { 
        if (triangles[j].p1[0] == -100000000.0f) {
            // moves pointers which are pointing to location after element about to be
            //    erased back by 1 to account for automatic vector reallocation to preserve
            //    contiguous-ness
            
            for (int v = 0; v < triangles.size(); v++) {
                 if (triangles[v].triangle_across_e1 > &(triangles[j]))
                     triangles[v].triangle_across_e1 = triangles[v].triangle_across_e1 - 1;
                 if (triangles[v].triangle_across_e2 > &(triangles[j]))
                     triangles[v].triangle_across_e2 = triangles[v].triangle_across_e2 - 1;
                 if (triangles[v].triangle_across_e3 > &(triangles[j]))
                     triangles[v].triangle_across_e3 = triangles[v].triangle_across_e3 - 1; 
             }

            triangles.erase(triangles.begin() + j);
        }
    }

}

void DelaunayTriangulation::EdgeFlip(int j, double* p4_og, int edge)
{
    /*
     Find the points that share an edge with the 4th point inside the circumcircle. Get this info by which if statement above returns 'true'.
     These points should no longer have an edge between them. Therefore, flip that edge to be between the 4th point and the other point in the triangle.
     
			      1								 1
			    / | \						       /   \
			  /   |   \						     /       \
			2     |     4    <-- Does not meet DT requirement          2 _________ 4	This new triangle does meet the DT condition 
			  \   |   /			Flip to instead    --->      \       /
			    \ | /						       \   /
			      3								 3
     Then update points, edges, etc. of affected triangles to keep DS up to date.
    */
    int new_edge;
    double *p4 = new double[2];
    p4[0] = p4_og[0];
    p4[1] = p4_og[1];

    if (edge == 1) {
        
        OneTriangle *t;
        OneTriangle *ot;
        edge = WhatEdge(triangles[j].p1, triangles[j].p2, triangles[j].triangle_across_e1);
        
            new_edge = (edge % 3) + 1;
            if (new_edge == 1)
                t = triangles[j].triangle_across_e1->triangle_across_e1;
            else if (new_edge == 2)
                t = triangles[j].triangle_across_e1->triangle_across_e2;
            else if (new_edge == 3 )
                t = triangles[j].triangle_across_e1->triangle_across_e3;
    
            new_edge = ((edge + 1) % 3) + 1;
            if (new_edge == 1)
                ot = triangles[j].triangle_across_e1->triangle_across_e1;
            else if (new_edge == 2)
                ot = triangles[j].triangle_across_e1->triangle_across_e2;
            else if (new_edge == 3)
                ot = triangles[j].triangle_across_e1->triangle_across_e3;


        triangles[j].triangle_across_e1->p1[0] = p4[0];
        triangles[j].triangle_across_e1->p1[1] = p4[1];
        triangles[j].triangle_across_e1->p2[0] = triangles[j].p2[0];
        triangles[j].triangle_across_e1->p2[1] = triangles[j].p2[1];
        triangles[j].triangle_across_e1->p3[0] = triangles[j].p3[0];
        triangles[j].triangle_across_e1->p3[1] = triangles[j].p3[1];

        //FIX triangle across
        triangles[j].triangle_across_e1->triangle_across_e2 = triangles[j].triangle_across_e2;
        triangles[j].triangle_across_e1->triangle_across_e3 = &(triangles[j]);
        triangles[j].triangle_across_e1->triangle_across_e1 = ot;
       
        triangles[j].triangle_across_e2 = triangles[j].triangle_across_e1;
        triangles[j].triangle_across_e1 = t;

        //triangles[j].p1 and triangles[j].p3 are the same
        triangles[j].p2[0] = p4[0];
        triangles[j].p2[1] = p4[1];

        //Fix t across
        edge = WhatEdge(triangles[j].p1, triangles[j].p2, triangles[j].triangle_across_e1);
        if (edge == 1)
            triangles[j].triangle_across_e1 -> triangle_across_e1 = &(triangles[j]);
        else if (edge == 2)
            triangles[j].triangle_across_e1 -> triangle_across_e2 = &(triangles[j]);
        else if (edge == 3)
            triangles[j].triangle_across_e1 -> triangle_across_e3 = &(triangles[j]);
       
        //Fix x across
        edge = WhatEdge(triangles[j].triangle_across_e2 -> p2, triangles[j].triangle_across_e2 -> p3, triangles[j].triangle_across_e2 -> triangle_across_e2); 
        if (edge == 1)
            triangles[j].triangle_across_e2 -> triangle_across_e2 -> triangle_across_e1 = triangles[j].triangle_across_e2;
        else if (edge == 2)
            triangles[j].triangle_across_e2 -> triangle_across_e2 -> triangle_across_e2 = triangles[j].triangle_across_e2;
        else if (edge == 3)
            triangles[j].triangle_across_e2 -> triangle_across_e2 -> triangle_across_e3 = triangles[j].triangle_across_e2;
    } 
    
    else if (edge == 2) {

        OneTriangle *t;
        OneTriangle *ot;
        edge = WhatEdge(triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e2);
        
            new_edge = ( ( (edge % 3) + 1) % 3 ) + 1;
            if (new_edge == 1)
                t = triangles[j].triangle_across_e2->triangle_across_e1;
            else if (new_edge == 2)
                t = triangles[j].triangle_across_e2->triangle_across_e2;
            else if (new_edge == 3)
                t = triangles[j].triangle_across_e2->triangle_across_e3;

            new_edge = ((edge % 3) + 1);
            if (new_edge == 1)
                ot = triangles[j].triangle_across_e2->triangle_across_e1;
            else if (new_edge == 2)
                ot = triangles[j].triangle_across_e2->triangle_across_e2;
            else if (new_edge == 3)
                ot = triangles[j].triangle_across_e2->triangle_across_e3;

        triangles[j].triangle_across_e2->p1[0] = triangles[j].p1[0];
        triangles[j].triangle_across_e2->p1[1] = triangles[j].p1[1];
        triangles[j].triangle_across_e2->p2[0] = triangles[j].p2[0];
        triangles[j].triangle_across_e2->p2[1] = triangles[j].p2[1];
        triangles[j].triangle_across_e2->p3[0] = p4[0];
        triangles[j].triangle_across_e2->p3[1] = p4[1];
       
        //FIX triangle across
        triangles[j].triangle_across_e2->triangle_across_e1 = triangles[j].triangle_across_e1;
        triangles[j].triangle_across_e2->triangle_across_e3 = &(triangles[j]);
        triangles[j].triangle_across_e2->triangle_across_e2 = ot;
       
        triangles[j].triangle_across_e1 = triangles[j].triangle_across_e2;
        triangles[j].triangle_across_e2 = t; 

        //triangles[j].p1 and triangles[j].p3 are the same
        triangles[j].p2[0] = p4[0];
        triangles[j].p2[1] = p4[1];
       
        edge = WhatEdge(triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e2);
        if (edge == 1)
            triangles[j].triangle_across_e2 -> triangle_across_e1 = &(triangles[j]);
        else if (edge == 2)
            triangles[j].triangle_across_e2 -> triangle_across_e2 = &(triangles[j]);
        else if (edge == 3)
            triangles[j].triangle_across_e2 -> triangle_across_e3 = &(triangles[j]);
       
        edge = WhatEdge(triangles[j].triangle_across_e1 -> p1, triangles[j].triangle_across_e1 -> p2, triangles[j].triangle_across_e1 -> triangle_across_e1); 
        if (edge == 1)
            triangles[j].triangle_across_e1 -> triangle_across_e1 -> triangle_across_e1 = triangles[j].triangle_across_e1;
        else if (edge == 2)
            triangles[j].triangle_across_e1 -> triangle_across_e1 -> triangle_across_e2 = triangles[j].triangle_across_e1;
        else if (edge == 3)
            triangles[j].triangle_across_e1 -> triangle_across_e1 -> triangle_across_e3= triangles[j].triangle_across_e1;
    } 
    else if (edge == 3) {
        OneTriangle *t;
        OneTriangle *ot;
        edge = WhatEdge(triangles[j].p3, triangles[j].p1, triangles[j].triangle_across_e3);
       
            new_edge = (edge % 3) + 1;
            if (new_edge == 1)
                t = triangles[j].triangle_across_e3->triangle_across_e1;
            else if (new_edge == 2)
                t = triangles[j].triangle_across_e3->triangle_across_e2;
            else if (new_edge == 3)
                t = triangles[j].triangle_across_e3->triangle_across_e3;
        
            new_edge = (((edge + 1) % 3) + 1);

            if (new_edge == 1)
                ot = triangles[j].triangle_across_e3->triangle_across_e1;
            else if (new_edge == 2)
                ot = triangles[j].triangle_across_e3->triangle_across_e2;
            else if (new_edge == 3)
                ot = triangles[j].triangle_across_e3->triangle_across_e3;

        triangles[j].triangle_across_e3->p1[0] = triangles[j].p1[0];
        triangles[j].triangle_across_e3->p1[1] = triangles[j].p1[1];
        triangles[j].triangle_across_e3->p2[0] = triangles[j].p2[0];
        triangles[j].triangle_across_e3->p2[1] = triangles[j].p2[1];
        triangles[j].triangle_across_e3->p3[0] = p4[0];
        triangles[j].triangle_across_e3->p3[1] = p4[1];
       
        //FIX triangle across
        triangles[j].triangle_across_e3->triangle_across_e1 = triangles[j].triangle_across_e1;
        triangles[j].triangle_across_e3->triangle_across_e2 = &(triangles[j]);
        triangles[j].triangle_across_e3->triangle_across_e3 = ot;
       
        triangles[j].triangle_across_e1 = triangles[j].triangle_across_e2;
        triangles[j].triangle_across_e2 = t;

        triangles[j].p1[0] = triangles[j].p2[0];
        triangles[j].p1[1] = triangles[j].p2[1];
        triangles[j].p2[0] = triangles[j].p3[0];
        triangles[j].p2[1] = triangles[j].p3[1];
        triangles[j].p3[0] = p4[0];
        triangles[j].p3[1] = p4[1];
      
        //fix x 
        edge = WhatEdge(triangles[j].triangle_across_e3 -> p1, triangles[j].triangle_across_e3 -> p2, triangles[j].triangle_across_e3 -> triangle_across_e1);
        if (edge == 1)
            triangles[j].triangle_across_e3 -> triangle_across_e1 -> triangle_across_e1 = triangles[j].triangle_across_e3;
        else if (edge == 2)
            triangles[j].triangle_across_e3 -> triangle_across_e1 -> triangle_across_e2 = triangles[j].triangle_across_e3;
        else if (edge == 3)
            triangles[j].triangle_across_e3 -> triangle_across_e1 -> triangle_across_e3 = triangles[j].triangle_across_e3;
       
        //fix t
        edge = WhatEdge(triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e2); 
        if (edge == 1)
            triangles[j].triangle_across_e2 -> triangle_across_e1 = &(triangles[j]);
        else if (edge == 2)
            triangles[j].triangle_across_e2 -> triangle_across_e2 = &(triangles[j]);
        else if (edge == 3)
            triangles[j].triangle_across_e2 -> triangle_across_e3 = &(triangles[j]);
    }
    else printf("\n\n\n***edge error!***\n\n\n");/* */
}

//DO NOT EDIT THIS FUNCTION
void DelaunayTriangulation::WriteOutTriangle(char *filename)
{
    int ncells = triangles.size();
    cerr << "NUMBER OF TRIANGLES is " << ncells << endl;

    int *celltypes = new int[ncells];
    for (int i = 0 ; i < ncells ; i++)
        celltypes[i] = VISIT_TRIANGLE;

    int dimensions = 3; // always 3 for VTK
    int vertices_per_cell = 3;
    int npts = ncells*vertices_per_cell*dimensions;
    float *pts = new float[npts];
    int *conn = new int[ncells*vertices_per_cell];
    int offset = 0;
    for (int i = 0 ; i < ncells ; i++)
    {
        pts[offset+0] = triangles[i].p1[0];
        pts[offset+1] = triangles[i].p1[1];
        pts[offset+2] = 0;
        offset += 3;
        pts[offset+0] = triangles[i].p2[0];
        pts[offset+1] = triangles[i].p2[1];
        pts[offset+2] = 0;
        offset += 3;
        pts[offset+0] = triangles[i].p3[0];
        pts[offset+1] = triangles[i].p3[1];
        pts[offset+2] = 0;
        offset += 3;
    }

    for (int i = 0 ; i < 3*ncells ; i++)
    {
        conn[i] = i;
    }
    write_unstructured_mesh(filename, 0, npts/3, pts,
                            ncells, celltypes, conn, 0,
                            NULL, NULL, NULL, NULL);
}
    
double *
DelaunayTriangulation::Initialize(double *bounding_box, int point_count)
{
    double *bounding_tri = new double[6];

    /* Create triangle enclosing bounding box
                      x1, y1 
                       /@&
                      @@.@@
                     @@   @@
                   *@@     #@%
                  @@.       .@@
                 @@           @@
              ,@@&&&&&&&&&&&&&@@#
              @@@@             @@@@
             @@ &@             @& @@
           .@@  &@             @&  &@/
          &@/   &@             @&   /@&
         @@     &@             @&     @@
        @@      &@             @&      @@*
      &@(       &@             @&       (@&
     @@,........&@@@@@@@@@@@@@@@&.........@@
 x2, y2                                    x3, y3 
    */
    OneTriangle ot ;
    ot.p1[0] = (bounding_box[0] + bounding_box[1]) / 2;    
    ot.p1[1] = 2 * bounding_box[3] - bounding_box[2];    
    ot.p2[0] = 2 * bounding_box[0] - ot.p1[0];    
    ot.p2[1] = bounding_box[2];    
    ot.p3[0] = 2 * bounding_box[1] - ot.p1[0];    
    ot.p3[1] = bounding_box[2];    

    bounding_tri[0] = ot.p1[0];
    bounding_tri[1] = ot.p1[1];
    bounding_tri[2] = ot.p2[0];
    bounding_tri[3] = ot.p2[1];
    bounding_tri[4] = ot.p3[0];
    bounding_tri[5] = ot.p3[1];
    
    //Reserve a contiguous portion of memory for that triangles.
    //Note: Could be source of error is point number is very large.
    triangles.reserve(4 * point_count + 1);

    triangles.emplace_back(ot);

    return bounding_tri;
}

void
DelaunayTriangulation::AddPoint(double x1, double y1)
{
    double EPSILON = 0.000000001f;
    for (int i = 0 ; i < triangles.size() ; i++) {
        if (triangles[i].ContainsPoint(x1, y1)) {
            if ((fabs(x1 - triangles[i].p1[0]) < EPSILON && fabs(y1 - triangles[i].p1[1]) < EPSILON) ||
                (fabs(x1 - triangles[i].p2[0]) < EPSILON && fabs(y1 - triangles[i].p2[1]) < EPSILON) ||
                (fabs(x1 - triangles[i].p3[0]) < EPSILON && fabs(y1 - triangles[i].p3[1]) < EPSILON)) {                               //There are redundant Points, just dont add it.
                printf("Redundant Points\n");
                return;
            }
            else if (isCollinear(x1, y1, triangles[i].p1[0], triangles[i].p1[1], triangles[i].p2[0], triangles[i].p2[1])) { //SIDE 1
                
                OneTriangle tri_B, tri_C;
                OneTriangle *Q1, *Q2, *Q3, *Q4, *A, *B, *C, *D;
                int edge, edge_Do;
                double *p1 = new double[2];
                double *p2 = new double[2];
                double *p3 = new double[2];
                double *p4 = new double[2];
                double *p5 = new double[2];


                A = &(triangles[i]);
                D = A -> triangle_across_e1;

                edge_Do = WhatEdge(A -> p1, A -> p2, D);

                
                
                memcpy(p1, A -> p1, sizeof(double) * 2);
                memcpy(p2, A -> p2, sizeof(double) * 2);
                memcpy(p3, A -> p3, sizeof(double) * 2);
                
                p5[0] = x1;
                p5[1] = y1;


                Q2 = A -> triangle_across_e2;
                Q1 = A -> triangle_across_e3;
               
                //Do edge specific assignment
                if (edge_Do == 1) {
                    Q3 = D -> triangle_across_e3;
                    Q4 = D -> triangle_across_e2;

                    memcpy(p4, D -> p3, sizeof(double) * 2);
                }
                else if (edge_Do == 2) {
                    Q3 = D -> triangle_across_e1;
                    Q4 = D -> triangle_across_e3;

                    memcpy(p4, D -> p1, sizeof(double) * 2);
                }
                else if (edge_Do == 3) {
                    Q3 = D -> triangle_across_e2;
                    Q4 = D -> triangle_across_e1;

                    memcpy(p4, D -> p2, sizeof(double) * 2);
                }
                else {
                    printf("Something went wrong, you called what edge on either NULL Triangle or one that is not adjacent\n");
                }
                
                memcpy(A -> p1, p1, sizeof(double) * 2);
                memcpy(A -> p2, p5, sizeof(double) * 2);
                memcpy(A -> p3, p3, sizeof(double) * 2);

                memcpy(tri_B.p1, p5, sizeof(double) * 2);
                memcpy(tri_B.p2, p2, sizeof(double) * 2);
                memcpy(tri_B.p3, p3, sizeof(double) * 2);

                memcpy(tri_C.p1, p5, sizeof(double) * 2);
                memcpy(tri_C.p2, p4, sizeof(double) * 2);
                memcpy(tri_C.p3, p2, sizeof(double) * 2);

                memcpy(D -> p1, p1, sizeof(double) * 2);
                memcpy(D -> p2, p4, sizeof(double) * 2);
                memcpy(D -> p3, p5, sizeof(double) * 2);
               
                triangles.emplace_back(tri_B);
                triangles.emplace_back(tri_C);
                
                int index = triangles.size() - 1;
                B = &(triangles[index - 1]);
                C = &(triangles[index]);

                A -> triangle_across_e1 = D;
                A -> triangle_across_e2 = B;
                A -> triangle_across_e3 = Q1;

                B -> triangle_across_e1 = C;
                B -> triangle_across_e2 = Q2;
                B -> triangle_across_e3 = A;
                
                C -> triangle_across_e1 = D;
                C -> triangle_across_e2 = Q3;
                C -> triangle_across_e3 = B;
                
                D -> triangle_across_e1 = Q4;
                D -> triangle_across_e2 = C;
                D -> triangle_across_e3 = A;

                edge = WhatEdge(p1, p3, Q1);                //Fix Triangle Across for Q1
                if (edge == 1)
                    Q1 -> triangle_across_e1 = A;
                else if (edge == 2)
                    Q1 -> triangle_across_e2 = A;
                else if (edge == 3)
                    Q1 -> triangle_across_e3 = A;


                edge = WhatEdge(p2, p3, Q2);                //Fix Triangle Across for Q2
                if (edge == 1)
                    Q2 -> triangle_across_e1 = B;
                else if (edge == 2)
                    Q2 -> triangle_across_e2 = B;
                else if (edge == 3)
                    Q2 -> triangle_across_e3 = B;

                edge = WhatEdge(p2, p4, Q3);                //Fix Triangle Across for Q3
                if (edge == 1)
                    Q3 -> triangle_across_e1 = C;
                else if (edge == 2)
                    Q3 -> triangle_across_e2 = C;
                else if (edge == 3)
                    Q3 -> triangle_across_e3 = C;

                edge = WhatEdge(p4, p1, Q4);                //Fix Triangle Across for Q4
                if (edge == 1)
                    Q4 -> triangle_across_e1 = D;
                else if (edge == 2)
                    Q4 -> triangle_across_e2 = D;
                else if (edge == 3)
                    Q4 -> triangle_across_e3 = D;

                delete [] p1;
                delete [] p2;
                delete [] p3;
                delete [] p4;
                delete [] p5;

                return;
            }
            else if (isCollinear(x1, y1, triangles[i].p3[0], triangles[i].p3[1], triangles[i].p2[0], triangles[i].p2[1])) { //SIDE 2 

                OneTriangle tri_B, tri_C;
                OneTriangle *Q1, *Q2, *Q3, *Q4, *A, *B, *C, *D;
                int edge, edge_Do;
                double *p1 = new double[2];
                double *p2 = new double[2];
                double *p3 = new double[2];
                double *p4 = new double[2];
                double *p5 = new double[2];


                A = &(triangles[i]);
                D = A -> triangle_across_e2;

                edge_Do = WhatEdge(A -> p2, A -> p3, D);

                
                
                memcpy(p1, A -> p1, sizeof(double) * 2);
                memcpy(p2, A -> p2, sizeof(double) * 2);
                memcpy(p3, A -> p3, sizeof(double) * 2);
                
                p5[0] = x1;
                p5[1] = y1;


                Q1 = A -> triangle_across_e1;
                Q2 = A -> triangle_across_e3;
               
                //Do edge specific assignment
                if (edge_Do == 1) {
                    Q3 = D -> triangle_across_e3;
                    Q4 = D -> triangle_across_e2;

                    memcpy(p4, D -> p3, sizeof(double) * 2);
                }
                else if (edge_Do == 2) {
                    Q3 = D -> triangle_across_e1;
                    Q4 = D -> triangle_across_e3;

                    memcpy(p4, D -> p1, sizeof(double) * 2);
                }
                else if (edge_Do == 3) {
                    Q3 = D -> triangle_across_e2;
                    Q4 = D -> triangle_across_e1;

                    memcpy(p4, D -> p2, sizeof(double) * 2);
                }
                else {
                    printf("Something went wrong, you called what edge on either NULL Triangle or one that is not adjacent\n");
                }
                
                memcpy(A -> p1, p1, sizeof(double) * 2);
                memcpy(A -> p2, p2, sizeof(double) * 2);
                memcpy(A -> p3, p5, sizeof(double) * 2);

                memcpy(tri_B.p1, p1, sizeof(double) * 2);
                memcpy(tri_B.p2, p5, sizeof(double) * 2);
                memcpy(tri_B.p3, p3, sizeof(double) * 2);

                memcpy(tri_C.p1, p5, sizeof(double) * 2);
                memcpy(tri_C.p2, p4, sizeof(double) * 2);
                memcpy(tri_C.p3, p3, sizeof(double) * 2);

                memcpy(D -> p1, p2, sizeof(double) * 2);
                memcpy(D -> p2, p4, sizeof(double) * 2);
                memcpy(D -> p3, p5, sizeof(double) * 2);
               
                triangles.emplace_back(tri_B);
                triangles.emplace_back(tri_C);
                
                int index = triangles.size() - 1;
                B = &(triangles[index - 1]);
                C = &(triangles[index]);

                A -> triangle_across_e1 = Q1;
                A -> triangle_across_e2 = D;
                A -> triangle_across_e3 = B;

                B -> triangle_across_e1 = A;
                B -> triangle_across_e2 = C;
                B -> triangle_across_e3 = Q2;
                
                C -> triangle_across_e1 = D;
                C -> triangle_across_e2 = Q3;
                C -> triangle_across_e3 = B;
                
                D -> triangle_across_e1 = Q4;
                D -> triangle_across_e2 = C;
                D -> triangle_across_e3 = A;

                edge = WhatEdge(p1, p2, Q1);                //Fix Triangle Across for Q1
                if (edge == 1)
                    Q1 -> triangle_across_e1 = A;
                else if (edge == 2)
                    Q1 -> triangle_across_e2 = A;
                else if (edge == 3)
                    Q1 -> triangle_across_e3 = A;


                edge = WhatEdge(p1, p3, Q2);                //Fix Triangle Across for Q2
                if (edge == 1)
                    Q2 -> triangle_across_e1 = B;
                else if (edge == 2)
                    Q2 -> triangle_across_e2 = B;
                else if (edge == 3)
                    Q2 -> triangle_across_e3 = B;

                edge = WhatEdge(p3, p4, Q3);                //Fix Triangle Across for Q3
                if (edge == 1)
                    Q3 -> triangle_across_e1 = C;
                else if (edge == 2)
                    Q3 -> triangle_across_e2 = C;
                else if (edge == 3)
                    Q3 -> triangle_across_e3 = C;

                edge = WhatEdge(p4, p2, Q4);                //Fix Triangle Across for Q4
                if (edge == 1)
                    Q4 -> triangle_across_e1 = D;
                else if (edge == 2)
                    Q4 -> triangle_across_e2 = D;
                else if (edge == 3)
                    Q4 -> triangle_across_e3 = D;

                delete [] p1;
                delete [] p2;
                delete [] p3;
                delete [] p4;
                delete [] p5;
                
                return;
            }
            else if (isCollinear(x1, y1, triangles[i].p1[0], triangles[i].p1[1], triangles[i].p3[0], triangles[i].p3[1])) { //SIDE 3
                
                OneTriangle tri_B, tri_C;
                OneTriangle *Q1, *Q2, *Q3, *Q4, *A, *B, *C, *D;
                int edge, edge_Do;
                double *p1 = new double[2];
                double *p2 = new double[2];
                double *p3 = new double[2];
                double *p4 = new double[2];
                double *p5 = new double[2];


                A = &(triangles[i]);
                D = A -> triangle_across_e3;

                edge_Do = WhatEdge(A -> p1, A -> p3, D);

                
                
                memcpy(p1, A -> p1, sizeof(double) * 2);
                memcpy(p2, A -> p2, sizeof(double) * 2);
                memcpy(p3, A -> p3, sizeof(double) * 2);
                
                p5[0] = x1;
                p5[1] = y1;


                Q1 = A -> triangle_across_e1;
                Q2 = A -> triangle_across_e2;
               
                //Do edge specific assignment
                if (edge_Do == 1) {
                    Q3 = D -> triangle_across_e2;
                    Q4 = D -> triangle_across_e3;

                    memcpy(p4, D -> p3, sizeof(double) * 2);
                }
                else if (edge_Do == 2) {
                    Q3 = D -> triangle_across_e3;
                    Q4 = D -> triangle_across_e1;

                    memcpy(p4, D -> p1, sizeof(double) * 2);
                }
                else if (edge_Do == 3) {
                    Q3 = D -> triangle_across_e1;
                    Q4 = D -> triangle_across_e2;

                    memcpy(p4, D -> p2, sizeof(double) * 2);
                }
                else {
                    printf("Something went wrong, you called what edge on either NULL Triangle or one that is not adjacent\n");
                }
                
                memcpy(A -> p1, p1, sizeof(double) * 2);
                memcpy(A -> p2, p2, sizeof(double) * 2);
                memcpy(A -> p3, p5, sizeof(double) * 2);

                memcpy(tri_B.p1, p2, sizeof(double) * 2);
                memcpy(tri_B.p2, p3, sizeof(double) * 2);
                memcpy(tri_B.p3, p5, sizeof(double) * 2);

                memcpy(tri_C.p1, p5, sizeof(double) * 2);
                memcpy(tri_C.p2, p3, sizeof(double) * 2);
                memcpy(tri_C.p3, p4, sizeof(double) * 2);

                memcpy(D -> p1, p1, sizeof(double) * 2);
                memcpy(D -> p2, p5, sizeof(double) * 2);
                memcpy(D -> p3, p4, sizeof(double) * 2);
               
                triangles.emplace_back(tri_B);
                triangles.emplace_back(tri_C);
                
                int index = triangles.size() - 1;
                B = &(triangles[index - 1]);
                C = &(triangles[index]);

                A -> triangle_across_e1 = Q1;
                A -> triangle_across_e2 = B;
                A -> triangle_across_e3 = D;

                B -> triangle_across_e1 = Q2;
                B -> triangle_across_e2 = C;
                B -> triangle_across_e3 = A;
                
                C -> triangle_across_e1 = B;
                C -> triangle_across_e2 = Q3;
                C -> triangle_across_e3 = D;
                
                D -> triangle_across_e1 = A;
                D -> triangle_across_e2 = C;
                D -> triangle_across_e3 = Q4;

                edge = WhatEdge(p1, p2, Q1);                //Fix Triangle Across for Q1
                if (edge == 1)
                    Q1 -> triangle_across_e1 = A;
                else if (edge == 2)
                    Q1 -> triangle_across_e2 = A;
                else if (edge == 3)
                    Q1 -> triangle_across_e3 = A;


                edge = WhatEdge(p2, p3, Q2);                //Fix Triangle Across for Q2
                if (edge == 1)
                    Q2 -> triangle_across_e1 = B;
                else if (edge == 2)
                    Q2 -> triangle_across_e2 = B;
                else if (edge == 3)
                    Q2 -> triangle_across_e3 = B;

                edge = WhatEdge(p3, p4, Q3);                //Fix Triangle Across for Q3
                if (edge == 1)
                    Q3 -> triangle_across_e1 = C;
                else if (edge == 2)
                    Q3 -> triangle_across_e2 = C;
                else if (edge == 3)
                    Q3 -> triangle_across_e3 = C;

                edge = WhatEdge(p4, p1, Q4);                //Fix Triangle Across for Q4
                if (edge == 1)
                    Q4 -> triangle_across_e1 = D;
                else if (edge == 2)
                    Q4 -> triangle_across_e2 = D;
                else if (edge == 3)
                    Q4 -> triangle_across_e3 = D;

                delete [] p1;
                delete [] p2;
                delete [] p3;
                delete [] p4;
                delete [] p5;

                return;
            }
//
// T0
//      p1
//     /  \
// e1 /    \ e3
//   /      \
//  p2-------p3
//      e2
//      
//      -->
//
//         p1
//         /|\
//        / | \
//  TA   /  |T2\
// e1   /T1 p4  \ e3  TB
//     /   / \   \
//    /  / T3  \  \
//   / /         \ \
//  //             \\
// p2------e2-------p3
//       TC
//
//  Relationships:
//  Triangle T0 gets split into T1, T2, T3
//  T0 had edges e1, e2, e3 with triangles TA, TB, TC
//  T1 will have points: p1, p2, p4 and triangles across e1 is TA, triangle across e2 is T3, and triangle across e3 is T2
//  T2 will have points: p1, p4, p3 and triangles across e1 is T1, triangle across e2 is T3, and triangle across e3 is TB
//  T3 will have points: p4, p2, p3 and triangles across e1 is T1, triangle across e2 is TC, and triangle across e3 is T2
//
            OneTriangle original_triangle = triangles[i];
            OneTriangle *TA = original_triangle.triangle_across_e1;
            OneTriangle *TC = original_triangle.triangle_across_e2; 
            OneTriangle *TB = original_triangle.triangle_across_e3;

            // split triangle i into three triangles
            // note: no edge flipping or Delaunay business.
            // start by replacing triangle in the current list.
            triangles[i].p3[0] = x1;
            triangles[i].p3[1] = y1;
            OneTriangle *T1 = &(triangles[i]);

            // now add two more triangles.
            OneTriangle new_triangle1;
            new_triangle1.p1[0] = x1; 
            new_triangle1.p1[1] = y1; 
            new_triangle1.p2[0] = original_triangle.p2[0]; 
            new_triangle1.p2[1] = original_triangle.p2[1];
            new_triangle1.p3[0] = original_triangle.p3[0]; 
            new_triangle1.p3[1] = original_triangle.p3[1];
            triangles.emplace_back(new_triangle1);
            int index = triangles.size()-1;
            OneTriangle *T3 = &(triangles[index]);

            OneTriangle new_triangle2;
            new_triangle2.p1[0] = original_triangle.p1[0];
            new_triangle2.p1[1] = original_triangle.p1[1];
            new_triangle2.p2[0] = x1;
            new_triangle2.p2[1] = y1;
            new_triangle2.p3[0] = original_triangle.p3[0]; 
            new_triangle2.p3[1] = original_triangle.p3[1]; 
            triangles.emplace_back(new_triangle2);
            OneTriangle *T2 = &(triangles[index+1]);
        
            T1->triangle_across_e1 = TA; 
            T1->triangle_across_e2 = T3;
            T1->triangle_across_e3 = T2;
            T2->triangle_across_e1 = T1;
            T2->triangle_across_e2 = T3;
            T2->triangle_across_e3 = TB;
            T3->triangle_across_e1 = T1;
            T3->triangle_across_e2 = TC; 
            T3->triangle_across_e3 = T2;
            if (TA != NULL) {
                int edge = WhatEdge(original_triangle.p1, original_triangle.p2, TA);
                if (edge == 1) {
                    TA->triangle_across_e1 = T1;
                }
                else if (edge == 2) {
                    TA->triangle_across_e2 = T1;
                }
                else if (edge == 3) {
                    TA->triangle_across_e3 = T1;
                }
            }

            if (TB != NULL) {
                int edge = WhatEdge(original_triangle.p1, original_triangle.p3, TB);
                if (edge == 1) {
                    TB->triangle_across_e1 = T2;
                }
                else if (edge == 2) {
                    TB->triangle_across_e2 = T2;
                }
                else if (edge == 3) {
                    TB->triangle_across_e3 = T2;
                }
            }

            if (TC != NULL) {
                int edge = WhatEdge(original_triangle.p2, original_triangle.p3, TC);
                if (edge == 1) {
                    TC->triangle_across_e1 = T3;
                }
                else if (edge == 2) {
                    TC->triangle_across_e2 = T3;
                }
                else if (edge == 3) {
                    TC->triangle_across_e3 = T3;
                }
            }
            return;
        }
    }
}

bool
DelaunayTriangulation::AltCircumcircleCheck(double *ptA, double *ptB, double *ptC, double *ptD)
{
    double ax_ = ptA[0]-ptD[0];
    double ay_ = ptA[1]-ptD[1];
    double bx_ = ptB[0]-ptD[0];
    double by_ = ptB[1]-ptD[1];
    double cx_ = ptC[0]-ptD[0];
    double cy_ = ptC[1]-ptD[1];


    double ccw = (ptB[0] - ptA[0]) * (ptC[1] - ptA[1]) - (ptC[0] - ptA[0]) * (ptB[1] - ptA[1]);
    if (ccw <= 0.0f) {
        printf("Clockwise\n");
        printf("PtA:\t%f\t%f\n", ptA[0], ptA[1]);
        printf("PtB:\t%f\t%f\n", ptB[0], ptB[1]);
        printf("PtC:\t%f\t%f\n", ptC[0], ptC[1]);
    }

    double result = ((ax_ * ax_ + ay_ * ay_) * (bx_ * cy_ - cx_ * by_) -
                    (bx_ * bx_ + by_ * by_) * (ax_ * cy_ - cx_ * ay_) +
                    (cx_ * cx_ + cy_ * cy_) * (ax_ * by_ - bx_ * ay_));
    return result > 0;
}

//Returns the side of the triangle that is composed of these two points 
int
DelaunayTriangulation::WhatEdge(double *pt1, double *pt2, OneTriangle *tri)
{
    int total = 0;
    double EPSILON = 0.000001f;

    if (tri == NULL) { //Triangle Pointer was NULL, happens for triangles around perimeter
        return 0;
    }

    if ((fabs(pt1[0] - tri->p1[0]) < EPSILON && fabs(pt1[1] - tri->p1[1]) < EPSILON) || (fabs(pt2[0] - tri->p1[0]) < EPSILON && fabs(pt2[1] - tri->p1[1]) < EPSILON)) {
        total += 1;
    }
    if ((fabs(pt1[0] - tri->p2[0]) < EPSILON && fabs(pt1[1] - tri->p2[1]) < EPSILON) || (fabs(pt2[0] - tri->p2[0]) < EPSILON && fabs(pt2[1] - tri->p2[1]) < EPSILON)) {
        total += 2;
    }
    if ((fabs(pt1[0] - tri->p3[0]) < EPSILON && fabs(pt1[1] - tri->p3[1]) < EPSILON) || (fabs(pt2[0] - tri->p3[0]) < EPSILON && fabs(pt2[1] - tri->p3[1]) < EPSILON)) {
        total += 3;
    }

    if (total == 3) {        //Points 1 and 2
        return 1;
    }
    else if (total == 4) {   //Points 1 and 3
        return 3;
    }
    else if (total == 5) {   //Points 2 and 3
        return 2;
    }
    else {
        return 0;
    }
}

//Function for debugging. Prints info about triangle pointed to by *t
void
DelaunayTriangulation::PrintTri(OneTriangle *t)
{
    printf("Triangle at addr:\t%p\n", t );
    printf("P1: %f\t%f\n", t->p1[0], t->p1[1]);
    printf("P2: %f\t%f\n", t->p2[0], t->p2[1]);
    printf("P3: %f\t%f\n", t->p3[0], t->p3[1]);
    printf("Triangle across e1:\t%p\n", t->triangle_across_e1);
    printf("Triangle across e2:\t%p\n", t->triangle_across_e2);
    printf("Triangle across e3:\t%p\n", t->triangle_across_e3);
    printf("**********************************************\n\n");
}

//        creates vectors (outsidevectors) for use in sumangles in 2nd verify function
double *
DelaunayTriangulation::FindVectors(int edgeOfEdgeTri, OneTriangle * overEdge) 
{
    // will have 2 vectors, vector 1 is vectors[0-1], vector 2 is vectors[2-3]
    double * vectors = new double[4];
    for (int i = 0; i < 4; i++)     // init for error checking
        vectors[i] = 999.9;
    if (edgeOfEdgeTri == 1) {
        // create vector 1 (x2-x1, y2-y1)
        vectors[0] = overEdge->p1[0] - overEdge->p3[0];
        vectors[1] = overEdge->p1[1] - overEdge->p3[1];
        // create vector 2 (x3-x1, y3-y1)
        vectors[2] = overEdge->p2[0] - overEdge->p3[0];
        vectors[3] = overEdge->p2[1] - overEdge->p3[1];
    }
    else if (edgeOfEdgeTri == 2) {
        vectors[0] = overEdge->p2[0] - overEdge->p1[0];
        vectors[1] = overEdge->p2[1] - overEdge->p1[1];
        vectors[2] = overEdge->p3[0] - overEdge->p1[0];
        vectors[3] = overEdge->p3[1] - overEdge->p1[1];
    }
    else if (edgeOfEdgeTri == 3) {
        vectors[0] = overEdge->p1[0] - overEdge->p2[0];
        vectors[1] = overEdge->p1[1] - overEdge->p2[1];
        vectors[2] = overEdge->p3[0] - overEdge->p2[0];
        vectors[3] = overEdge->p3[1] - overEdge->p2[1];
    }
    else {
        cerr << "in findvectors, edge is not 1 2 or 3" << endl;
    }
    return vectors;
}

//       sums angles from vectors given in insideV (vectors composed of points in main triangle
//       and outsideV (vectors composed of points in edge triangle)
bool   DelaunayTriangulation::SumAngles(double * insideV, double * outsideV) 
{ 
    double m1 = sqrt((insideV[0] * insideV[0]) + (insideV[1] * insideV[1]));
    double m2 = sqrt((insideV[2] * insideV[2]) + (insideV[3] * insideV[3]));
    double insideRadians = acos((insideV[0]*insideV[2] + insideV[1] * insideV[3]) / (m1 * m2) );
    m1 = sqrt((outsideV[0] * outsideV[0]) + (outsideV[1] * outsideV[1]));
    m2 = sqrt((outsideV[2] * outsideV[2]) + (outsideV[3] * outsideV[3]));
    double outsideRadians = acos((outsideV[0]*outsideV[2] + outsideV[1] * outsideV[3]) / (m1 * m2) );
    if ((insideRadians + outsideRadians) > 3.14159265359) { 
        cerr << "---------\nsum over 180 degrees: " << (insideRadians + outsideRadians) * (180.0/3.14159) << endl;
        cerr << "inside degree: " << (insideRadians * (180.0/3.14159)) << "\toutside degree: " << (outsideRadians * (180.0/3.14159)) << endl;
        cerr << "\ninside: \tvector1 (" << insideV[0] << ", " << insideV[1] << ")";
        cerr << "\tvector2 (" << insideV[2] << ", " << insideV[3] << ")" << endl;
        cerr << "outside:\tvector1 (" << outsideV[0] << ", " << outsideV[1] << ")";
        cerr << "\tvector2 (" << outsideV[2] << ", " << outsideV[3] << ")" << endl;
        return true;
    }
    else 
        return false;
}

//       function to check whether meets delaunay condition - sums angle of opposite 
//       points' vectors (4th point to check and point opposite from evaluated edge)
//       and determines whether sum is over 180 degrees
void   DelaunayTriangulation::VerifyMeetDC() 
{
    int numTriangles = triangles.size();
    // outside vectors are vectors made up of triangle from 4th point, inside vectors from current index triangle
    //     arrays organized as [v1x, v1y, v2x, v2y], where
    //     (v1x, v2y) = 1st vector
    //     (v2x, v2y) = 2nd vector
    double * outsideVectors;
    double * insideVectors = new double[4];
    for (int i = 0; i < numTriangles; i++) {
        if (triangles[i].p1[0] != -100000000.0f) {
            if (triangles[i].triangle_across_e1) {    // check over edge 1
                insideVectors[0] = triangles[i].p1[0] - triangles[i].p3[0];
                insideVectors[1] = triangles[i].p1[1] - triangles[i].p3[1];
                insideVectors[2] = triangles[i].p2[0] - triangles[i].p3[0];
                insideVectors[3] = triangles[i].p2[1] - triangles[i].p3[1];
                int outEdge = WhatEdge(triangles[i].p1, triangles[i].p2, triangles[i].triangle_across_e1);
                outsideVectors = FindVectors(outEdge, triangles[i].triangle_across_e1);
                if (SumAngles(insideVectors, outsideVectors)) {
                    cerr << "SUM > 180 DEGREES OVER EDGE 1 FOR TRI INDEX: " << i << "\n------------\n" << endl;
                    PrintTri(&triangles[i]);
                    PrintTri(triangles[i].triangle_across_e1);
                }
            }
            if (triangles[i].triangle_across_e2) {    // check over edge 2
                insideVectors[0] = triangles[i].p2[0] - triangles[i].p1[0];
                insideVectors[1] = triangles[i].p2[1] - triangles[i].p1[1];
                insideVectors[2] = triangles[i].p3[0] - triangles[i].p1[0];
                insideVectors[3] = triangles[i].p3[1] - triangles[i].p1[1];
                int outEdge = WhatEdge(triangles[i].p2, triangles[i].p3, triangles[i].triangle_across_e2);
                outsideVectors = FindVectors(outEdge, triangles[i].triangle_across_e2);
                if (SumAngles(insideVectors, outsideVectors)) {
                    cerr << "SUM > 180 DEGREES OVER EDGE 2 FOR TRI INDEX: " << i << "\n------------\n" << endl;
                    PrintTri(&triangles[i]);
                    PrintTri(triangles[i].triangle_across_e2);
                }
            }
            if (triangles[i].triangle_across_e3) {   // check over edge 3
                insideVectors[0] = triangles[i].p1[0] - triangles[i].p2[0];
                insideVectors[1] = triangles[i].p1[1] - triangles[i].p2[1];
                insideVectors[2] = triangles[i].p3[0] - triangles[i].p2[0];
                insideVectors[3] = triangles[i].p3[1] - triangles[i].p2[1];
                int outEdge = WhatEdge(triangles[i].p3, triangles[i].p1, triangles[i].triangle_across_e3);
                outsideVectors = FindVectors(outEdge, triangles[i].triangle_across_e3);
                if (SumAngles(insideVectors, outsideVectors)) {
                    cerr << "SUM > 180 DEGREES OVER EDGE 3 FOR TRI INDEX: " << i << "\n-----------\n" << endl;
                    PrintTri(&triangles[i]);
                    PrintTri(triangles[i].triangle_across_e3);
                }
            }
        }
    }
    if (outsideVectors)
        delete [] outsideVectors;
    if (insideVectors)
        delete [] insideVectors;
}

double *
DelaunayTriangulation::FindBoundingBox(double *points, int num_points) 
{
    double *bounding_box = new double[4];
    double x_min = 0.0f;
    double x_max = 0.0f;
    double y_min = 0.0f;
    double y_max = 0.0f;

    int i;

    for (i = 0; i < num_points; i++) {
        if (points[2 * i] < x_min) {
            x_min = points[2 * i];
        }
        if (points[2 * i] > x_max) {
            x_max = points[2 * i];
        }
        if (points[2 * i + 1] < y_min) {
            y_min = points[2 * i + 1];
        }
        if (points[2 * i + 1] > y_max) {
            y_max = points[2 * i + 1];
        }
    }

    /*
        TODO
        Add a condition for determining if their are no points to use.
    */
    
    bounding_box[0] = x_min - 1.0f;
    bounding_box[1] = x_max + 1.0f;
    bounding_box[2] = y_min - 1.0f;
    bounding_box[3] = y_max + 1.0f;
   

    return bounding_box;
}

bool
DelaunayTriangulation::isCollinear(double x1, double y1, double x2, double y2, double x3, double y3) 
{
    double EPSILON = 0.00000001f;
    double slope_a = (y2 - y1) / (x2 - x1);
    double slope_b = (y3 - y2) / (x3 - x2);
    double slope_dif = slope_a - slope_b;


    if (isinf(slope_a) && isinf(slope_b)) {
        printf("The slopes were infinite : %f\n", slope_dif);
        return true;
    }

    return (fabs(slope_dif) < EPSILON);
}

class ParallelDelaunayTriangulation 
{
    public:

    DelaunayTriangulation DT, DT1, DT2, DT3, DT4;

    double *bounding_box;
    double *pts;

    double median_right;
    double median_left;
    double width;
    double height;
    
    void
    Initialize(double *pts) 
    {
        this->pts = pts;
        this->bounding_box = DT.FindBoundingBox(pts, NUM_POINTS);
        this->width = bounding_box[1] - bounding_box[0];
        this->height = bounding_box[3] - bounding_box[2];
        
        SortPoints();
    }

    void
    AddAndVerify()
    {
        
        int i = NUM_POINTS / 2, j;
        double *bb, *bt;
        
        bb = DT1.FindBoundingBox(pts, i + 1000);
        bt = DT1.Initialize(bb, i + 1000);
        for (j = 0; j < i + 1000; j++) {
             DT1.AddPoint(pts[2 * j], pts[2 * j + 1]);
        }
        DT1.Verify();
        DT1.DelBoundingTri(bt);
        DT1.VerifyMeetDC();
        
        
        bb = DT2.FindBoundingBox(pts + i - 1000, NUM_POINTS - i + 1000);
        bt = DT2.Initialize(bb, NUM_POINTS - i + 1000);
        for (j = j - 1000; j < NUM_POINTS; j++) {
             DT2.AddPoint(pts[2 * j], pts[2 * j + 1]);
        }
        DT2.Verify();
        DT2.DelBoundingTri(bt);
        DT2.VerifyMeetDC();
    }

    void
    SortPoints() 
    {
        int i = 1, j;
        double temp_point[2];

        while (i < NUM_POINTS) {
            j = i;
            
            while (j > 0 && compareGreater(pts + 2 * j, pts + 2 * j - 2)) {
                memcpy(temp_point, pts + 2 * j, sizeof(double) * 2);
                memcpy(pts + 2 * j, pts + 2 * j - 2, sizeof(double) * 2);
                memcpy(pts + 2 * j - 2, temp_point, sizeof(double) * 2);
                j--;
            }
            i++;
        }
        
        //for (i = 0; i < NUM_POINTS; i++) {
        //    printf("%f\t%f\n", *(pts + 2 * i), *(pts + 2 * i + 1));
        //}
    }
    bool
    compareGreater(double *p1, double *p2)
    {
        return mapToBBX(p1[0]) + (mapToBBY(p1[1]) * width) < mapToBBX(p2[0]) + (mapToBBY(p2[1]) * width);
    }

    double
    mapToBBX(double input) 
    {
        double input_start = bounding_box[0];
        double input_end = bounding_box[1];
        double output_start = 0;
        double output_end = 100;
        return output_start + ((output_end - output_start) / (input_end - input_start)) * (input - input_start);
    }
    double
    mapToBBY(double input) 
    {
        double input_start = bounding_box[2];
        double input_end = bounding_box[3];
        double output_start = 0;
        double output_end = 100;
        return output_start + ((output_end - output_start) / (input_end - input_start)) * (input - input_start);
    }

};

int main()
{
    double *pts = PointsGenerator(NUM_POINTS, 2);
    DelaunayTriangulation DT;
    ParallelDelaunayTriangulation PDT;

    PDT.Initialize(pts);
    printf("Finished Sorting\n");
    PDT.AddAndVerify();
    
    double *bounding_box = DT.FindBoundingBox(pts, NUM_POINTS);
    double *bounding_tri = DT.Initialize(bounding_box, NUM_POINTS);
    

    for (int i = 0 ; i < NUM_POINTS ; i++) {
        DT.AddPoint(pts[2*i], pts[2*i+1]);
    }
    
    DT.Verify(); 
    DT.DelBoundingTri(bounding_tri);
    DT.VerifyMeetDC();

    char *filename = (char *)"kristi.vtk";
    char *filename1 = (char *)"alonzoA.vtk";
    char *filename2 = (char *)"alonzoB.vtk";
    DT.WriteOutTriangle(filename);
    PDT.DT1.WriteOutTriangle(filename1);
    PDT.DT2.WriteOutTriangle(filename2);

    delete [] pts;
    delete [] bounding_box;
    delete [] bounding_tri;
    return 0;
}
