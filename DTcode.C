#include <iostream>
#include <cmath>
#include "visit_writer.c"
#include <vector>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <random>

#define NUM_POINTS 20000
using std::vector;
using std::cerr;
using std::endl;

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

//DO NOT EDIT THIS FUNCTION
float * 
PointsGenerator(int numPoints, int dim = 2)
{
    float *array = new float[numPoints*dim];
    for (int i = 0 ; i < numPoints ; i++)
    {
        for (int j = 0 ; j < dim ; j++)
        {
            float rand_value = rand() % 1000000 / 1000000.0;
            array[dim*i+j] = rand_value;
        }
    }

    return array;
}

//DO NOT EDIT THIS FUNCTION
bool IsOnSameSide(float *endPoint1, float *endPoint2, 
                  float *referencePoint, float *newPoint)
{
    // see: http://doubleroot.in/lessons/straight-line/position-of-a-point-relative-to-a-line/#.Wt5H7ZPwalM
    float m, b;
    // need to solve equation y = mx + b for endPoint1 
    // and endPoint2.

    if (endPoint1[0] == endPoint2[0])
    {
        // infinite slope ... fail
        return false;
    }
    m = (endPoint2[1] - endPoint1[1])/(endPoint2[0] - endPoint1[0]);
    // y = mx+b
    // a'x+b'y+c' = 0
    // mx-y+b = 0;
    // a' = m, b' = -1, c' = b
    b = endPoint2[1]-m*endPoint2[0];
    float a_formula = m;
    float b_formula = -1;
    float c_formula = b;

    float val1 = referencePoint[0]*a_formula + referencePoint[1]*b_formula + c_formula;
    float val2 = newPoint[0]*a_formula + newPoint[1]*b_formula + c_formula;

    float product = val1*val2;
    return (product < 0 ? false : true);
}

class OneTriangle
{
  public:
    
    float     p1[2]; 
    float     p2[2]; 
    float     p3[2]; 
    OneTriangle  *ta_e1;           // ta is abbreviation of ta
    OneTriangle  *ta_e2;
    OneTriangle  *ta_e3;

    bool      ContainsPoint(float x, float y);
    
    OneTriangle(const OneTriangle &c)  //Copy Constructor
    {
        memcpy(p1, c.p1, sizeof(float) * 2);
        memcpy(p2, c.p2, sizeof(float) * 2);
        memcpy(p3, c.p3, sizeof(float) * 2);
        ta_e1 = c.ta_e1;
        ta_e2 = c.ta_e2;
        ta_e3 = c.ta_e3;
    }
    
    OneTriangle()
    {
        ta_e1 = NULL;
        ta_e2 = NULL;
        ta_e3 = NULL;
    }
};

//DO NOT EDIT THIS FUNCTION
bool
OneTriangle::ContainsPoint(float x, float y)
{
    float p4[2];
    p4[0] = x;
    p4[1] = y;
    bool p3_and_p4 = IsOnSameSide(p1, p2, p3, p4);
    bool p1_and_p4 = IsOnSameSide(p3, p2, p1, p4);
    bool p2_and_p4 = IsOnSameSide(p3, p1, p2, p4);
    if (p3_and_p4 && p1_and_p4 && p2_and_p4)
        return true;
    return false;
}

class DelaunayTriangulation
{
  public:
    float *Initialize(float *, int);
    float *FindBoundingBox(float *);                             // alonzo added
    float *FindVectors(int, OneTriangle *);                      // john added
    bool   CircumcircleCheck(float*, float*, float*, float*);
    bool   AltCircumcircleCheck(float*, float*, float*, float*); // alonzo added
    bool   SumAngles(float *, float *);                          // john added
    void   AddPoint(float, float);
    void   Verify();
    void   VerifyMeetDC();                                       // john added
    void   DelBoundingTri(float *);
    void   WriteOutTriangle(char *filename);
    void   PrintTri(OneTriangle *);                              // alonzo added
    void   CopyPoint(float *, float *);                          // john added
    int    WhatEdge(float *, float *, OneTriangle *);            // alonzo added

  private:
    std::vector<OneTriangle>  triangles;
    float DetHelp(float, float, float, float);
    void EdgeFlip(int, float*, int);
    
};

void DelaunayTriangulation::Verify()
{
    int ncells = triangles.size();
    int iteration = 0;
    int totalFlips = 0;
    int numTrianglesFlipped = 0;
    bool done = false;
    FILE *flipLog = fopen("flip.log", "w");     //Create file for logging flips per iteration

    while (!done) {
        numTrianglesFlipped = 0;
        for (int j = 0; j < ncells; j++) {   
            if (triangles[j].ta_e1 != NULL) {
                if (AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].ta_e1->p2)) {
                    numTrianglesFlipped++; 
                    EdgeFlip(j,triangles[j].ta_e1->p2, 1);
                }
                else if (AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].ta_e1->p1)) {
                    numTrianglesFlipped++; 
                    EdgeFlip(j,triangles[j].ta_e1->p1, 1);
                }
                else if (AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].ta_e1->p3)) {
                    numTrianglesFlipped++; 
                    EdgeFlip(j,triangles[j].ta_e1->p3, 1);
                }
            }

            if (triangles[j].ta_e2 != NULL) {
                if (AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].ta_e2->p1)) { 
                    numTrianglesFlipped++;
                    EdgeFlip(j,triangles[j].ta_e2->p1, 2);
                }
                else if (AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].ta_e2->p2)) { 
                    numTrianglesFlipped++;
                    EdgeFlip(j,triangles[j].ta_e2->p2, 2);
                }
                else if (AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].ta_e2->p3)) { 
                    numTrianglesFlipped++;
                    EdgeFlip(j,triangles[j].ta_e2->p3, 2);
                }
            } 

            if (triangles[j].ta_e3 != NULL) {
                if(AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].ta_e3->p3)) { 
                    numTrianglesFlipped++;
                    EdgeFlip(j, triangles[j].ta_e3->p3, 3);
                }
                else if(AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].ta_e3->p2)) { 
                    numTrianglesFlipped++;
                    EdgeFlip(j, triangles[j].ta_e3->p2, 3);
                }
                else if(AltCircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].ta_e3->p1)) { 
                    numTrianglesFlipped++;
                    EdgeFlip(j, triangles[j].ta_e3->p1, 3);
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

void DelaunayTriangulation::DelBoundingTri(float *bounding_tri) 
{
    /*
      Here is where I delete the first, bounding triangle - update any triangles who have a ta_e* that is this bounding triangle. The DT should now be complete.
    */
    int ncells = triangles.size(), edge, j;
    float EPSILON = 0.00001f;

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
            if (triangles[j].ta_e1 != NULL) {
                edge = WhatEdge(triangles[j].p1, triangles[j].p2, triangles[j].ta_e1);
                if (edge == 1) 
                    triangles[j].ta_e1->ta_e1 = NULL;
                else if (edge == 2)
                    triangles[j].ta_e1->ta_e2 = NULL;
                else if (edge == 3)
                    triangles[j].ta_e1->ta_e3 = NULL;
            }
            
            if (triangles[j].ta_e2 != NULL){
                edge = WhatEdge(triangles[j].p2, triangles[j].p3, triangles[j].ta_e2);
                if (edge == 1)
                    triangles[j].ta_e2->ta_e1 = NULL;
                else if (edge == 2)
                    triangles[j].ta_e2->ta_e2 = NULL;
                else if (edge == 3)
                    triangles[j].ta_e2->ta_e3 = NULL;
            }
            
            if (triangles[j].ta_e3 != NULL) {
                edge = WhatEdge(triangles[j].p3, triangles[j].p1, triangles[j].ta_e3);
                if (edge == 1)
                    triangles[j].ta_e3->ta_e1 = NULL;
                else if (edge == 2)
                    triangles[j].ta_e3->ta_e2 = NULL;
                else if (edge == 3)
                    triangles[j].ta_e3->ta_e3 = NULL;
            }
            
            triangles[j].p1[0] = -100000000.0f;             //Mark triangle for deletion
        }
    }
    for (j = ncells - 1; j >= 0; j--) { 
        if (triangles[j].p1[0] == -100000000.0f) {
            // preserves pointer locations of edge triangles for tri with index < j
            for (int v = 0; v < triangles.size(); v++) {
                 if (triangles[v].ta_e1 > &(triangles[j]))
                     triangles[v].ta_e1 = triangles[v].ta_e1 - 1;
                 if (triangles[v].ta_e2 > &(triangles[j]))
                     triangles[v].ta_e2 = triangles[v].ta_e2 - 1;
                 if (triangles[v].ta_e3 > &(triangles[j]))
                     triangles[v].ta_e3 = triangles[v].ta_e3 - 1; 
             }
            triangles.erase(triangles.begin() + j);
        }
    }
}

void DelaunayTriangulation::CopyPoint(float * to, float * from) {
    to[0] = from[0];
    to[1] = from[1];
}

void DelaunayTriangulation::EdgeFlip(int j, float* p4_og, int edge)
{
    /*
     Find the points that share an edge with the 4th point inside the circumcircle. Get this info by which if statement above returns 'true'.
     These points should no longer have an edge between them. Therefore, flip that edge to be between the 4th point and the other point in the triangle.
     
			      1								 1
			    / | \						       /   \
			  /   |   \						     /       \
			2     |     4    <-- Does not meet DT requirement          2 _________ 4	
			  \   |   /			Flip to instead    --->      \       /
			    \ | /						       \   /
			      3								 3
     Then update points, edges, etc. of affected triangles to keep DS up to date.
    */
    int new_edge;
    float *p4 = new float[2];
    p4[0] = p4_og[0];
    p4[1] = p4_og[1];
    OneTriangle *t, *ot;
    
    if (edge == 1) {
        edge = WhatEdge(triangles[j].p1, triangles[j].p2, triangles[j].ta_e1);
            new_edge = (edge % 3) + 1;
            if (new_edge == 1)
                t = triangles[j].ta_e1->ta_e1;
            else if (new_edge == 2)
                t = triangles[j].ta_e1->ta_e2;
            else if (new_edge == 3 )
                t = triangles[j].ta_e1->ta_e3;
    
            new_edge = ((edge + 1) % 3) + 1;
            if (new_edge == 1)
                ot = triangles[j].ta_e1->ta_e1;
            else if (new_edge == 2)
                ot = triangles[j].ta_e1->ta_e2;
            else if (new_edge == 3)
                ot = triangles[j].ta_e1->ta_e3;

        CopyPoint(triangles[j].ta_e1->p1, p4);
        CopyPoint(triangles[j].ta_e1->p2, triangles[j].p2);
        CopyPoint(triangles[j].ta_e1->p3, triangles[j].p3);

        //FIX triangle across
        triangles[j].ta_e1->ta_e2 = triangles[j].ta_e2;
        triangles[j].ta_e1->ta_e3 = &(triangles[j]);
        triangles[j].ta_e1->ta_e1 = ot;
       
        triangles[j].ta_e2 = triangles[j].ta_e1;
        triangles[j].ta_e1 = t;

        //triangles[j].p1 and triangles[j].p3 are the same
        triangles[j].p2[0] = p4[0];
        triangles[j].p2[1] = p4[1];

        //Fix t across
        edge = WhatEdge(triangles[j].p1, triangles[j].p2, triangles[j].ta_e1);
        if (edge == 1)
            triangles[j].ta_e1 -> ta_e1 = &(triangles[j]);
        else if (edge == 2)
            triangles[j].ta_e1 -> ta_e2 = &(triangles[j]);
        else if (edge == 3)
            triangles[j].ta_e1 -> ta_e3 = &(triangles[j]);
       
        //Fix x across
        edge = WhatEdge(triangles[j].ta_e2 -> p2, triangles[j].ta_e2 -> p3, triangles[j].ta_e2 -> ta_e2); 
        if (edge == 1)
            triangles[j].ta_e2 -> ta_e2 -> ta_e1 = triangles[j].ta_e2;
        else if (edge == 2)
            triangles[j].ta_e2 -> ta_e2 -> ta_e2 = triangles[j].ta_e2;
        else if (edge == 3)
            triangles[j].ta_e2 -> ta_e2 -> ta_e3 = triangles[j].ta_e2;
    } 
    
    else if (edge == 2) {
        edge = WhatEdge(triangles[j].p2, triangles[j].p3, triangles[j].ta_e2);
            new_edge = ( ( (edge % 3) + 1) % 3 ) + 1;
            if (new_edge == 1)
                t = triangles[j].ta_e2->ta_e1;
            else if (new_edge == 2)
                t = triangles[j].ta_e2->ta_e2;
            else if (new_edge == 3)
                t = triangles[j].ta_e2->ta_e3;

            new_edge = ((edge % 3) + 1);
            if (new_edge == 1)
                ot = triangles[j].ta_e2->ta_e1;
            else if (new_edge == 2)
                ot = triangles[j].ta_e2->ta_e2;
            else if (new_edge == 3)
                ot = triangles[j].ta_e2->ta_e3;

        triangles[j].ta_e2->p1[0] = triangles[j].p1[0];
        triangles[j].ta_e2->p1[1] = triangles[j].p1[1];
        triangles[j].ta_e2->p2[0] = triangles[j].p2[0];
        triangles[j].ta_e2->p2[1] = triangles[j].p2[1];
        triangles[j].ta_e2->p3[0] = p4[0];
        triangles[j].ta_e2->p3[1] = p4[1];
       
        //FIX triangle across
        triangles[j].ta_e2->ta_e1 = triangles[j].ta_e1;
        triangles[j].ta_e2->ta_e3 = &(triangles[j]);
        triangles[j].ta_e2->ta_e2 = ot;
       
        triangles[j].ta_e1 = triangles[j].ta_e2;
        triangles[j].ta_e2 = t; 

        //triangles[j].p1 and triangles[j].p3 are the same
        triangles[j].p2[0] = p4[0];
        triangles[j].p2[1] = p4[1];
       
        edge = WhatEdge(triangles[j].p2, triangles[j].p3, triangles[j].ta_e2);
        if (edge == 1)
            triangles[j].ta_e2 -> ta_e1 = &(triangles[j]);
        else if (edge == 2)
            triangles[j].ta_e2 -> ta_e2 = &(triangles[j]);
        else if (edge == 3)
            triangles[j].ta_e2 -> ta_e3 = &(triangles[j]);
       
        edge = WhatEdge(triangles[j].ta_e1 -> p1, triangles[j].ta_e1 -> p2, triangles[j].ta_e1 -> ta_e1); 
        if (edge == 1)
            triangles[j].ta_e1 -> ta_e1 -> ta_e1 = triangles[j].ta_e1;
        else if (edge == 2)
            triangles[j].ta_e1 -> ta_e1 -> ta_e2 = triangles[j].ta_e1;
        else if (edge == 3)
            triangles[j].ta_e1 -> ta_e1 -> ta_e3= triangles[j].ta_e1;
    } 
    else if (edge == 3) {
        edge = WhatEdge(triangles[j].p3, triangles[j].p1, triangles[j].ta_e3);
       
            new_edge = (edge % 3) + 1;
            if (new_edge == 1)
                t = triangles[j].ta_e3->ta_e1;
            else if (new_edge == 2)
                t = triangles[j].ta_e3->ta_e2;
            else if (new_edge == 3)
                t = triangles[j].ta_e3->ta_e3;
        
            new_edge = (((edge + 1) % 3) + 1);
            if (new_edge == 1)
                ot = triangles[j].ta_e3->ta_e1;
            else if (new_edge == 2)
                ot = triangles[j].ta_e3->ta_e2;
            else if (new_edge == 3)
                ot = triangles[j].ta_e3->ta_e3;

        triangles[j].ta_e3->p1[0] = triangles[j].p1[0];
        triangles[j].ta_e3->p1[1] = triangles[j].p1[1];
        triangles[j].ta_e3->p2[0] = triangles[j].p2[0];
        triangles[j].ta_e3->p2[1] = triangles[j].p2[1];
        triangles[j].ta_e3->p3[0] = p4[0];
        triangles[j].ta_e3->p3[1] = p4[1];
       
        //FIX triangle across
        triangles[j].ta_e3->ta_e1 = triangles[j].ta_e1;
        triangles[j].ta_e3->ta_e2 = &(triangles[j]);
        triangles[j].ta_e3->ta_e3 = ot;
       
        triangles[j].ta_e1 = triangles[j].ta_e2;
        triangles[j].ta_e2 = t;

        triangles[j].p1[0] = triangles[j].p2[0];
        triangles[j].p1[1] = triangles[j].p2[1];
        triangles[j].p2[0] = triangles[j].p3[0];
        triangles[j].p2[1] = triangles[j].p3[1];
        triangles[j].p3[0] = p4[0];
        triangles[j].p3[1] = p4[1];
      
        //fix x 
        edge = WhatEdge(triangles[j].ta_e3 -> p1, triangles[j].ta_e3 -> p2, triangles[j].ta_e3 -> ta_e1);
        if (edge == 1)
            triangles[j].ta_e3 -> ta_e1 -> ta_e1 = triangles[j].ta_e3;
        else if (edge == 2)
            triangles[j].ta_e3 -> ta_e1 -> ta_e2 = triangles[j].ta_e3;
        else if (edge == 3)
            triangles[j].ta_e3 -> ta_e1 -> ta_e3 = triangles[j].ta_e3;
       
        //fix t
        edge = WhatEdge(triangles[j].p2, triangles[j].p3, triangles[j].ta_e2); 
        if (edge == 1)
            triangles[j].ta_e2 -> ta_e1 = &(triangles[j]);
        else if (edge == 2)
            triangles[j].ta_e2 -> ta_e2 = &(triangles[j]);
        else if (edge == 3)
            triangles[j].ta_e2 -> ta_e3 = &(triangles[j]);
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
    
float *
DelaunayTriangulation::Initialize(float *bounding_box, int point_count)
{
    float *bounding_tri = new float[6];
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
    
    //Reserve a contiguous portion of memory for that triangles, note: could be source of error if large
    triangles.reserve(4 * point_count + 1);
    triangles.emplace_back(ot);
    return bounding_tri;
}

void
DelaunayTriangulation::AddPoint(float x1, float y1)
{
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

    for (int i = 0 ; i < triangles.size() ; i++) {
        if (triangles[i].ContainsPoint(x1, y1)) {
            OneTriangle original_triangle = triangles[i];
            OneTriangle *TA = original_triangle.ta_e1;
            OneTriangle *TC = original_triangle.ta_e2; 
            OneTriangle *TB = original_triangle.ta_e3;

            // split triangle i into three triangles
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
        
            T1->ta_e1 = TA; 
            T1->ta_e2 = T3;
            T1->ta_e3 = T2;
            T2->ta_e1 = T1;
            T2->ta_e2 = T3;
            T2->ta_e3 = TB;
            T3->ta_e1 = T1;
            T3->ta_e2 = TC; 
            T3->ta_e3 = T2;
            if (TA != NULL) {
                int edge = WhatEdge(original_triangle.p1, original_triangle.p2, TA);
                if (edge == 1) {
                    TA->ta_e1 = T1;
                }
                else if (edge == 2) {
                    TA->ta_e2 = T1;
                }
                else if (edge == 3) {
                    TA->ta_e3 = T1;
                }
            }
            if (TB != NULL) {
                int edge = WhatEdge(original_triangle.p1, original_triangle.p3, TB);
                if (edge == 1) {
                    TB->ta_e1 = T2;
                }
                else if (edge == 2) {
                    TB->ta_e2 = T2;
                }
                else if (edge == 3) {
                    TB->ta_e3 = T2;
                }
            }
            if (TC != NULL) {
                int edge = WhatEdge(original_triangle.p2, original_triangle.p3, TC);
                if (edge == 1) {
                    TC->ta_e1 = T3;
                }
                else if (edge == 2) {
                    TC->ta_e2 = T3;
                }
                else if (edge == 3) {
                    TC->ta_e3 = T3;
                }
            }
            return;
        }
    }
}

// Inputs: 3 points, each with x and y coordinates (or z in case of 3D) and a 4th point with x and y coordinates. This makes up points A, B, C (of the
// triangle) and the 4th point, D.
// Output: Boolean value. True if 4th point is inside circle. This means we have to flip. False if 4th point is outside circle. This means we're ok.
bool 
DelaunayTriangulation::CircumcircleCheck(float* ptA, float* ptB, float* ptC, float* ptD)
{
    float result = 0.0;

    //find the Determinant
    float Part1 = (ptB[1] - ptD[1]) * (DetHelp(ptC[0], ptD[0], ptC[1], ptD[1])) - (ptC[1] - ptD[1]) * (DetHelp(ptB[0], ptD[0], ptB[1], ptD[1]));
    float Part2 = (ptB[0] - ptD[0]) * (DetHelp(ptC[0], ptD[0], ptC[1], ptD[1])) - (ptC[0] - ptD[0]) * (DetHelp(ptB[0], ptD[0], ptB[1], ptD[1]));
    float Part3 = (ptB[0] - ptD[0]) * (ptC[1] - ptD[1]) - (ptC[0] - ptD[0]) * (ptB[1] - ptD[1]);

    float A1 = (ptA[0] - ptD[0]) * Part1;
    float A2 = (ptA[1] - ptD[1]) * Part2;
    float A3 = DetHelp(ptA[0], ptD[0], ptA[1], ptD[1]) * Part3;
    result = A1 - A2 + A3;

    if (result > 0)
        return true;                   // ptD lies inside circle
    else
        return false;                  // ptD lies outside or on the circle
}

bool
DelaunayTriangulation::AltCircumcircleCheck(float *ptA, float *ptB, float *ptC, float *ptD)
{
    float ax_ = ptA[0]-ptD[0];
    float ay_ = ptA[1]-ptD[1];
    float bx_ = ptB[0]-ptD[0];
    float by_ = ptB[1]-ptD[1];
    float cx_ = ptC[0]-ptD[0];
    float cy_ = ptC[1]-ptD[1];

    //This calculates radius of circumcircle defined by A, B, C
    float a = sqrt( pow(ptA[0] - ptB[0], 2) + pow(ptA[1] - ptB[1], 2) );
    float b = sqrt( pow(ptB[0] - ptC[0], 2) + pow(ptB[1] - ptC[1], 2) );
    float c = sqrt( pow(ptC[0] - ptA[0], 2) + pow(ptC[1] - ptA[1], 2) );
    
    float radius = a * b * c;
    radius /= sqrt( (a+b+c) * (b+c-a) * (c+a-b) * (a+b-c));

    float result = ((ax_ * ax_ + ay_ * ay_) * (bx_ * cy_ - cx_ * by_) -
                    (bx_ * bx_ + by_ * by_) * (ax_ * cy_ - cx_ * ay_) +
                    (cx_ * cx_ + cy_ * cy_) * (ax_ * by_ - bx_ * ay_));
    
    return result > 0;
}

//Returns the side of the triangle that is composed of these two points 
int
DelaunayTriangulation::WhatEdge(float *pt1, float *pt2, OneTriangle *tri)
{
    int total = 0;
    float EPSILON = 0.000001f;

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
    else {                   //Not in Triangle, Shouldn't return this, function should not be called on non adjacent triangles
        printf("Triangle didn't have the point\n");
        printf("%f\t%f\n", pt1[0], pt1[1]);
        printf("%f\t%f\n", pt2[0], pt2[1]);
        PrintTri(tri);
        return 0;
    }
}

//Helper function for CircumcirlceCheck function. This helps the readability of the code. Does a simple calculation on 4 values, returns result.
float
DelaunayTriangulation::DetHelp(float pt1, float pt2, float pt3, float pt4)
{
    return ((pow((pt1 - pt2), 2)) + (pow((pt3 - pt4), 2)));
}

//Function for debugging. Prints info about triangle pointed to by *t
void
DelaunayTriangulation::PrintTri(OneTriangle *t)
{
    printf("Triangle at addr:\t%p\n", t );
    printf("P1: %f\t%f\n", t->p1[0], t->p1[1]);
    printf("P2: %f\t%f\n", t->p2[0], t->p2[1]);
    printf("P3: %f\t%f\n", t->p3[0], t->p3[1]);
    printf("Triangle across e1:\t%p\n", t->ta_e1);
    printf("Triangle across e2:\t%p\n", t->ta_e2);
    printf("Triangle across e3:\t%p\n", t->ta_e3);
    printf("**********************************************\n\n");
}

float *
DelaunayTriangulation::FindBoundingBox(float *points) 
{
    float *bounding_box = new float[4];
    float x_min = 0.0f;
    float x_max = 0.0f;
    float y_min = 0.0f;
    float y_max = 0.0f;

    for (int i = 0; i < NUM_POINTS; i++) {
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

//        creates vectors (outsidevectors) for use in sumangles in 2nd verify function
float *
DelaunayTriangulation::FindVectors(int edgeOfEdgeTri, OneTriangle * overEdge) 
{
    // will have 2 vectors, vector 1 is vectors[0-1], vector 2 is vectors[2-3]
    float * vectors = new float[4];
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
bool   DelaunayTriangulation::SumAngles(float * insideV, float * outsideV) 
{ 
    float m1 = sqrt((insideV[0] * insideV[0]) + (insideV[1] * insideV[1]));
    float m2 = sqrt((insideV[2] * insideV[2]) + (insideV[3] * insideV[3]));
    float insideRadians = acos((insideV[0]*insideV[2] + insideV[1] * insideV[3]) / (m1 * m2) );
    m1 = sqrt((outsideV[0] * outsideV[0]) + (outsideV[1] * outsideV[1]));
    m2 = sqrt((outsideV[2] * outsideV[2]) + (outsideV[3] * outsideV[3]));
    float outsideRadians = acos((outsideV[0]*outsideV[2] + outsideV[1] * outsideV[3]) / (m1 * m2) );
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
    float * outsideVectors;
    float * insideVectors = new float[4];
    for (int i = 0; i < numTriangles; i++) {
        if (triangles[i].p1[0] != -100000000.0f) {
            if (triangles[i].ta_e1) {    // check over edge 1
                insideVectors[0] = triangles[i].p1[0] - triangles[i].p3[0];
                insideVectors[1] = triangles[i].p1[1] - triangles[i].p3[1];
                insideVectors[2] = triangles[i].p2[0] - triangles[i].p3[0];
                insideVectors[3] = triangles[i].p2[1] - triangles[i].p3[1];
                int outEdge = WhatEdge(triangles[i].p1, triangles[i].p2, triangles[i].ta_e1);
                outsideVectors = FindVectors(outEdge, triangles[i].ta_e1);
                if (SumAngles(insideVectors, outsideVectors)) {
                    cerr << "SUM > 180 DEGREES OVER EDGE 1 FOR TRI INDEX: " << i << "\n------------" << endl;
                    PrintTri(&triangles[i]);
                    PrintTri(triangles[i].ta_e1);
                }
            }
            if (triangles[i].ta_e2) {    // check over edge 2
                insideVectors[0] = triangles[i].p2[0] - triangles[i].p1[0];
                insideVectors[1] = triangles[i].p2[1] - triangles[i].p1[1];
                insideVectors[2] = triangles[i].p3[0] - triangles[i].p1[0];
                insideVectors[3] = triangles[i].p3[1] - triangles[i].p1[1];
                int outEdge = WhatEdge(triangles[i].p2, triangles[i].p3, triangles[i].ta_e2);
                outsideVectors = FindVectors(outEdge, triangles[i].ta_e2);
                if (SumAngles(insideVectors, outsideVectors)) {
                    cerr << "SUM > 180 DEGREES OVER EDGE 2 FOR TRI INDEX: " << i << "\n------------" << endl;
                    PrintTri(&triangles[i]);
                    PrintTri(triangles[i].ta_e2);
                }
            }
            if (triangles[i].ta_e3) {   // check over edge 3
                insideVectors[0] = triangles[i].p1[0] - triangles[i].p2[0];
                insideVectors[1] = triangles[i].p1[1] - triangles[i].p2[1];
                insideVectors[2] = triangles[i].p3[0] - triangles[i].p2[0];
                insideVectors[3] = triangles[i].p3[1] - triangles[i].p2[1];
                int outEdge = WhatEdge(triangles[i].p3, triangles[i].p1, triangles[i].ta_e3);
                outsideVectors = FindVectors(outEdge, triangles[i].ta_e3);
                if (SumAngles(insideVectors, outsideVectors)) {
                    cerr << "SUM > 180 DEGREES OVER EDGE 3 FOR TRI INDEX: " << i << "\n-----------" << endl;
                    PrintTri(&triangles[i]);
                    PrintTri(triangles[i].ta_e3);
                }
            }
        }
    }
    delete [] outsideVectors;
    delete [] insideVectors;
}

int main()
{
    float *pts = PointsGenerator(NUM_POINTS, 2);
    DelaunayTriangulation DT;
    
    //Declare timers
    struct timeval start, end;
    
    //Determine bounding box
    float *bounding_box = DT.FindBoundingBox(pts);

    //Make initial triangulation.  Allocate memory for vector
    float *bounding_tri = DT.Initialize(bounding_box, NUM_POINTS);
    
    //AddPoints to triangulation.  Produces and initial tesselation
    gettimeofday(&start, NULL);
    for (int i = 0 ; i < NUM_POINTS ; i++) {
        DT.AddPoint(pts[2*i], pts[2*i+1]);
    }
    gettimeofday(&end, NULL);
    
    //AddPoint time report
    double runtime = end.tv_sec + (end.tv_usec / 1000000.0) - start.tv_sec - (start.tv_usec / 1000000.0);
    printf("\nComputation time for %d AddPoint Operations: %.4f s\n\n", NUM_POINTS, runtime);

    //Make tessellation meet Delaunay condition
    gettimeofday(&start, NULL);
    DT.Verify(); 
    gettimeofday(&end, NULL);

    //Verify time report
    runtime = end.tv_sec + (end.tv_usec / 1000000.0) - start.tv_sec - (start.tv_usec / 1000000.0);
    printf("\nComputation time for Verify method: %.4f s\n\n", runtime);

    //Delete triangles connected to bounding triangle vertices
    gettimeofday(&start, NULL);
    DT.DelBoundingTri(bounding_tri);
    gettimeofday(&end, NULL);

    // check for delaunay condition being met using sumangles method  
    DT.VerifyMeetDC();
 
    //DelBoundingTriangle time report
    runtime = end.tv_sec + (end.tv_usec / 1000000.0) - start.tv_sec - (start.tv_usec / 1000000.0);
    printf("\nComputation time for DelBoundingTriangle method: %.4f s\n\n", runtime);
    
    char *filename = (char *)"kristi.vtk";
    DT.WriteOutTriangle(filename);

    delete [] pts;
    delete [] bounding_box;
    delete [] bounding_tri;
    return 0;
}
