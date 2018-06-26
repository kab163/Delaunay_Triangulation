#include <iostream>
#include <cmath>

#include "visit_writer.c"
#include <vector>
#include <logging.h>

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
            float rand_value = rand() % 100000 / 100000.0;
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
    OneTriangle  *triangle_across_e1;
    OneTriangle  *triangle_across_e2;
    OneTriangle  *triangle_across_e3;

    bool      ContainsPoint(float x, float y);
    OneTriangle()
    {
        triangle_across_e1 = NULL;
        triangle_across_e2 = NULL;
        triangle_across_e3 = NULL;
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
    void   Initialize(float, float, float, float, float, float);
    void   AddPoint(float, float);
    bool   CircumcircleCheck(float*, float*, float*, float*);
    void   Verify();
    void   DelBoundingTri();
    void   WriteOutTriangle(char *filename);
    int    WhatEdge(float *, float *, OneTriangle *);
    //int numA = 0, numB = 0, numC = 0, numD = 0, numE = 0, numF = 0;
  private:
    std::vector<OneTriangle>  triangles;
    float DetHelp(float, float, float, float);
    float * determineFourthPoint(OneTriangle, OneTriangle *, int *);
    void EdgeFlip(int, float*, OneTriangle * edgeTriangle, int edge, int *);
};

/*
    Function to determine fourth point to used in edge flip
    p1 to p1 - p2 to p2 -> p3
    p1 to p2 - p2 to p1 -> p3
    p1 to p3 - p2 to p1 -> p2
    p2 to p3 - p3 to p2 -> p1
     
     
			  opx   bpx							 
			    / |  \
			  /   |   \
			opy   |    bpy  
			  \   |   /
			    \ |  /
			  opz   bpz
    works by comparing bpx and bpz to opx and opz (respectively)
        to determine the point of opy (the fourth point for determining 
        whether in circumcircle)
    NOTE: NEED TO MODIFY TO RETURN WHAT POINT NUMBER IT IS AS WELL - WILL BE USED IN EDGEFLIP
*/
float *
DelaunayTriangulation::determineFourthPoint(OneTriangle base, OneTriangle * overEdge, int * pLoc) {
    if ((base.p1[0] == overEdge->p1[0]) && (base.p1[1] == overEdge->p1[1])) {        // base p1 matches overedge p1
        if ((base.p2[0] == overEdge->p2[0]) && (base.p2[1] == overEdge->p2[1])) {      // base p2 matches overedge p2
            *pLoc = 3;
             //numA++;
             return overEdge->p3;
        }
        else {
            *pLoc = 2;
            //numB++;
            return overEdge->p2;
        }
    }
    else if ((base.p1[0] == overEdge->p2[0]) && (base.p1[1] == overEdge->p2[1])) {   // base p1 matches overedge p2
        if ((base.p2[0] == overEdge->p1[0]) && (base.p2[1] == overEdge->p1[1])) {
            *pLoc = 3;
            //numC++;
            return overEdge->p3;
        }
        else {
            *pLoc = 1;
            //numD++;
            return overEdge->p1;
        }
    }
    else if ((base.p1[0] == overEdge->p3[0]) && (base.p1[1] == overEdge->p3[1])) {
        if ((base.p2[0] == overEdge->p2[0]) && (base.p2[1] == overEdge->p2[1])) {
            *pLoc = 3;
            //numE++;
            return overEdge->p3;
        }
        else {
            *pLoc = 2;
            //numF++;
            return overEdge->p2;
        }
    }
    
}

void DelaunayTriangulation::Verify()
{
    int ncells = triangles.size();
    //// ADDED
    char msg[1024];
    //
    int iteration = 0;
    int totalFlips = 0;
    int numTrianglesFlipped = 0;
    int pointLocation = 0;
    float *fourthPoint;
    bool done = false;
    
    while (!done) {
      numTrianglesFlipped = 0;
      for(int j = 0; j < ncells; j++) {   
        if (triangles[j].triangle_across_e1 != NULL) {
            //// ADDED
            // logging
            sprintf(msg, "verify:\ttriangle %d - TAe1 not null", j);
            Logger::LogEvent(msg);
            // finding 4th point
            fourthPoint = determineFourthPoint(triangles[j], triangles[j].triangle_across_e1, &pointLocation);
            //if(CircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e1->p2)) {
            if(CircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, fourthPoint)) {
	        numTrianglesFlipped++; 
	        //EdgeFlip(j,triangles[j].triangle_across_e1->p2, 1);
	        EdgeFlip(j, fourthPoint, triangles[j].triangle_across_e1, 1, &pointLocation);
                //// END ADDED
	    }
        }
        if (triangles[j].triangle_across_e2 != NULL) {
            //// ADDED
            // logging
            sprintf(msg, "verify:\ttriangle %d - TAe2 not null", j);
            Logger::LogEvent(msg);
            // finding 4th point
            fourthPoint = determineFourthPoint(triangles[j], triangles[j].triangle_across_e2, &pointLocation);
            //if(CircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e2->p2)) { 
            if(CircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, fourthPoint)) { 
	        numTrianglesFlipped++;
		//EdgeFlip(j,triangles[j].triangle_across_e2->p2, 2);
		EdgeFlip(j, fourthPoint, triangles[j].triangle_across_e2, 2, &pointLocation);
                //// END ADDED
            }
        } 
        if (triangles[j].triangle_across_e3 != NULL) {
            //// ADDED
            // logging
            sprintf(msg, "verify:\ttriangle %d - TAe3 not null", j);
            Logger::LogEvent(msg);
            //
            // finding 4th point
            fourthPoint = determineFourthPoint(triangles[j], triangles[j].triangle_across_e3, &pointLocation);
            //if(CircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e3->p3)) { 
            if(CircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, fourthPoint)) { 
	        numTrianglesFlipped++;
	        //EdgeFlip(j, triangles[j].triangle_across_e3->p3, 3);
	        EdgeFlip(j, fourthPoint, triangles[j].triangle_across_e3, 3, &pointLocation);
                //// END ADDED
	    }
        }
      }
      totalFlips += numTrianglesFlipped;
      done = (numTrianglesFlipped == 0 ? true : false);
      iteration++;
    }
    Logger::LogEvent("size of vector: ", triangles.size());
    //printf("A: %d\n B: %d\n C: %d\n D: %d\n E: %d\n F: %d\n", numA, numB, numC, numD, numE, numF);
    printf("Iteration count: %d\n", iteration);
    printf("Total flips: %d\n", totalFlips);
}

void DelaunayTriangulation::DelBoundingTri() 
{
    /*
      TODO:
      Here is where I should delete the first, bounding triangle - update any triangles who have a triangle_across_e*
      that is this bounding triangle. The DT should now be complete.
    */
}

/*
//// ADDED 
// function for point comparisons
bool DelaunayTriangulation::PointCompare(float ax, float ay, float bx, float by) {
    if ((ax = bx) && (ay == by))
        return true;
    else
        return false;
}
//// END ADDED
*/
void DelaunayTriangulation::EdgeFlip(int j, float* p4, OneTriangle * edgeTriangle, int edge, int * pLoc)
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

    /*
    if (PointCompare(edgeTriangle->p1[0], triangles[j].p1[0], edgeTriangle->p1[1], triangles[j].p1[1])) {
        if (PointCompare(edgeTriangle->p2[0], triangles[j].p3[0], edgeTriangle->p2[1], triangles[j].p3[1])) {
            // then fourth point is edgetriangle p3
            
        edgeTriangle->p1[0] = p4[0];
        edgeTriangle->p1[1] = p4[1];
        edgeTriangle->p2[0] = triangles[j].p2[0];
        edgeTriangle->p2[1] = triangles[j].p2[1];
        edgeTriangle->p3[0] = triangles[j].p3[0];
        edgeTriangle->p3[1] = triangles[j].p3[1];
        triangles[j].p3[0] = p4[0];
        triangles[j].p3[1] = p4[1];
    } else {
        edgeTriangle->p1[0] = triangles[j].p1[0];
        edgeTriangle->p1[1] = triangles[j].p1[1];
        edgeTriangle->p2[0] = triangles[j].p3[0];
        edgeTriangle->p2[1] = triangles[j].p3[1];
        triangles[j].p3[0] = edgeTriangle->p3[0];
        triangles[j].p3[1] = edgeTriangle->p3[1];
        edgeTriangle->p3[0] = p4[0];
        edgeTriangle->p3[1] = p4[1];
        
    }

    */
    float * bP1 = triangles[j].p1;
    float * bP2 = triangles[j].p2;
    float * bP3 = triangles[j].p3;
    float * eP1 = edgeTriangle->p1;
    float * eP2 = edgeTriangle->p2;
    float * eP3 = edgeTriangle-> p3;
    if (edge == 1) {
        if (*pLoc == 1) {
            triangles[j].p2 = p4;
            int edgeTedge = WhatEdge(bP1, bP2, edgeTriangle);
            if (eP2 == bP1) {		// then next point is clockwise if orientation not mirred as in addpoint convention
        }
    }
            
    /*
    if (edge == 1) {
       triangles[j].triangle_across_e1->p1[0] = p4[0];
       triangles[j].triangle_across_e1->p1[1] = p4[1];
       triangles[j].triangle_across_e1->p2[0] = triangles[j].p2[0];
       triangles[j].triangle_across_e1->p2[1] = triangles[j].p2[1];
       triangles[j].triangle_across_e1->p3[0] = triangles[j].p3[0];
       triangles[j].triangle_across_e1->p3[1] = triangles[j].p3[1];
       //triangles[j].p1 and triangles[j].p3 are the same
       triangles[j].p2[0] = p4[0];
       triangles[j].p2[1] = p4[1];
       //// ADDED:
       OneTriangle * original_e2 = triangles[j].triangle_across_e2;
       triangles[j].triangle_across_e2 = triangles[j].triangle_across_e3;
       triangles[j].triangle_across_e3 = edgeTriangle->triangle_across_e1;
       edgeTriangle->triangle_across_e1 = edgeTriangle->triangle_across_e2;
       edgeTriangle->triangle_across_e2 = original_e2;
       ////
    } else if (edge == 2) {
       triangles[j].triangle_across_e2->p1[0] = triangles[j].p1[0];
       triangles[j].triangle_across_e2->p1[1] = triangles[j].p1[1];
       triangles[j].triangle_across_e2->p2[0] = triangles[j].p2[0];
       triangles[j].triangle_across_e2->p2[1] = triangles[j].p2[1];
       triangles[j].triangle_across_e2->p3[0] = p4[0];
       triangles[j].triangle_across_e2->p3[1] = p4[1];
       //triangles[j].p1 and triangles[j].p3 are the same
       triangles[j].p2[0] = p4[0];
       triangles[j].p2[1] = p4[1];
       //// ADDED:
       OneTriangle * original_e1 = triangles[j].triangle_across_e1;
       triangles[j].triangle_across_e1 = triangles[j].triangle_across_e2;
       triangles[j].triangle_across_e2 = edgeTriangle->triangle_across_e2;
       edgeTriangle->triangle_across_e2 = edgeTriangle->triangle_across_e1;
       edgeTriangle->triangle_across_e1 = original_e1;
       ////
    } else if (edge == 3) {
       triangles[j].triangle_across_e3->p1[0] = triangles[j].p1[0];
       triangles[j].triangle_across_e3->p1[1] = triangles[j].p1[1];
       triangles[j].triangle_across_e3->p2[0] = triangles[j].p2[0];
       triangles[j].triangle_across_e3->p2[1] = triangles[j].p2[1];
       triangles[j].triangle_across_e3->p3[0] = p4[0];
       triangles[j].triangle_across_e3->p3[1] = p4[1];

       triangles[j].p1[0] = triangles[j].p2[0];
       triangles[j].p1[1] = triangles[j].p2[1];
       triangles[j].p2[0] = triangles[j].p3[0];
       triangles[j].p2[1] = triangles[j].p3[1];
       triangles[j].p3[0] = p4[0];
       triangles[j].p3[1] = p4[1];
       //triangles[j].p1 and triangles[j].p3 are the same
       //// ADDED:
       OneTriangle * orig_e2 = edgeTriangle->triangle_across_e2;
       edgeTriangle->triangle_across_e1 = triangles[j].triangle_across_e1;
       triangles[j].triangle_across_e1 = triangles[j].triangle_across_e2;
       edgeTriangle->triangle_across_e2 = &(triangles[j]);
       triangles[j].triangle_across_e2 = orig_e2;
       ////

    } else printf("\n\n\n***edge error!***\n\n\n");
    */
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

    // added removal of dynamically allocated memory
    if (celltypes)
        delete[] celltypes;
    if (pts)
        delete[] pts;
    if (conn)
        delete[] conn;
    //

}
    
void
DelaunayTriangulation::Initialize(float x1, float y1, float x2, float y2, float x3, float y3)
{
    OneTriangle ot;
    ot.p1[0] = x1;
    ot.p1[1] = y1;
    ot.p2[0] = x2;
    ot.p2[1] = y2;
    ot.p3[0] = x3;
    ot.p3[1] = y3;
    triangles.reserve(2000);
    triangles.push_back(ot);
}

void
DelaunayTriangulation::AddPoint(float x1, float y1)
{
    for (int i = 0 ; i < triangles.size() ; i++)
    {
        if (triangles[i].ContainsPoint(x1, y1))
        {
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
            //// ADDED
            /*
            if (triangles[i].triangle_across_e1 != original_triangle.triangle_across_e1)
                cerr << "ERROR! pointers not equal!" << endl; 
	    cerr << "original p1: " << triangles[i].p1[0] << ", " << triangles[i].p1[1] << endl
                 << "p2: " << triangles[i].p2[0] << ", " << triangles[i].p2[1] << endl
                 << "p3: " << triangles[i].p3[0] << ", " << triangles[i].p3[1] << endl;
	    cerr << "copy of original p1: " << original_triangle.p1[0] << ", " << original_triangle.p1[1] << endl
                 << "p2: " << original_triangle.p2[0] << ", " << original_triangle.p2[1] << endl
                 << "p3: " << original_triangle.p3[0] << ", " << original_triangle.p3[1] << endl;
            //// END ADDED
            */
            OneTriangle *TA = original_triangle.triangle_across_e1;
	    OneTriangle *TC = original_triangle.triangle_across_e2; //KB
	    OneTriangle *TB = original_triangle.triangle_across_e3; //KB
            // split triangle i into three triangles
            // note: no edge flipping or Delaunay business.
            // start by replacing triangle in the current list.
            triangles[i].p3[0] = x1;
            triangles[i].p3[1] = y1;

            // now add two more triangles
            OneTriangle new_triangle1;
            new_triangle1.p1[0] = original_triangle.p1[0]; //KB
            new_triangle1.p1[1] = original_triangle.p1[1]; //KB
            new_triangle1.p2[0] = x1; //KB
            new_triangle1.p2[1] = y1; //KB
            new_triangle1.p3[0] = original_triangle.p3[0]; //KB
            new_triangle1.p3[1] = original_triangle.p3[1]; //KB

            OneTriangle new_triangle2;
            new_triangle2.p1[0] = x1; //KB
            new_triangle2.p1[1] = y1; //KB
            new_triangle2.p2[0] = original_triangle.p2[0]; //KB
            new_triangle2.p2[1] = original_triangle.p2[1]; //KB
            new_triangle2.p3[0] = original_triangle.p3[0]; //KB
            new_triangle2.p3[1] = original_triangle.p3[1]; //KB

            // add to vector and triangles used in convention (used for comparisons)
            triangles.push_back(new_triangle1);
            int index = triangles.size() - 1; 
            OneTriangle *T2 = &(triangles[index]);
            triangles.push_back(new_triangle2);
            OneTriangle *T3 = &(triangles[index + 1]);
            OneTriangle *T1 = &(triangles[i]);
            
            T1->triangle_across_e2 = T3;
            T1->triangle_across_e3 = T2;
            T2->triangle_across_e1 = T1;
            T2->triangle_across_e2 = T3;
            T3->triangle_across_e1 = T1;
            T3->triangle_across_e3 = T2;
            /*
            //// ADDED:
            //
            cerr << "--T1 p1: " << T1->p1[0] << ", " << T1->p1[1] << endl
                 << "p2: " << T1->p2[0] << ", " << T1->p2[1] << endl
                 << "p3: " << T1->p3[0] << ", " << T1->p3[1] << endl;
	    cerr << "--T2 p1: " << T2->p1[0] << ", " << T2->p1[1] << endl
                 << "p2: " << T2->p2[0] << ", " << T2->p2[1] << endl
                 << "p3: " << T2->p3[0] << ", " << T2->p3[1] << endl;
	    cerr << "--T3 p1: " << T3->p1[0] << ", " << T3->p1[1] << endl
                 << "p2: " << T3->p2[0] << ", " << T3->p2[1] << endl
                 << "p3: " << T3->p3[0] << ", " << T3->p3[1] << endl
                 << "-----------------------------------" << endl;
            //// END ADDED
            */
          if (TA == NULL) {
              T1->triangle_across_e1 = NULL;
          } 
          else { 
              T1->triangle_across_e1 = TA;
              int edge = WhatEdge(original_triangle.p1, original_triangle.p2, TA);
              if (edge == 1)
                  TA->triangle_across_e1 = T1;
              else if (edge == 2)
                  TA->triangle_across_e2 = T1;
              else if (edge == 3)
                  TA->triangle_across_e3 = T1;
          } 
          if (TB == NULL) {
              T2->triangle_across_e3 = NULL;
          }
          else {
              T2->triangle_across_e3 = TB;
              int edge = WhatEdge(original_triangle.p1, original_triangle.p3, TB);
              if (edge == 1)
                  TB->triangle_across_e1 = T2;
              else if (edge == 2)
                  TB->triangle_across_e2 = T2;
              else if (edge == 3)
                  TB->triangle_across_e3 = T2;
          }
          if (TC == NULL) {
              T3->triangle_across_e2 = NULL;
          }
          else{
              T3->triangle_across_e2 = TC;
              int edge = WhatEdge(original_triangle.p2, original_triangle.p3, TC);
              if (edge == 1)
                  TC->triangle_across_e1 = T3;
              else if (edge == 2)
                  TC->triangle_across_e2 = T3;
              else if (edge == 3)
                  TC->triangle_across_e3 = T3;
          }
          break;
        }
    }
}

// Function to calculate the determinant of 3 points. This gives us the circumcircle of the triangle. If a 4th point is inside the triangle,
// the result will be negative. If the 4th point lies outside the circle, result will be positive. A result equal to zero means that the 4th
// point lies on the circle exactly. In this case, the DT of the set of points is not unique. It's like drawing two triangles in a square. 
// Whether the 3rd point and the 0th point make up the hypotenous or the 2nd and the 1st point make up the hypotenus, it's equivalent and
// therefore you could do either one and make a valid DT. Will need to call a seperate function to handle that case later. For now I return false,
// meaning that I do not flip it and leave it as is...
//
// Inputs: 3 points, each with x and y coordinates (or z in case of 3D) and a 4th point with x and y coordinates. This makes up points A, B, C (of the
// triangle) and the 4th point, D.
// Output: Boolean value. True if 4th point is inside circle. This means we have to flip. False if 4th point is outside circle. This means we're ok.
//
// TODO: create another case where we call a function to handle if the result is equal to zero - for now, return false  
bool 
DelaunayTriangulation::CircumcircleCheck(float* ptA, float* ptB, float* ptC, float* ptD)
{
    float result = 0.0;

    //find the Determinant
    ///*
    float Part1 = (ptB[1] - ptD[1]) * (DetHelp(ptC[0], ptD[0], ptC[1], ptD[1])) - (ptC[1] - ptD[1]) * (DetHelp(ptB[0], ptD[0], ptB[1], ptD[1]));
    float Part2 = (ptB[0] - ptD[0]) * (DetHelp(ptC[0], ptD[0], ptC[1], ptD[1])) - (ptC[0] - ptD[0]) * (DetHelp(ptB[0], ptD[0], ptB[1], ptD[1]));
    float Part3 = (ptB[0] - ptD[0]) * (ptC[1] - ptD[1]) - (ptC[0] - ptD[0]) * (ptB[1] - ptD[1]);

    float A1 = (ptA[0] - ptD[0]) * Part1;
    float A2 = (ptA[1] - ptD[1]) * Part2;
    float A3 = DetHelp(ptA[0], ptD[0], ptA[1], ptD[1]) * Part3;
    
    result = A1 - A2 + A3;
    /*
    float firstX = ptA[0] - ptD[0];
    float firstY = ptA[1] - ptD[1];
    float secondX = ptB[0] - ptD[0];
    float secondY = ptB[1] - ptD[1];
    float thirdX = ptC[0] - ptD[0];
    float thirdY = ptC[1] - ptD[1];
    result = ((firstX*firstX) + (firstY*firstY))*((secondX*thirdY)-(thirdX*secondY))-((secondX*secondX)+(secondY*secondY))*((firstX*thirdY)-(thirdX*firstY))+((thirdX*thirdX)+(thirdY*thirdY))*((firstX*secondY)-(secondX*firstY));
    */

    //cerr << "RESULT IS: " << result << endl;
    if (result < 0.0) return false; //ptD lies outside circumcircle
    else if (result > 0.0) return true; //ptD lies inside circumcircle
    else if (result == 0.0) return false; // ptD lies ON the circle, for now acting like it lies outside...
    //what about if it equals zero? Lies *on* the circle...
}

//Helper function for CircumcirlceCheck function. This helps the readability of the code. Does a simple calculation on 4 values, returns result.
float
DelaunayTriangulation::DetHelp(float pt1, float pt2, float pt3, float pt4)
{
    return ((pow((pt1 - pt2), 2)) + (pow((pt3 - pt4), 2)));
}


//Returns the side of the triangle that is composed of these two points
int
DelaunayTriangulation::WhatEdge(float *pt1, float *pt2, OneTriangle *tri)
{
    int total = 0;
    if (tri == NULL)
        return 0;
    if ((pt1[0] == tri->p1[0] && pt1[1] == tri->p1[1]) || (pt2[0] == tri->p1[0] && pt2[1] == tri->p1[1])) // bp1 matches op1 or bp2 matches op1
        total += 1;
    if ((pt1[0] == tri->p2[0] && pt1[1] == tri->p2[1]) || (pt2[0] == tri->p2[0] && pt2[1] == tri->p2[1])) // bp1 matches op2 or bp2 matches op2 
        total += 2;
    if ((pt1[0] == tri->p3[0] && pt1[1] == tri->p3[1]) || (pt2[0] == tri->p3[0] && pt2[1] == tri->p3[1]))
        total += 3;
    if (total == 3)         //Points 1 and 2
        return 1;
    else if (total == 4)    //Points 1 and 3
        return 3;
    else if (total == 5)    //Points 2 and 3
        return 2;
    else                    //Not in Triangle, Shouldn't return this
        return 0;
}

int main()
{
    int points = 100;
    float *pts = PointsGenerator(points, 2);
    DelaunayTriangulation DT;

    DT.Initialize(-1, -1, 2, -1, .5, 2);
    for (int i = 0 ; i < points ; i++)
        DT.AddPoint(pts[2*i], pts[2*i+1]);
 
    DT.Verify(); 
    //DT.DelBoundingTri();
    
    char *filename = (char *)"kristi.vtk";
    DT.WriteOutTriangle(filename);
    if (pts)
        delete [] pts; 
    return 0;
}
