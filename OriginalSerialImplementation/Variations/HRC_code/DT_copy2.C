
#include <iostream>
#include <cmath>

#include "visit_writer.c"
#include <vector>

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
    int numTrianglesFlipped;
    bool done = false;

    while (!done) {
      for(int j = 1; j < ncells; j++) {
	numTrianglesFlipped = 0;
        if (triangles[j].triangle_across_e1 != NULL) {
            if(CircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e1->p2)) {
	        numTrianglesFlipped++; 
	        EdgeFlip(j,triangles[j].triangle_across_e1->p2, 1);
	    }
        }
        if (triangles[j].triangle_across_e2 != NULL) {
            if(CircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e2->p3)) { 
	        numTrianglesFlipped++;
		EdgeFlip(j,triangles[j].triangle_across_e2->p3, 2);
            }
        } 
        if (triangles[j].triangle_across_e3 != NULL) {
            if(CircumcircleCheck(triangles[j].p1, triangles[j].p2, triangles[j].p3, triangles[j].triangle_across_e3->p3)) { 
	        numTrianglesFlipped++;
	        EdgeFlip(j, triangles[j].triangle_across_e3->p3, 3);
	    }
        }
        totalFlips += numTrianglesFlipped;
      }
      done = (numTrianglesFlipped == 0 ? true : false);
      iteration++;
    }

    printf("Iteration count: %d\n", iteration);
    printf("Total flips: %d\n", totalFlips);
}

void DelaunayTriangulation::DelBoundingTri() 
{
    /*
      Here is where I should delete the first, bounding triangle - update any triangles who have a triangle_across_e*
      that is this bounding triangle. the DT should now be complete.
    */
}

void DelaunayTriangulation::EdgeFlip(int j, float* p4, int edge)
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

    //TODO create function to do the edge flips

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

    } else printf("\n\n\n***edge error!***\n\n\n");/* */
}

void DelaunayTriangulation::WriteOutTriangle(char *filename)
{
    int ncells = triangles.size();
cerr << "NUMBER OF TRIANGLE is " << ncells << endl;
    int *celltypes = new int[ncells];
    for (int i = 0 ; i < ncells ; i++)
        celltypes[i] = VISIT_TRIANGLE;

    int dimensions = 3; // always 3 for VTK
    int vertices_per_cell = 3;
    int npts = ncells*vertices_per_cell*dimensions;
    float *pts = new float[npts];
    int *conn = new int[ncells*vertices_per_cell];
    int offset = 0;
    float *ghost = new float[ncells];
    for (int i = 0 ; i < ncells ; i++)
    {
        bool isGhost = false;
        pts[offset+0] = triangles[i].p1[0];
        pts[offset+1] = triangles[i].p1[1];
        pts[offset+2] = 0;
        if (pts[offset+0] == X-value of V0-from-big-triangle && 
            pts[offset+1] == Y-value of V0-from-big-triangle && 
            pts[offset+2] == Z-value of V0-from-big-triangle)
           isGhost = true;
        if (pts[offset+0] == X-value of V1-from-big-triangle && 
            pts[offset+1] == Y-value of V1-from-big-triangle && 
            pts[offset+2] == Z-value of V1-from-big-triangle)
           isGhost = true;
        if (pts[offset+0] == X-value of V2-from-big-triangle && 
            pts[offset+1] == Y-value of V2-from-big-triangle && 
            pts[offset+2] == Z-value of V2-from-big-triangle)
           isGhost = true;
        offset += 3;
        pts[offset+0] = triangles[i].p2[0];
        pts[offset+1] = triangles[i].p2[1];
        pts[offset+2] = 0;
        if (pts[offset+0] == X-value of V1-from-big-triangle && 
            pts[offset+1] == Y-value of V1-from-big-triangle && 
            pts[offset+2] == Z-value of V1-from-big-triangle)
           isGhost = true;
        offset += 3;
        pts[offset+0] = triangles[i].p3[0];
        pts[offset+1] = triangles[i].p3[1];
        pts[offset+2] = 0;
        if (pts[offset+0] == X-value of V1-from-big-triangle && 
            pts[offset+1] == Y-value of V1-from-big-triangle && 
            pts[offset+2] == Z-value of V1-from-big-triangle)
           isGhost = true;
        offset += 3;
        ghost[i] = isGhost ? 1.0 : 0.0;
    }

    for (int i = 0 ; i < 3*ncells ; i++)
    {
        conn[i] = i;
    }
    write_unstructured_mesh(filename, 0, npts/3, pts,
                            ncells, celltypes, conn, 0,
                            NULL, NULL, NULL, NULL); <<-- update for ghost
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
//  T1 will have points: p1, p4, p2 and triangles across e1 is T2, triangle across e2 is T3, and triangle across e3 is TA
//  T2 will have points: p1, p3, p4 and triangles across e1 is TB, triangle across e2 is T3, and triangle across e3 is T1
//  T3 will have points: p2, p4, p3 and triangles across e1 is T1, triangle across e2 is T2, and triangle across e3 is TC
//
            OneTriangle original_triangle = triangles[i];
	    OneTriangle *TA = original_triangle.triangle_across_e1;
	    OneTriangle *TB = original_triangle.triangle_across_e2;
	    OneTriangle *TC = original_triangle.triangle_across_e3;

            // split triangle i into three triangles
            // note: no edge flipping or Delaunay business.
            // start by replacing triangle in the current list.
            triangles[i].p3[0] = x1;
            triangles[i].p3[1] = y1;
            OneTriangle *T1 = &(triangles[i]);

            // now add two more triangles.
            OneTriangle new_triangle1;
            new_triangle1.p1[0] = original_triangle.p2[0];
            new_triangle1.p1[1] = original_triangle.p2[1];
            new_triangle1.p2[0] = original_triangle.p3[0];
            new_triangle1.p2[1] = original_triangle.p3[1];
            new_triangle1.p3[0] = x1;
            new_triangle1.p3[1] = y1;
            triangles.push_back(new_triangle1);
	    int index = triangles.size()-1;
	    OneTriangle *T3 = &(triangles[index]);

            OneTriangle new_triangle2;
            new_triangle2.p1[0] = original_triangle.p3[0];
            new_triangle2.p1[1] = original_triangle.p3[1];
            new_triangle2.p2[0] = original_triangle.p1[0];
            new_triangle2.p2[1] = original_triangle.p1[1];
            new_triangle2.p3[0] = x1;
            new_triangle2.p3[1] = y1;
            triangles.push_back(new_triangle2);
	    OneTriangle *T2 = &(triangles[index+1]);

	    if (T1 == NULL) {
              T2->triangle_across_e3 = NULL;
              T3->triangle_across_e1 = NULL;
            } else {
              T2->triangle_across_e3 = T1;
              T3->triangle_across_e1 = T1;
            }
	    if (T2 == NULL) {
	      T1->triangle_across_e1 = NULL;
	      T3->triangle_across_e2 = NULL;
	    } else {
	      T1->triangle_across_e1 = T2;
	      T3->triangle_across_e2 = T2;
	    }
	    if (T3 == NULL) {
	      T1->triangle_across_e2 = NULL;
	      T2->triangle_across_e2 = NULL;
	    } else {
	      T1->triangle_across_e2 = T3;
	      T2->triangle_across_e2 = T3;
	    }

	    if (TA == NULL)
	      T1->triangle_across_e3 = NULL;
	    else
	      T1->triangle_across_e3 = TA;
	    if (TB == NULL)
	      T2->triangle_across_e1 = NULL;
	    else
	      T2->triangle_across_e1 = TB;
	    if (TC == NULL)
	      T3->triangle_across_e3 = NULL;
	    else
	      T3->triangle_across_e3 = TC;

            break;
        }
    }
}

// Function to calculate the determinant of 3 points. This gives us the circumcircle of the triangle. If a 4th point is inside the triangle,
// the result will be negative. If the 4th point lies outside the circle, result will be positive. A result equal to zero means that the 4th
// point lies on the circle exactly. In this case, the DT of the set of points is not unique. It's like drawing two triangles in a square. 
// Whether the 3rd point and the 0th point make up the hypotenous or the 2nd and the 1st point make up the hypotenus, it's equivalent and
// therefore you could do either one and make a valid DT. Will need to call a seperate function to handle that case later.
//
// Inputs: 3 points, each with x and y coordinates (or z in case of 3D) and a 4th point with x and y coordinates. This makes up points A, B, C (of the
// triangle) and the 4th point, D.
// Output: Boolean value. True if 4th point is inside circle. This means we have to split. False if 4th point is outside circle. This means we're ok.
//
// TODO: create another case where we call a function to handle if the result is equal to zero  
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

    if (result < 0) return false; //ptD lies outside circumcircle
    else if (result > 0) return true; //ptD lies inside circumcircle

    //what about if it equals zero? Lies *on* the circle...
}

//Helper function for CircumcirlceCheck function. This helps the readability of the code. Does a simple calculation on 4 values, returns result.
float
DelaunayTriangulation::DetHelp(float pt1, float pt2, float pt3, float pt4)
{
    return ((pow((pt1 - pt2), 2)) + (pow((pt3 - pt4), 2)));
}

int main()
{
    float *pts = PointsGenerator(100, 2);
    DelaunayTriangulation DT;

    DT.Initialize(-1, -1, 2, -1, .5, 2);
    for (int i = 0 ; i < 100 ; i++)
        DT.AddPoint(pts[2*i], pts[2*i+1]);
 
    DT.Verify(); 
    DT.DelBoundingTri();

    DT.WriteOutTriangle("kristi.vtk");
    return 0;
}
