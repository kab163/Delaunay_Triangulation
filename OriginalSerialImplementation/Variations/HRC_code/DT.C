
#include <iostream>

#include "visit_writer.c"
#include <vector>

using std::vector;

using std::cerr;
using std::endl;


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

    bool      ContainsPoint(float x, float y);
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

    void   WriteOutTriangle(char *filename);

  private:
    std::vector<OneTriangle>  triangles;
};

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
            OneTriangle original_triangle = triangles[i];
            // split triangle i into three triangles
            // note: no edge flipping or Delaunay business.
            // start by replacing triangle in the current list.
            triangles[i].p3[0] = x1;
            triangles[i].p3[1] = y1;

            // now add two more triangles.
            OneTriangle new_triangle1;
            new_triangle1.p1[0] = original_triangle.p2[0];
            new_triangle1.p1[1] = original_triangle.p2[1];
            new_triangle1.p2[0] = original_triangle.p3[0];
            new_triangle1.p2[1] = original_triangle.p3[1];
            new_triangle1.p3[0] = x1;
            new_triangle1.p3[1] = y1;
            triangles.push_back(new_triangle1);

            OneTriangle new_triangle2;
            new_triangle2.p1[0] = original_triangle.p3[0];
            new_triangle2.p1[1] = original_triangle.p3[1];
            new_triangle2.p2[0] = original_triangle.p1[0];
            new_triangle2.p2[1] = original_triangle.p1[1];
            new_triangle2.p3[0] = x1;
            new_triangle2.p3[1] = y1;
            triangles.push_back(new_triangle2);

            break;
        }
    }
}

int main()
{
    float *pts = PointsGenerator(100, 2);
/*
    for (int i = 0 ; i < 100 ; i++)
        cerr << "Pt[" << i << "]: (" << pts[2*i] << "," << pts[2*i+1] << ")" << endl;
 */
    DelaunayTriangulation DT;
    DT.Initialize(pts[0], pts[1], pts[2], pts[3], pts[4], pts[5]);
    for (int i = 3 ; i < 100 ; i++)
        DT.AddPoint(pts[2*i], pts[2*i+1]);
    DT.WriteOutTriangle("kristi.vtk");
}
