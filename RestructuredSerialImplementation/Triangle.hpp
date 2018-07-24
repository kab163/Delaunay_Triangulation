#include <stdlib.h>
#include <string.h>

#ifndef TRIANGLE_H
#define TRIANGLE_H

class Triangle {
    public:
        //MARK: Properties
        double p1[2];
        double p2[2];
        double p3[2];

        int triangle_across_e1;
        int triangle_across_e2;
        int triangle_across_e3;

        //MARK: Methods
                Triangle();
                Triangle(const Triangle&);
        double  Sign(double *, double *, double *);
        bool    ContainsPoint(double, double);

};


Triangle::Triangle() {
    triangle_across_e1 = -1;
    triangle_across_e2 = -1;
    triangle_across_e3 = -1;
}

Triangle::Triangle(const Triangle &c) {
    memcpy(p1, c.p1, sizeof(double) * 2);
    memcpy(p2, c.p2, sizeof(double) * 2);
    memcpy(p3, c.p3, sizeof(double) * 2);
    triangle_across_e1 = c.triangle_across_e1;
    triangle_across_e2 = c.triangle_across_e2;
    triangle_across_e3 = c.triangle_across_e3;
}

double Triangle::Sign(double *pA, double *pB, double *pC) {
    return (pA[0] - pC[0]) * (pB[1] - pC[1]) - (pB[0] - pC[0]) * (pA[1] - pC[1]);
}

bool Triangle::ContainsPoint(double x, double y) {
    double p4[2];
    bool b1, b2, b3;

    p4[0] = x;
    p4[1] = y;

    b1 = Sign(p4, p1, p2) < 0.0;
    b2 = Sign(p4, p2, p3) < 0.0;
    b3 = Sign(p4, p3, p1) < 0.0;

    return ((b1 == b2) && (b2 == b3));
}

#endif
