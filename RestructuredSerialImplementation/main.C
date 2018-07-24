#include "Triangle.hpp"
#include "Triangulation.hpp"

#include <stdio.h>

int main(int argc, char **argv) {
    Triangle *t = new Triangle();
    printf("%d\n", t -> triangle_across_e1);
    printf("%d\n", t -> triangle_across_e2);
    printf("%d\n", t -> triangle_across_e3);
}
