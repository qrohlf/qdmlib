#include <FPT.h>
/*
 * QUINN, DANI, & MILES' FANTASTIC LIBRARY
 * Licensed under the WTFPL
 */

int read_points_from_file(FILE*, double*, double*);
void scale_polygon(double*, double*, int, double, double);
void translate_polygon(double*, double*, int, double, double);
void rotate_polygon(double*, double*, int, double);
void sort_asc(double*, int);
double max(double*, int);
double min(double*, int);
int in_range(double, double, double);
void scale_polygon(double*, double*, int, double, double);
void scale_polygon_about(double*, double*, int, double, double, double, double);
void rotate_polygon(double*, double*, int, double);
void rotate_polygon_about(double*, double*, int, double, double, double);
void pivot_point_about(double*, double*, double, double, double);
double hyp(double, double);