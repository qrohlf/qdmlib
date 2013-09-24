/*
 * QUINN & MILES' FANTASTIC LIBRARY
 * Licensed under the WTFPL
 */

int read_points_from_file(FILE*, double*, double*);
void scale_polygon(double*, double*, int, double, double);
void translate_polygon(double*, double*, int, double, double);
void rotate_polygon(double*, double*, int, double);
double* bounding_box(double*, double*, int);
double* center_point(double*, double*, int);