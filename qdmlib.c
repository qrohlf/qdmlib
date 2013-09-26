/*
 * QUINN, DANI, & MILES' FANTASTIC LIBRARY
 * Licensed under the WTFPL
 */

#include "qdmlib.h"
#include <math.h>

/*
 * Takes a file pointer located at the number of vertices,
 * a double array to fill with the x coordinates, and a
 * double array to fill with x coordinates. Returns the number
 * of points read, and increments the file pointer past the
 * last point read.
 */
int read_points_from_file(FILE* f, double* x, double* y) {
    int vertices;

    fscanf(f, "%d\n", &vertices);
    printf("%d points read\n", vertices);
    for (int i=0; i<vertices; i++) {
        fscanf(f, "%lf %lf", &x[i], &y[i]);
    }
    return vertices;
}

/*
 * Sort an array of doubles containing n elements in ascending order
 */
void sort_asc(double* arr, int n) {
    int min;
    double temp;
    for(int i=0; i<n; i++) {
        min = i; //min holds the location of the smallest known item
        // If there are any smaller values past i, swap them with i
        for(int j=i+1; j<n; j++) {
            if (arr[j] < arr[min]) {
                min = j;
                temp = arr[i];
                arr[i] = arr[min];
                arr[min] = temp;
            }
        }
    }
}

/*
 * Return the smallest element of an array of doubles containing n elements.
 * n must be >0.
 */
double min(double* arr, int n) {
    double min = arr[0];
    for(int i=0; i<n; i++) {
        if (arr[i] < min) min = arr[i];
    }
    return min;
}

/*
 * Return the largest element of an array of doubles containing n elements.
 * n must be >0.
 */
double max(double* arr, int n) {
    double max = arr[0];
    for(int i=0; i<n; i++) {
        if (arr[i] > max) max = arr[i];
    }
    return max;
}

/*
 * Return whether x is included in [a, b]
 */
int in_range(double x, double a, double b) {
    double range[] = {a, b};
    double lower = min(range, 2);
    double upper = max(range, 2);
    return (x <= upper && x >= lower);
}

/*
 * Takes a polygon specified by x and y arrays and scales its
 * coordinates about the polygon's center point according to x_scale and y_scale
 */
void scale_polygon(double* x, double* y, int n, double x_scale, double y_scale) {
    double centerx = (max(x, n) + min(x, n))/2;
    double centery = (max(y, n) + min(y, n))/2;
    scale_polygon_about(x, y, n, x_scale, y_scale, centerx, centery);
}

void scale_polygon_about(double* x, double* y, int n, double x_scale, double y_scale, double aboutx, double abouty) {
    for (int i=0; i<n; i++) {
        x[i] = x_scale * (x[i] - aboutx) + aboutx;
        y[i] = y_scale * (y[i] - abouty) + abouty;
    }
}

/*
 * Takes a polygon specified by x and y arrays and translates its
 * coordinates according to x_trans and y_trans
 */
void translate_polygon(double* x, double* y, int n, double x_trans, double y_trans) {
    for (int i=0; i<n; i++) {
        x[i] += x_trans;
        y[i] += y_trans;
    }
}

/*
 * Takes a polygon and rotates its coordinates about the center point
 */
void rotate_polygon(double* x, double* y, int n, double radians) {
    double centerx = (max(x, n) + min(x, n))/2;
    double centery = (max(y, n) + min(y, n))/2;
    rotate_polygon_about(x, y, n, radians, centerx, centery);
}

void rotate_polygon_about(double* x, double* y, int n, double radians, double aboutx, double abouty) {
    double h;
    for (int i=0; i<n; i++) {
        h = hyp(x[i] - aboutx, y[i] - abouty);
        pivot_point_about(&x[i], &y[i], aboutx, abouty, radians);
    }
}

/*
 * Rotate a single point about another point
 */
void pivot_point_about(double* x, double* y, double aboutx, double abouty, double radians) {
    *x -= aboutx;
    *y -= abouty;

    double xnew = *x*cos(radians) - *y*sin(radians);
    double ynew = *x*sin(radians) + *y*cos(radians);

    *x = xnew + aboutx;
    *y = ynew + abouty;
}

/*
 * apply the Pythagorean theorem
 */
double hyp(double x, double y) {
    return sqrt(x*x + y*y);
}

/*
 * Print two arrays x[] and y[] of size n
 */
void print_x_y(double *x, double *y, int n){
    int i;
    for(i = 0; i < n; i++){
        printf("x: %lf, y: %lf\n", x[i], y[i]);
    }
}

