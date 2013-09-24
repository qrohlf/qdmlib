#include "qlib.h"

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
 * Takes a polygon specified by x and y arrays and scales its 
 * coordinates according to x_scale and y_scale
 */
void scale_polygon(double* x, double* y, int vertices, double x_scale, double y_scale) {

}
/* 
 * Takes a polygon specified by x and y arrays and translates its 
 * coordinates according to x_trans and y_trans
 */
void translate_polygon(double* x, double* y, int vertices, double x_trans, double y_trans) {

}

void rotate_polygon(double* x, double* y, int vertices, double rot_radian) {

}

/*
 * Finds the bounding box for a polygon specified by x and y arrays.
 * Returns a double array containing the x_min, y_min, x_max, and y_max
 * (in that order)
 */
double* bounding_box(double* x, double* y, int vertices) {

}

double* center_point(double* x, double* y, int vertices) {
    
}


