#ifndef Drawframework
#define Drawframework
#include <FPT.h>
#include <qdmlib.h>
#include <D3d_matrix.h>
#define true 1
#define false 0

/*
 * shape definition. Needs to be wrapped in a 2dObject or 3dObject to work.
 */
typedef struct {
    int n;                  //Number of points in the shape
    int vertices[20];     //List of points to connect in the parent 2dObject to make the shape
    double R;
    double G;
    double B;
} shape;

/*
 * polygon definition. Standalone way to represent a single 2d polygon
 */
 typedef struct {
    int n;                  //Number of points in the shape
    double xs[10000];
    double ys[10000];
    double R;
    double G;
    double B;
} polygon;

/*
 * 2d object definition. Maxes out at 1000 points/100 shapes.
 */
typedef struct  {
    int num_shapes;
    int n;                  //Number of points in the object
    double xs[10000];
    double ys[10000];
    shape shapes[10000];
} object2d;

/*
 * 2d point definition. Useful for comparisons and return values.
 */
typedef struct {
    double x;
    double y;
} point2d;

typedef struct {
    int num_shapes;
    double xs[10000];
    double ys[10000];
    double zs[10000];
    shape shapes[10000];
    int n;
} object3d;

void read_object3d_from_file(FILE* f, object3d* obj);
void read_object2d_from_file(FILE* f, object2d *obj);
void print_object3d(object3d* obj);
void print_object2d(object2d* obj);
void print_shape(shape* shape);
void read_shape_from_file(FILE* f, shape* shape);
void draw_object3d(object3d* obj, double fov, double viewdistance);
void draw_object2d(object2d* obj);
void draw_object2d_wireframe(object2d* obj);
void render_object3d(object3d* obj, object2d* result, double fov, double viewdistance);
void transform_object3d(object3d* obj, double mat[4][4]);
void clip_object2d(object2d* obj, polygon* win);
void clip_line(object2d* fig, point2d l1, point2d l2);
void clip_shape(shape* fig, object2d* parent, point2d l1, point2d l2);
int isRight(point2d a, point2d b, point2d c);
void draw_shape(shape* s, object2d* obj, int fill);

#endif