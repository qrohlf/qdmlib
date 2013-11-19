#include <Drawframework.h>

void read_object3d_from_file(FILE* f, object3d* obj) {
    fscanf(f, "%d\n", &obj->n);
    printf("Reading %d points\n", obj->n);
    for (int i=0; i<obj->n; i++) {
        fscanf(f, "%lf %lf %lf", &obj->xs[i], &obj->ys[i], &obj->zs[i]);
    }
    printf("done reading points\n");
    fscanf(f, "%d\n", &obj->num_shapes);
    printf("Reading %d shapes\n", obj->num_shapes);
    shape* shape;
    for (int i=0; i<obj->num_shapes; i++) {
        shape = &obj->shapes[i];
        shape->R = 1;
        shape->G = 1;
        shape->B = 1;
        shape->parent = obj;
        read_shape_from_file(f, shape);
    }
    center_of_mass(obj, &obj->xs[obj->n], &obj->ys[obj->n], &obj->zs[obj->n]);
    obj->center = obj->n;
    obj->n++;
    printf("done reading object\n");
}

/*
 * constructs an object from a list of points, connections, and colors
 */
void read_object2d_from_file(FILE* f, object2d *obj) {
    // Get the points
    fscanf(f, "%d", &obj->n);
    for (int i=0; i<obj->n; i++) {
        fscanf(f, "%lf %lf", &obj->xs[i], &obj->ys[i]);
    }


    // Get the shapes
    fscanf(f, "%d", &obj->num_shapes);
    for (int i=0; i<obj->num_shapes; i++) {
        fscanf(f, "%d", &obj->shapes[i].n);
        for (int j=0; j<obj->shapes[i].n; j++) {
            fscanf(f, "%d", &obj->shapes[i].vertices[j]);
        }
    }

    // Get the colors
    for (int i=0; i<obj->num_shapes; i++) {
        fscanf(f, "%lf %lf %lf", &obj->shapes[i].R, &obj->shapes[i].G, &obj->shapes[i].B);
    }
}

void read_shape_from_file(FILE* f, shape* shape) {
    fscanf(f, "%d\n", &shape->n);
    for (int i=0; i<shape->n; i++) {
        fscanf(f, "%d", &shape->vertices[i]);
    }
}

void print_object3d(object3d* obj) {
    printf("%d\n", obj->n);
    for (int i=0; i<obj->n; i++) {
        printf("\t%lf %lf %lf\n", obj->xs[i], obj->ys[i], obj->zs[i]);
    }
    printf("%d\n", obj->num_shapes);
    shape* shape;
    for (int i=0; i<obj->num_shapes; i++) {
        shape = &obj->shapes[i];
        print_shape(shape);
    }
}

void print_object2d(object2d* obj) {
    printf("%d\n", obj->n);
    for (int i=0; i<obj->n; i++) {
        printf("\t%lf %lf\n", obj->xs[i], obj->ys[i]);
    }
    printf("%d\n", obj->num_shapes);
    shape* shape;
    for (int i=0; i<obj->num_shapes; i++) {
        shape = &obj->shapes[i];
        print_shape(shape);
    }
}

void print_shape(shape* shape) {
    printf("%d\t", shape->n);
    for (int i=0; i<shape->n; i++) {
        printf("%d ", shape->vertices[i]);
    }
    printf("\n");
}

void draw_object3d(object3d* obj, double fov, double viewdistance) {
    object2d renderResult;
    render_object3d(obj, &renderResult, fov, viewdistance);
    double xs[20];
    double ys[20];
    shape* s;
    for (int i=0; i<renderResult.num_shapes; i++) {
        s = &renderResult.shapes[i];
        draw_shape(s, &renderResult, true);
        s->R = 0;
        s->G = 0;
        s->B = 0;
        draw_shape(s, &renderResult, false);
    }
}

void draw_object2d(object2d* obj) {
    for (int i=0; i<obj->num_shapes; i++) {
        draw_shape(&obj->shapes[i], obj, true);
    }
}

void draw_object2d_wireframe(object2d* obj) {
    for (int i=0; i<obj->num_shapes; i++) {
        draw_shape(&obj->shapes[i], obj, false);
    }
}

void draw_shape(shape* s, object2d* obj, int fill) {
    double x[1000];
    double y[1000];
    int point_index;
    for (int j=0; j<s->n; j++) {
        point_index = s->vertices[j];
        x[j] = obj->xs[point_index];
        y[j] = obj->ys[point_index];
    }
    G_rgb(s->R, s->G, s->B);
    if (fill) G_fill_polygon(x, y, s->n);
    else G_polygon(x, y, s->n);
}

void render_object3d(object3d* obj, object2d* result, double fov, double viewdistance) {
    double x, y, z;
    result->n = obj->n;
    for (int i=0; i<obj->n; i++) {
        x = obj->xs[i];
        y = obj->ys[i];
        z = obj->zs[i];
        result->xs[i] = x*(300/tan(fov))/(z + viewdistance) + 300;
        result->ys[i] = y*(300/tan(fov))/(z + viewdistance) + 300;
    }

    double r[3];
    double v[3];
    result->num_shapes = 0;
    shape* s;
    for (int i=0; i<obj->num_shapes; i++) {
        s = &obj->shapes[i];
        r[0] = obj->xs[s->vertices[0]];
        r[1] = obj->ys[s->vertices[0]];
        r[2] = obj->zs[s->vertices[0]] + viewdistance;
        normal_vector(obj, &obj->shapes[i], v); //not the problem
        double a = angle_between(r, v); //problem somewhere here
        if (a < M_PI/2) {
            //debug information
            //draw_vector(r, v, fov, viewdistance);
            //printf("angle for shape %d* is %f (normal: %.2f %.2f %.2f) (vec: %.2f %.2f %.2f)\n", i, a, v[0], v[1], v[2], r[0], r[1], r[2]);
            result->shapes[result->num_shapes] = obj->shapes[i];
            result->num_shapes++;
        } else {  
            //printf("angle for shape %d  is %f (normal: %.2f %.2f %.2f) (vec: %.2f %.2f %.2f)\n", i, a, v[0], v[1], v[2], r[0], r[1], r[2]);
        }
    }
}

/*
 * Return the intersection of line p1-p2
 * and line j1-j2
 */
void intersection(
    point2d p1,
    point2d p2,
    point2d j1,
    point2d j2,
    point2d* intersect) {
    printf("line 1: {%f, %f} to {%f, %f}\n", p1.x, p1.y, p2.x, p2.y);
    // Slope of line p (mp)
    double mp = (p2.y - p1.y) / (p2.x - p1.x);
    printf("Slope of line 1 is %f\n", mp);
    // Slope of line j (mj)
    double mj = (j2.y - j1.y) / (j2.x - j1.x);
    printf("Slope of line 2 is %f\n", mj);
    //p1.y = mp (p1.x) + bp
    double bp = p1.y - mp*p1.x;
    printf("Intercept of line 1 is %f\n", bp);
    //j1.y = mj (j1.x) + bj
    double bj = j1.y - mj*j1.x;
    printf("Intercept of line 2 is %f\n", bj);
    if (!isfinite(mp)) {
        //p1-p2 is a vertical line
        intersect->x = p1.x;
        intersect->y = mj * intersect->x + bj;
        return;
    }
    if (!isfinite(mj)) {
        //j1-j2 is a vertical line
        intersect->x = j1.x;
        intersect->y = mp * intersect->x + bp;
        return;
    }
    //y = mp * x + bp
    //y = mj * x + bj
    //mj*x + bj = mp*x + bp
    //mj*x - mp*x = bp - bj
    //x(mj - mp) = (bp - bj)
    intersect->x = (bp - bj)/(mj - mp);
    printf("Intersection x is %f\n", intersect->x);
    intersect->y = mp * intersect->x + bp;
    printf("Intersection y is %f\n", intersect->y);
}

void transform_object3d(object3d* obj, double mat[4][4]) {
    D3d_mat_mult_points(obj->xs, obj->ys, obj->zs, mat, obj->xs, obj->ys, obj->zs, obj->n);
}

void clip_object2d(object2d* obj, polygon* win) {
    for (int j=0; j<win->n; j++) {
        point2d p3 = {win->xs[j], win->ys[j]};
        point2d p4 = {win->xs[(j+1)%win->n], win->ys[(j+1)%win->n]};
        clip_line(obj, p3, p4);
    }
}

// Clip object2d fig by a single line
void clip_line(object2d* fig, point2d l1, point2d l2) {
    shape* s;
    for (int i=0; i<fig->num_shapes; i++) {
        clip_shape(&fig->shapes[i], fig, l1, l2);
    }
}

void clip_shape(shape* fig, object2d* parent, point2d l1, point2d l2) {
    shape out = {0, {}, fig->R, fig->G, fig->B};
    int i1, i2;
    for (int i=0; i<fig->n; i++) {
        //For each line in the shape...
        i1 = fig->vertices[i];
        i2 = fig->vertices[(i+1)%fig->n];
        point2d p1 = {parent->xs[i1], parent->ys[i1]};
        point2d p2 = {parent->xs[i2], parent->ys[i2]};
        if (isRight(p1, l1, l2)) {
            //If point2d p1 is inside the clipping line, add it to the output shape
            out.vertices[out.n] = i1;
            out.n++;
        }
        point2d intersect;
        intersection(p1, p2, l1, l2, &intersect);
        if (in_range(intersect.x, p1.x, p2.x) && in_range(intersect.y, p1.y, p2.y)) {
            
            parent->xs[parent->n] = intersect.x;
            parent->ys[parent->n] = intersect.y;
            out.vertices[out.n] = parent->n;
            parent->n++;
            out.n++;
        }
    }
    *fig = out;
}

// Use the cross product to tell if a point2d is on the right or colinear of the line bc
int isRight(point2d a, point2d b, point2d c){
     return ((b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x)) <= 0;
}

// Use the cross product to get the normal vector of a plane.
// Assumes that the shape has >= 3 points defined and that all points are on the same plane.
void normal_vector(object3d* parent, shape* shape, double r[3]) {
    double u[3]; //First vector
    u[0] = parent->xs[shape->vertices[0]] - parent->xs[shape->vertices[1]];
    u[1] = parent->ys[shape->vertices[0]] - parent->ys[shape->vertices[1]];
    u[2] = parent->zs[shape->vertices[0]] - parent->zs[shape->vertices[1]];

    double v[3]; //Second vector
    v[0] = parent->xs[shape->vertices[1]] - parent->xs[shape->vertices[2]];
    v[1] = parent->ys[shape->vertices[1]] - parent->ys[shape->vertices[2]];
    v[2] = parent->zs[shape->vertices[1]] - parent->zs[shape->vertices[2]];

    r[0] = u[1]*v[2] - u[2]*v[1];
    r[1] = u[2]*v[0] - u[0]*v[2];
    r[2] = u[0]*v[1] - u[1]*v[0];
}

// Use the inverse cosine of the dot product to find the angle between
// Two vectors
double angle_between(double i[3], double j[3]) {
    return acos(D3d_dot(i, j)/(magnitude(i) * magnitude(j)));
}

double magnitude(double i[3]) {
    return sqrt(pow(i[0], 2) + pow(i[1], 2) + pow(i[2], 2));
}

void draw_vector(double loc[3], double vec[3], double fov, double viewdistance) {
    double x1 = loc[0];
    double y1 = loc[1];
    double z1 = loc[2];
    double x2 = x1 + vec[0];
    double y2 = y1 + vec[1];
    double z2 = z1 + vec[2];
    G_rgb(0, 1, 1);
    G_fill_circle(
        x1*(300/tan(fov))/(z1) + 300,
        y1*(300/tan(fov))/(z1) + 300, 5);
    G_line(
        x1*(300/tan(fov))/(z1) + 300,
        y1*(300/tan(fov))/(z1) + 300,
        x2*(300/tan(fov))/(z2) + 300,
        y2*(300/tan(fov))/(z2) + 300);
}

void center_of_mass(object3d* obj, double* x, double* y, double* z) {
    for (int i=0; i<obj->n; i++) {
        *x += obj->xs[i];
        *y += obj->ys[i];
        *z += obj->zs[i];
    }
    *x = *x/obj->n;
    *y = *y/obj->n;
    *z = *z/obj->n;
}

/*
 * Rotate an object3d about its center of mass by x, y, z radians on the x, y z, axes (respectively)
 */
void object3d_rotate(object3d* obj, double x, double y, double z) {
    double transform[4][4], useless[4][4];
    D3d_make_identity(transform);
    int c = obj->center;
    //center_of_mass(obj, &center);
    //printf("center of mass is %f %f %f\n", center.x, center.y, center.z);
    // translate to the origin
    D3d_translate(transform, useless, -obj->xs[c], -obj->ys[c], -obj->zs[c]);
    D3d_rotate_x(transform, useless, x);
    D3d_rotate_y(transform, useless, y);
    D3d_rotate_z(transform, useless, z);
    // translate back from the origin
    D3d_translate(transform, useless, obj->xs[c], obj->ys[c], obj->zs[c]);
    // apply!
    transform_object3d(obj, transform);
}

/*
 * Scale an object3d about its center of mass by x, y, z scalefactors
 * not yet tested
 */
void object3d_scale(object3d* obj, double x, double y, double z) {
    double transform[4][4], useless[4][4];
    D3d_make_identity(transform);
    int c = obj->center;
    //center_of_mass(obj, &center);
    //printf("center of mass is %f %f %f\n", center.x, center.y, center.z);
    // translate to the origin
    D3d_translate(transform, useless, -obj->xs[c], -obj->ys[c], -obj->zs[c]);
    D3d_scale(transform, useless, x, y, z);
    // translate back from the origin
    D3d_translate(transform, useless, obj->xs[c], obj->ys[c], obj->zs[c]);
    // apply!
    transform_object3d(obj, transform);
}

void sort_shapes_by_z(object3d* obj) {
    qsort(obj->shapes, obj->num_shapes, sizeof(shape), shape_compare_distance);
}

double distance(shape* s) {
    double x = s->parent->xs[s->vertices[0]];
    double y = s->parent->ys[s->vertices[0]];
    double z = s->parent->zs[s->vertices[0]];
    return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
}

int shape_compare_distance (const void* s1, const void* s2) {
    return distance((shape*)s1) - distance((shape*)s2);
}
