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
        read_shape_from_file(f, shape);
    }
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
    draw_object2d(&renderResult);
}

void draw_object2d(object2d* obj) {
    shape *s;
    double x[1000];
    double y[1000];
    int point_index;
    for (int i=0; i<obj->num_shapes; i++) {
        s = &obj->shapes[i];
        for (int j=0; j<s->n; j++) {
            point_index = s->vertices[j];
            x[j] = obj->xs[point_index];
            y[j] = obj->ys[point_index];
        }
        G_rgb(s->R, s->G, s->B);
        G_polygon(x, y, s->n);
    }
}

void draw_shape(shape* shape, object2d* parent) {

}

void render_object3d(object3d* obj, object2d* result, double fov, double viewdistance) {
    double x, y, z;
    result->n = obj->n;
    for (int i=0; i<obj->n; i++) {
        x = obj->xs[i];
        y = obj->ys[i];
        z = obj->zs[i];
        
        // orthagonal
        // result->xs[i] = x*fov + 300;
        // result->ys[i] = y*fov + 300;
        result->xs[i] = x*fov/(z - viewdistance) + 300;
        result->ys[i] = y*fov/(z - viewdistance) + 300;
    }
    result->num_shapes = obj->num_shapes;
    for (int i=0; i<obj->num_shapes; i++) {
        result->shapes[i] = obj->shapes[i];
    }
}

void transform_object3d(object3d* obj, double mat[4][4]) {
    D3d_mat_mult_points(obj->xs, obj->ys, obj->zs, mat, obj->xs, obj->ys, obj->zs, obj->n);
}
