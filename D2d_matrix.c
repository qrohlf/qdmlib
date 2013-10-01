#include <D2d_matrix.h>



int D2d_print_mat (double a[3][3]) {
    return 0; //TODO: implement this function and remove this comment
}

/* Copy the contents of b into a */
int D2d_copy_mat (double a[3][3], double b[3][3]) {
    for (int i = 0; i < 3; i++) {
        for(int j=0; j < 3; j++) {
            a[i][j] = b[i][j];
        }
    }
    return 0; //What are we supposed to return here?
}

int D2d_mat_mult (double res[3][3], double a[3][3], double b[3][3]) {
    return 0; //TODO: implement this function and remove this comment
}

int D2d_make_identity (double a[3][3]) {
    return 0; //TODO: implement this function and remove this comment
}

int D2d_translate (double a[3][3], double b[3][3], double dx, double dy) {
    return 0; //TODO: implement this function and remove this comment
}

int D2d_scale (double a[3][3], double b[3][3], double sx, double sy) {
    return 0; //TODO: implement this function and remove this comment
}

int D2d_rotate (double a[3][3], double b[3][3], double radians) {
    return 0; //TODO: implement this function and remove this comment
}

int D2d_negate_x (double a[3][3], double b[3][3]) {
    return 0; //TODO: implement this function and remove this comment
}

int D2d_negate_y (double a[3][3], double b[3][3]) {
    return 0; //TODO: implement this function and remove this comment
}

int D2d_mat_mult_points (double *X, double *Y,
                         double m[3][3],
                         double *x, double *y, int numpoints) {

    return 0; //TODO: implement this function and remove this comment
}

int D2d_transpose (double a[3][3]) {
    double b[3][3];
    for (int i = 0; i < 3; i++) {
        for(int j=0; j < 3; j++) {
            b[i][j] = a[j][i];
        }
    }
    D2d_copy_mat(a, b);
    return 0; //What are we supposed to return here?
}

int D2d_dot (double a[3], double b[3]) {
    double dotprod = 0;
    for (int i=0; i<3; i++) {
        dotprod += a[i] * b[i];
    }
    return dotprod;
}
