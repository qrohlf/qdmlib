#include <D2d_matrix.h>


/*
 * Print a 3x3 matrix to stdout
 * Tested, working
 */
int D2d_print_mat (double a[3][3]) {
    for (int i = 0; i < 3; i++) {
        printf("%10.3f %10.3f %10.3f \n", a[i][0], a[i][1], a[i][2]);
    }
    return 0;
}

/*
 * Copy the contents of b into a
 * Tested, working
 */
int D2d_copy_mat (double a[3][3], double b[3][3]) {
    for (int i = 0; i < 3; i++) {
        for(int j=0; j < 3; j++) {
            a[i][j] = b[i][j];
        }
    }
    return 0; //What are we supposed to return here?
}

/*
 * res = a.b
 * SAFE: will not modify a or b.
 * Tested, working
 */
int D2d_mat_mult (double res[3][3], double a[3][3], double b[3][3]) {
    double B[3][3];
    D2d_transpose(b, B);
    double A[3][3];
    D2d_copy_mat(A, a);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            res[i][j] = D2d_dot(A[i], B[j]);
        }
    }
    return 0; //What are we supposed to return here?
}

/*
 * Write an identity matrix to a
 * Tested, working
 */
int D2d_make_identity (double a[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            a[i][j] = (i == j);
        }
    }
    return 0; //TODO: implement this function and remove this comment
}

/*
 * Generate a translation matrix.
 * a = translation*a
 * b = b*translation_inverse
 * DANGER WILL ROBINSON: THIS FUNCTION HAS NOT BEEN TESTED
 */
int D2d_translate (double a[3][3], double b[3][3], double dx, double dy) {
    double trans[3][3] = {
        {1, 0, dx},
        {0, 1, dy},
        {0, 0, 1}
    };
    double trans_inv[3][3] = {
        {1, 0, -dx},
        {0, 1, -dy},
        {0, 0, 1}
    };
    D2d_mat_mult(a, trans, a);
    D2d_mat_mult(b, b, trans_inv);
    return 0; //TODO: implement this function and remove this comment
}

int D2d_scale (double a[3][3], double b[3][3], double sx, double sy) {
    double scale[3][3] = {
        {sx, 0, 0},
        {0, sy, 0},
        {0, 0, 1}
    };
    double scale_inv[3][3] = {
        {1/sx, 0, 0},
        {0, 1/sy, 0},
        {0, 0, 1}
    };
    D2d_mat_mult(a, scale, a);
    D2d_mat_mult(b, b, scale_inv);
    return 0; //TODO: implement this function and remove this comment
}

int D2d_rotate (double a[3][3], double b[3][3], double radians) {
    double rotate[3][3] = {
        {cos(radians), sin(radians), 0},
        {-sin(radians), cos(radians), 0},
        {0, 0, 1}
    };
    double rotate_inv[3][3] = {
        {cos(-radians), sin(-radians), 0},
        {-sin(-radians), cos(-radians), 0},
        {0, 0, 1}
    };
    D2d_mat_mult(a, rotate, a);
    D2d_mat_mult(b, b, rotate_inv);
    return 0; //TODO: implement this function and remove this comment
}

// DANGER WILL ROBINSON, NOT TESTED
int D2d_negate_x (double a[3][3], double b[3][3]) {
    double negate[3][3] = {
        {-1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };
    D2d_mat_mult(a, negate, a);
    D2d_mat_mult(b, b, negate);
    return 0; //TODO: implement this function and remove this comment
}

// DANGER WILL ROBINSON, NOT TESTED
int D2d_negate_y (double a[3][3], double b[3][3]) {
    double negate[3][3] = {
        {1, 0, 0},
        {0, -1, 0},
        {0, 0, 1}
    };
    D2d_mat_mult(a, negate, a);
    D2d_mat_mult(b, b, negate);
    return 0; //TODO: implement this function and remove this comment
}

int D2d_mat_mult_points (double *X, double *Y,
                         double m[3][3],
                         double *x, double *y, int numpoints) {
    double coord[3] = {0, 0, 1};

    for (int i = 0; i < numpoints; i++) {
        coord[0] = x[i];
        coord[1] = y[i];
        X[i] = D2d_dot(m[0], coord);
        Y[i] = D2d_dot(m[1], coord);
    }

    return 0; //TODO: implement this function and remove this comment
}

/*
 * Transpose a into b
 */
int D2d_transpose (double a[3][3], double b[3][3]) {
    for (int i = 0; i < 3; i++) {
        for(int j=0; j < 3; j++) {
            b[i][j] = a[j][i];
        }
    }
    return 0; //What are we supposed to return here?
}

double D2d_dot (double a[3], double b[3]) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
