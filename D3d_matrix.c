#include <D3d_matrix.h>


/*
 * Print a 3x3 matrix to stdout
 * Tested, working
 */
int D3d_print_mat (double a[4][4]) {
    for (int i = 0; i < 4; i++) {
        printf("%10.3f %10.3f %10.3f %10.3f \n", a[i][0], a[i][1], a[i][2], a[i][3]);
    }
    return 0;
}

/*
 * Copy the contents of b into a
 * Tested, working
 */
int D3d_copy_mat (double a[4][4], double b[4][4]) {
    for (int i = 0; i < 4; i++) {
        for(int j=0; j < 4; j++) {
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
int D3d_mat_mult (double res[4][4], double a[4][4], double b[4][4]) {
    double B[4][4];
    D3d_transpose(b, B);
    double A[4][4];
    D3d_copy_mat(A, a);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            res[i][j] = D3d_dot(A[i], B[j]);
        }
    }
    return 0; //What are we supposed to return here?
}

/*
 * Write an identity matrix to a
 * Tested, working
 */
int D3d_make_identity (double a[4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
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
int D3d_translate (double a[4][4], double b[4][4], double dx, double dy, double dz) {
    double trans[4][4] = {
        {1, 0, 0, dx},
        {0, 1, 0, dy},
        {0, 0, 1, dz},
        {0, 0, 0, 1}
    };
    double trans_inv[4][4] = {
        {1, 0, 0, -dx},
        {0, 1, 0, -dy},
        {0, 0, 1, -dz},
        {0, 0, 0, 1}
    };
    D3d_mat_mult(a, trans, a);
    D3d_mat_mult(b, b, trans_inv);
    return 0; //TODO: implement this function and remove this comment
}

int D3d_scale (double a[4][4], double b[4][4], double sx, double sy, double sz) {
    double scale[4][4] = {
        {sx, 0, 0, 0},
        {0, sy, 0, 0},
        {0, 0, sz, 0},
        {0, 0,  0, 1}
    };
    double scale_inv[4][4] = {
        {1/sx, 0, 0, 0},
        {0, 1/sy, 0, 0},
        {0, 0, 1/sz, 0},
        {0, 0,  0, 1}
    };
    D3d_mat_mult(a, scale, a);
    D3d_mat_mult(b, b, scale_inv);
    return 0; //TODO: implement this function and remove this comment
}

int D3d_rotate_x (double a[4][4], double b[4][4], double r) {
    double rotate[4][4] = {
        {1, 0, 0, 0},
        {0, cos(r), -sin(r), 0},
        {0, sin(r), cos(r), 0},
        {0, 0, 0, 1}
    };
    double rotate_inv[4][4] = {
        {1, 0, 0, 0},
        {0, cos(r), sin(r), 0},
        {0, -sin(r), cos(r), 0},
        {0, 0, 0, 1}
    };
    D3d_mat_mult(a, rotate, a);
    D3d_mat_mult(b, b, rotate_inv);
    return 0; //TODO: implement this function and remove this comment
}



int D3d_rotate_y (double a[4][4], double b[4][4], double r) {
    double rotate[4][4] = {
        {cos(r), 0, sin(r), 0},
        {0, 1, 0, 0},
        {-sin(r), 0, cos(r), 0},
        {0, 0, 0, 1}
    };
    double rotate_inv[4][4] = {
        {cos(r), 0, -sin(r), 0},
        {0, 1, 0, 0},
        {sin(r), 0, cos(r), 0},
        {0, 0, 0, 1}
    };
    D3d_mat_mult(a, rotate, a);
    D3d_mat_mult(b, b, rotate_inv);
}
// a = rotate*a  
// b = b*rotate_inverse  



int D3d_rotate_z (double a[4][4], double b[4][4], double r) {
    double rotate[4][4] = {
        {cos(r), 0, sin(r), 0},
        {0, 1, 0, 0},
        {-sin(r), 0, cos(r), 0},
        {0, 0, 0, 1}
    };
    double rotate_inv[4][4] = {
        {cos(r), 0, -sin(r), 0},
        {0, 1, 0, 0},
        {sin(r), 0, cos(r), 0},
        {0, 0, 0, 1}
    };
    D3d_mat_mult(a, rotate, a);
    D3d_mat_mult(b, b, rotate_inv);
}
// a = rotate*a  
// b = b*rotate_inverse  



int D3d_cs_rotate_x (double a[4][4], double b[4][4], double cs, double sn) {
    double rotate[4][4] = {
        {1, 0, 0, 0},
        {0, cs, -sn, 0},
        {0, sn, cs, 0},
        {0, 0, 0, 1}
    };
    double rotate_inv[4][4] = {
        {1, 0, 0, 0},
        {0, cs, sn, 0},
        {0, -sn, cs, 0},
        {0, 0, 0, 1}
    };
    D3d_mat_mult(a, rotate, a);
    D3d_mat_mult(b, b, rotate_inv);
    return 0; //TODO: implement this function and remove this comment
}
// a = rotate*a  
// b = b*rotate_inverse  


int D3d_cs_rotate_y (double a[4][4], double b[4][4], double cs, double sn) {
    double rotate[4][4] = {
        {cs, 0, sn, 0},
        {0, 1, 0, 0},
        {-sn, 0, cs, 0},
        {0, 0, 0, 1}
    };
    double rotate_inv[4][4] = {
        {cs, 0, -sn, 0},
        {0, 1, 0, 0},
        {sn, 0, cs, 0},
        {0, 0, 0, 1}
    };
    D3d_mat_mult(a, rotate, a);
    D3d_mat_mult(b, b, rotate_inv);
}
// a = rotate*a  
// b = b*rotate_inverse  


int D3d_cs_rotate_z (double a[4][4], double b[4][4], double cs, double sn) {
    double rotate[4][4] = {
        {cs, 0, sn, 0},
        {0, 1, 0, 0},
        {-sn, 0, cs, 0},
        {0, 0, 0, 1}
    };
    double rotate_inv[4][4] = {
        {cs, 0, -sn, 0},
        {0, 1, 0, 0},
        {sn, 0, cs, 0},
        {0, 0, 0, 1}
    };
    D3d_mat_mult(a, rotate, a);
    D3d_mat_mult(b, b, rotate_inv);
}
// a = rotate*a  
// b = b*rotate_inverse  

// DANGER WILL ROBINSON, NOT TESTED
int D3d_negate_x (double a[4][4], double b[4][4]) {
    double negate[4][4] = {
        {-1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    D3d_mat_mult(a, negate, a);
    D3d_mat_mult(b, b, negate);
    return 0; //TODO: implement this function and remove this comment
}

// DANGER WILL ROBINSON, NOT TESTED
int D3d_negate_y (double a[4][4], double b[4][4]) {
    double negate[4][4] = {
        {1, 0, 0, 0},
        {0, -1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    D3d_mat_mult(a, negate, a);
    D3d_mat_mult(b, b, negate);
    return 0; //TODO: implement this function and remove this comment
}

int D3d_mat_mult_points (double *X, double *Y, double *Z,
                         double m[4][4],
                         double *x, double *y, double *z, int numpoints) {
    double coord[4] = {0, 0, 0, 1};

    for (int i = 0; i < numpoints; i++) {
        coord[0] = x[i];
        coord[1] = y[i];
        coord[2] = z[i];
        X[i] = D3d_dot(m[0], coord);
        Y[i] = D3d_dot(m[1], coord);
        Z[i] = D3d_dot(m[2], coord);
    }

    return 0; //TODO: implement this function and remove this comment
}

/*
 * Transpose a into b
 */
int D3d_transpose (double a[4][4], double b[4][4]) {
    for (int i = 0; i < 4; i++) {
        for(int j=0; j < 4; j++) {
            b[i][j] = a[j][i];
        }
    }
    return 0; //What are we supposed to return here?
}

double D3d_dot (double a[4], double b[4]) {
    return 
    a[0]*b[0]
    +a[1]*b[1]
    +a[2]*b[2]
    +a[3]*b[3];
}
