

#ifndef D2d_matrix_stuff
#define D2d_matrix_stuff

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


/*

 ( x')          (x)
 ( y')  =   M * (y)  
 ( 1 )          (1)


instead of (x',y',1) = (x,y,1) * M  

*/

int D2d_print_mat (double a[3][3]) ;

int D2d_copy_mat (double a[3][3], double b[3][3]) ;
// a = b

int D2d_mat_mult (double res[3][3], double a[3][3], double b[3][3]) ;
// res = a * b
// this is SAFE, i.e. the user can make a call such as 
// D2d_mat_mult(p,  p,q) or D2d_mat_mult(p,  q,p) or  D2d_mat_mult(p, p,p)

int D2d_make_identity (double a[3][3]) ;
// a = I

int D2d_translate (double a[3][3], double b[3][3], double dx, double dy) ; 
// a = translation*a  
// b = b*translation_inverse  

int D2d_scale (double a[3][3], double b[3][3], double sx, double sy) ;
// a = scale*a  
// b = b*scale_inverse  

int D2d_rotate (double a[3][3], double b[3][3], double radians) ;
// a = rotate*a  
// b = b*rotate_inverse  

int D2d_negate_x (double a[3][3], double b[3][3]) ;
// negate the x....reflects in the y-axis
// a = reflect*a 
// b = b*reflect_inverse  

int D2d_negate_y (double a[3][3], double b[3][3]) ;
// negate the y....reflects in the x-axis
// a = reflect*a 
// b = b*reflect_inverse  

int D2d_mat_mult_points (double *X, double *Y,
                         double m[3][3],
                         double *x, double *y, int numpoints) ;
// |X0 X1 X2 ...|       |x0 x1 x2 ...|
// |Y0 Y1 Y2 ...| = m * |y0 y1 y2 ...|
// | 1  1  1 ...|       | 1  1  1 ...|

// SAFE, user may make a call like D2d_mat_mult_points (x,y, m, x,y, n) ;

#endif
