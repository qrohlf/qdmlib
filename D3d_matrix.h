
#ifndef D3d_matrix_stuff
#define D3d_matrix_stuff

#include <stdio.h>
#include <math.h>
#include <stdlib.h>




/*

 ( x')          (x)
 ( y')  =   M * (y)
 ( z')          (z)  
 ( 1 )          (1)

instead of (x',y',z',1) = (x,y,z,1) * M  

*/




int D3d_print_mat (double a[4][4]) ;

int D3d_copy_mat (double a[4][4], double b[4][4]) ;
// a = b

int D3d_mat_mult (double res[4][4], double a[4][4], double b[4][4]) ;
// res = a * b
// this is SAFE, i.e. the user can make a call such as 
// D3d_mat_mult(p,  p,q) or D3d_mat_mult(p,  q,p) or  D3d_mat_mult(p, p,p)

int D3d_make_identity (double a[4][4]) ;
// a = I



int D3d_translate (double a[4][4], double b[4][4], 
                   double dx, double dy, double dz) ;
// a = translation*a  
// b = b*translation_inverse  



int D3d_scale (double a[4][4], double b[4][4],
               double sx, double sy, double sz) ;
// a = scale*a  
// b = b*scale_inverse  



int D3d_rotate_x (double a[4][4], double b[4][4], double radians) ;
// a = rotate*a  
// b = b*rotate_inverse  



int D3d_rotate_y (double a[4][4], double b[4][4], double radians) ;
// a = rotate*a  
// b = b*rotate_inverse  



int D3d_rotate_z (double a[4][4], double b[4][4], double radians) ;
// a = rotate*a  
// b = b*rotate_inverse  



int D3d_cs_rotate_x (double a[4][4], double b[4][4], double cs, double sn) ;
// a = rotate*a  
// b = b*rotate_inverse  


int D3d_cs_rotate_y (double a[4][4], double b[4][4], double cs, double sn) ;
// a = rotate*a  
// b = b*rotate_inverse  


int D3d_cs_rotate_z (double a[4][4], double b[4][4], double cs, double sn) ;
// a = rotate*a  
// b = b*rotate_inverse  


int D3d_negate_x (double a[4][4], double b[4][4]) ;
// negate the x....reflects in the yz-plane
// a = reflect*a 
// b = b*reflect_inverse  


int D3d_negate_y (double a[4][4], double b[4][4]) ;
// negate the y....reflects in the xz-plane
// a = reflect*a 
// b = b*reflect_inverse  


int D3d_negate_z (double a[4][4], double b[4][4]) ;
// negate the z....reflects in the xy-plane
// a = reflect*a 
// b = b*reflect_inverse  




int D3d_mat_mult_points (double *X, double *Y, double *Z,
                         double m[4][4],
                         double *x, double *y, double *z, int numpoints) ;
// |X0 X1 X2 ...|       |x0 x1 x2 ...|
// |Y0 Y1 Y2 ...| = m * |y0 y1 y2 ...|
// |Z0 Z1 Z2 ...|       |z0 z1 z2 ...|
// | 1  1  1 ...|       | 1  1  1 ...|

// SAFE, user may make a call like D3d_mat_mult_points (x,y,z, m, x,y,z, n) ;

double D3d_dot (double a[4], double b[4]);


#endif
