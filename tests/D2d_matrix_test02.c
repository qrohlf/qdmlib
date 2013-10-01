#include <D2d_matrix.h>
#include <FPT.h>

double x[13] = {175,225,225,300,225,225,250,200,150,175,175,100,175} ;
double y[13] = {300,300,250,225,225,200,100,175,100,200,225,225,250} ;
int  n = 13 ;


int main()
{
  double a[3][3];
  printf("Identity Matrix:\n");
  D2d_make_identity(a);
  D2d_print_mat(a);

  printf("\n\nCopy of Identity Matrix:\n");
  double b[3][3];
  D2d_copy_mat(b, a);
  D2d_print_mat(b);

  printf("\n\nTranspose Test:\n");
  double c[3][3] = {
        {1, 2, 3},
        {2, 3, 5},
        {4, 2, 1}
      };
  printf("Pre-transpose:\n");
  D2d_print_mat(c);
  D2d_transpose(c, b);
  printf("Post-transpose:\n");
  D2d_print_mat(b);

  printf("\nDot product test:\n");
  double thing1[] = {1, 2, 3};
  double thing2[] = {2, 1, 2};
  // (1*2 + 2*1 + 3*2) = 10
  printf("(1*2 + 2*1 + 3*2) = %.0f \n", D2d_dot(thing1, thing2));

  printf("\nMultiplying by Identity Matrix:\n");
  double res[3][3];
  printf("Pre-multiply:\n");
  D2d_print_mat(c);
  printf("\n");
  D2d_print_mat(a);
  D2d_mat_mult(res, a, c);
  printf("Post-multiply:\n");
  D2d_print_mat(c);

  // double m[3][3], minv[3][3] ;
  // int q ;

  // G_init_graphics(600,600) ;
  // G_rgb(0,0,0) ; G_clear() ;

  // G_rgb(1,0,0) ;
  // G_fill_polygon(x,y,n) ;
  // q = G_wait_key() ;

  // D2d_make_identity(m) ;
  // D2d_make_identity(minv) ;
  // D2d_translate(m,minv,  200,100) ;

  // D2d_mat_mult_points(x,y,  m, x,y,n) ;
  // G_rgb(1,1,0) ;
  // G_fill_polygon(x,y,n) ;
  // q = G_wait_key() ;

  // D2d_mat_mult_points(x,y,  minv, x,y,n) ;
  // G_rgb(0,0,1) ;
  // G_fill_polygon(x,y,n) ;
  // q = G_wait_key() ;


}

