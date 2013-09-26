#include <D2d_matrix.h>
#include <FPT.h>

double x[13] = {175,225,225,300,225,225,250,200,150,175,175,100,175} ;
double y[13] = {300,300,250,225,225,200,100,175,100,200,225,225,250} ;
int  n = 13 ;


int main()
{
  double m[3][3], minv[3][3] ;
  int q ;

  G_init_graphics(600,600) ;
  G_rgb(0,0,0) ; G_clear() ;

  G_rgb(1,0,0) ;
  G_fill_polygon(x,y,n) ;
  q = G_wait_key() ;

  D2d_make_identity(m) ;
  D2d_make_identity(minv) ;
  D2d_translate(m,minv,  200,100) ;

  D2d_mat_mult_points(x,y,  m, x,y,n) ;
  G_rgb(1,1,0) ;
  G_fill_polygon(x,y,n) ;
  q = G_wait_key() ;

  D2d_mat_mult_points(x,y,  minv, x,y,n) ;
  G_rgb(0,0,1) ;
  G_fill_polygon(x,y,n) ;
  q = G_wait_key() ;


}

