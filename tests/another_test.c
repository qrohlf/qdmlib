#include <D2d_matrix.h>
#include <FPT.h>

double x[13] = {2,6,8, 2, 6} ;
double y[13] = {0,0,0,-2,-2} ;
int  n = 5 ;



int print_data()
{
  int i ;
  printf("\n") ;
  for (i = 0 ; i < n ; i++) {
    printf("%10.3lf  %10.3lf\n",x[i],y[i]) ;
  }
  printf("\n") ;
}



int testA()
{
  double m[3][3], minv[3][3] ;

  print_data() ;

  D2d_make_identity(m) ;  D2d_make_identity(minv) ;
  D2d_translate(m,minv,  2,1) ;
  D2d_rotate   (m,minv,  225*M_PI/180) ;
  D2d_scale    (m,minv,  2,2) ;
  D2d_rotate   (m,minv, -45*M_PI/180) ;  
  D2d_translate(m,minv,  0,-2) ;


  D2d_mat_mult_points(x,y,  minv, x,y,n) ;
  print_data() ;


  D2d_print_mat(m) ;
  printf("\n") ;
  D2d_print_mat(minv) ;


}






int main()
{
  testA() ;
}

