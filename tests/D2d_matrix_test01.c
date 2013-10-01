#include <D2d_matrix.h>
#include <FPT.h>

int main(int argc, char const *argv[])
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

  printf("\nMultiplying by another Matrix:\n");
  printf("C:\n");
  D2d_print_mat(c);
  printf("\n");
  printf("D:\n");
    double d[3][3] = 
    {
    {1, 6, 6},
    {7, 8, 9},
    {1, 2, 4}
    };
  D2d_print_mat(c);
  D2d_mat_mult(res, c, d);
  printf("\nC*D:\n");
  D2d_print_mat(res);
    return 0;
}