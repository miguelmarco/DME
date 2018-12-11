#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

typedef uint_least64_t fq_elem;
typedef uint_least64_t u128[2];
typedef uint_least64_t u192[3];

#include "setup1.h"

int max(int a, int b) { return (a>b) ? a : b; }
int min(int a, int b) { return (a<b) ? a : b; }

/* Returns a+b in Fq */
fq_elem fq_add(fq_elem a, fq_elem b)
{
  return a ^ b;
}

/* Returns a*b in Fq */
fq_elem fq_mul(fq_elem a, fq_elem b)
{
  int i;
  fq_elem c;
  c = 0;
  for (i=0; i<48; i++)
  {
    c <<= 1;
    b <<= 1;
    if (c & (UINT64_C(1) << 48)) c ^= min_poly;
    if (b & (UINT64_C(1) << 48)) c ^= a;
  }
  return c;
}

/* Returns a^(-1) in Fq; exits if a=0 */
fq_elem fq_inv(fq_elem a)
{
  int i, j;
  fq_elem t, r, b, a2, b2;
  if (!a)
  {
    fprintf(stderr, "error: division by zero in Fq!\n");
    exit(-1);
  }
  t = 0;
  r = min_poly;
  b = 1;
  while (a)
  {
    for (i=1; a>>i; i++);
    for (j=1; r>>j; j++);
    b2 = (j>=i) ? t ^ (b << (j-i)) : t;
    a2 = (j>=i) ? r ^ (a << (j-i)) : r;
    t = b;
    r = a;
    b = b2;
    a = a2;
  }
  if (t & (UINT64_C(1) << 48)) t ^= min_poly;
  return t;
}

/* Returns a random element of Fq. */
fq_elem fq_rnd(void)
{
  fq_elem a;
  a  = (uint_least64_t)(rand() & 0xffffff);
  a <<= 24;
  a += (uint_least64_t)(rand() & 0xffffff);
  return a;
}

/* Returns a^n in Fq; 0<=n<2^48 */
fq_elem fq_pow(fq_elem b, uint_least64_t n)
{
  int j;
  fq_elem a = 1;
  for (j=47; j>=0; j--)
  {
    a = fq_mul(a, a);
    if ((n >> j) & 1) a = fq_mul(a, b);
  }
  return a;
}

/* Returns a^(2^n) in Fq */
fq_elem fq_pow_2exp(fq_elem a, unsigned int n)
{
  unsigned int i;
  fq_elem b;
  b = a;
  for (i=0; i<n; i++)
    b = fq_mul(b, b);
  return b;
}

/* Computes the product of the polynomials b[0]+b[1]*x+...+b[deg_b]*x^deg_b */
/* and c[0]+c[1]*x+...+c[deg_c]*x^deg_c, and puts the coefficients of the   */
/* result in the array a[0],...,a[deg_b+deg_c] */
void fq_poly_multiply(fq_elem *a, const fq_elem *b, const fq_elem *c, int deg_b,
int deg_c)
{
  int i, j;
  for (i=0; i<=deg_b+deg_c; i++)
  {
    a[i] = 0;
    for (j=max(0,i-deg_c); j<=min(i,deg_b); j++)
      a[i] = fq_add(a[i], fq_mul(b[j], c[i-j]));
  }
}

/* Computes the product of the matrices b (size nxm) and c (size mxl), and */
/* puts the result in a (size nxl). Matrices are linearized in row-major   */
/* order */
void fq_matrix_multiply(fq_elem *a, const fq_elem *b, const fq_elem *c, int n,
int m, int l)
{
  int i, j, k;
  fq_elem tmp;
  for (i=0; i<n; i++)
  {
    for (j=0; j<l; j++)
    {
      tmp = 0;
      for (k=0; k<m; k++)
        tmp = fq_add(tmp, fq_mul(b[i*m+k], c[k*l+j]));
      a[i*l+j] = tmp;
    }
  }
}

/* Computes the inverse of matrix b (size nxn) and puts the result in a.   */
/* Both matrices (a and b) are linearized in row-major order. Returns 0    */
/* if b in invertible and -1 otherwise (in which case a is left untouched).*/
/* n must be <= 64 */
int fq_matrix_inverse(fq_elem *a, const fq_elem *b, int n)
{
  int i, j, k;
  fq_elem tmp;
  fq_elem c[64][128];

  /* Fill c with a copy of (b|Id_n) */
  memset(c, 0, 64*128*sizeof(fq_elem));
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
      c[i][j] = b[i*n+j];
    c[i][i+n] = 1;
  }

  /* Eliminate the non-zero entries below the main diagonal (of c), pivoting */
  /* when necessary. */
  for (j=0; j<n; j++)
  {
    /* Find pivot in j-th column */
    for (i=j; i<n && !c[i][j]; i++);
    /* If no-pivot was found, the matrix is not invertible. */
    if (i==n)
      return -1;
    /* Swap i-th and j-th rows */
    for (k=j; k<2*n; k++)
    {
      tmp = c[i][k];
      c[i][k] = c[j][k];
      c[j][k] = tmp;
    }
    /* Multiply j-th row by C(j,j)^(-1) */
    tmp = fq_inv(c[j][j]);
    for (k=j; k<2*n; k++)
      c[j][k] = fq_mul(tmp, c[j][k]);
    /* Eliminate all non-zero entries below (j,j) */
    for (i=j+1; i<n; i++)
    {
      if (!c[i][j]) continue;
      for (k=j+1; k<2*n; k++)
        c[i][k] = fq_add(c[i][k], fq_mul(c[i][j], c[j][k]));
      c[i][j] = 0;
    }
  }

  /* Eliminate the non-zero entries above the main diagonal (no-pivoting is */
  /* is necessary, since we already have 1's in the diagonal. */
  for (j=n-1; j>=0; j--)
  {
    for (i=0; i<j; i++)
    {
      if (!c[i][j]) continue;
      for (k=0; k<j-1; k++)
        c[i][k] = fq_add(c[i][k], fq_mul(c[i][j], c[j][k]));
      for (k=n; k<2*n; k++)
        c[i][k] = fq_add(c[i][k], fq_mul(c[i][j], c[j][k]));
      c[i][j] = 0;
    }
  }

  /* At this point the matrix c contains (Id_n|b^-1), so we copy the second */
  /* block to array a (linearized in row-major order). */
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      a[i*n+j] = c[i][j+n];

  return 0;
}

/* Computes the kronecker product of b (size mb x nb) and c (size mc x nc), */
/* and puts the result in a (size mb*mc x nb*nc).                           */
void fq_matrix_tensor(fq_elem *a, const fq_elem *b, const fq_elem *c, int mb,
  int nb, int mc, int nc)
{
  int i, j, k, l;
  for (i=0; i<mb; i++)
    for (j=0; j<nb; j++)
      for (k=0; k<mc; k++)
        for (l=0; l<nc; l++)
          a[nb*nc*(i*mc+k)+(j*nc+l)] = fq_mul(b[i*nb+j], c[k*nc+l]);
}

/* Computes b^e coordinate-wise (BYTES = n) and put the result in a. */
void fq_vector_pow_2exp(fq_elem *a, const fq_elem *b, int n, unsigned int e)
{
  int i;
  for (i=0; i<n; i++)
    a[i] = fq_pow_2exp(b[i], e);
}

fq_elem M1[64][64];
fq_elem M1_inverse[64][64];
fq_elem M2[64][64];
fq_elem M2_inverse[64][64];
fq_elem pt_sec2pub[64][8];

void compute_monomials(fq_elem *vec, const fq_elem *pt)
{
  int i, j;
  fq_elem pt_pow[8][2], row1[8], row2[8], tmp1, tmp2, tmp3, tmp4;
  pt_pow[0][0] = fq_pow_2exp(pt[0], E11 % 48);
  pt_pow[0][1] = fq_pow_2exp(pt[0], E31 % 48);
  pt_pow[1][0] = fq_pow_2exp(pt[1], E11 % 48);
  pt_pow[1][1] = fq_pow_2exp(pt[1], E31 % 48);
  pt_pow[2][0] = fq_pow_2exp(pt[2], E12 % 48);
  pt_pow[2][1] = fq_pow_2exp(pt[2], E42 % 48);
  pt_pow[3][0] = fq_pow_2exp(pt[3], E12 % 48);
  pt_pow[3][1] = fq_pow_2exp(pt[3], E42 % 48);
  pt_pow[4][0] = fq_pow_2exp(pt[4], E23 % 48);
  pt_pow[4][1] = fq_pow_2exp(pt[4], E33 % 48);
  pt_pow[5][0] = fq_pow_2exp(pt[5], E23 % 48);
  pt_pow[5][1] = fq_pow_2exp(pt[5], E33 % 48);
  pt_pow[6][0] = fq_pow_2exp(pt[6], E24 % 48);
  pt_pow[6][1] = fq_pow_2exp(pt[6], E44 % 48);
  pt_pow[7][0] = fq_pow_2exp(pt[7], E24 % 48);
  pt_pow[7][1] = fq_pow_2exp(pt[7], E44 % 48);
  row1[0] = fq_mul(pt_pow[0][0], pt_pow[2][0]);
  row1[1] = fq_mul(pt_pow[0][0], pt_pow[3][0]);
  row1[2] = fq_mul(pt_pow[1][0], pt_pow[2][0]);
  row1[3] = fq_mul(pt_pow[1][0], pt_pow[3][0]);
  row1[4] = fq_mul(pt_pow[4][0], pt_pow[6][0]);
  row1[5] = fq_mul(pt_pow[4][0], pt_pow[7][0]);
  row1[6] = fq_mul(pt_pow[5][0], pt_pow[6][0]);
  row1[7] = fq_mul(pt_pow[5][0], pt_pow[7][0]);
  row2[0] = fq_mul(pt_pow[0][1], pt_pow[4][1]);
  row2[1] = fq_mul(pt_pow[0][1], pt_pow[5][1]);
  row2[2] = fq_mul(pt_pow[1][1], pt_pow[4][1]);
  row2[3] = fq_mul(pt_pow[1][1], pt_pow[5][1]);
  row2[4] = fq_mul(pt_pow[2][1], pt_pow[6][1]);
  row2[5] = fq_mul(pt_pow[2][1], pt_pow[7][1]);
  row2[6] = fq_mul(pt_pow[3][1], pt_pow[6][1]);
  row2[7] = fq_mul(pt_pow[3][1], pt_pow[7][1]);
  for (i=0; i<8; i++)
  {
    tmp1 = fq_pow_2exp(row1[i], F11 % 48);
    tmp3 = fq_pow_2exp(row1[i], F21 % 48);
    for (j=0; j<8; j++)
    {

      tmp2 = fq_pow_2exp(row2[j], F12 % 48);
      tmp4 = fq_pow_2exp(row2[j], F22 % 48);
      vec[8*i+j] = fq_mul(tmp1, tmp2);
      vec[8*i+j+64] = fq_mul(tmp3, tmp4);
    }
  }
}

void find_m1_m2(void)
{
  int i, j, done;
  fq_elem vec[2][64];
  done = 0;
  while (!done)
  {
    for (i=0; i<64; i++)
    {
      do
      {
        for (j=0; j<8; j++)
          pt_sec2pub[i][j] = fq_rnd();
      }
      while (!pt_sec2pub[i][1] || !pt_sec2pub[i][3] ||
             !pt_sec2pub[i][5] || !pt_sec2pub[i][7]);
      compute_monomials(&vec[0][0], &pt_sec2pub[i][0]);
      for (j=0; j<64; j++)
      {
        M1[j][i] = vec[0][j];
        M2[j][i] = vec[1][j];
      }
    }
    if (!fq_matrix_inverse(&M1_inverse[0][0], &M1[0][0], 64) &&
        !fq_matrix_inverse(&M2_inverse[0][0], &M2[0][0], 64)) done = 1;
  }
}

int main(void)
{
  int i, j;
  find_m1_m2();
  printf("const fq_elem M1_inverse[64][64] = {\n");
  for (i=0; i<64; i++)
  {
    printf("  {\n");
    for (j=0; j<64; j++)
    {
      printf("    UINT64_C(0x%012" PRIxLEAST64 ")", M1_inverse[i][j]);
      if (j != 63) printf(",");
      printf("\n");
    }
    printf("  }");
    if (i != 63) printf(",");
    printf("\n");
  }
  printf("};\n\n");

  printf("const fq_elem M2_inverse[64][64] = {\n");
  for (i=0; i<64; i++)
  {
    printf("  {\n");
    for (j=0; j<64; j++)
    {
      printf("    UINT64_C(0x%012" PRIxLEAST64 ")", M2_inverse[i][j]);
      if (j != 63) printf(",");
      printf("\n");
    }
    printf("  }");
    if (i != 63) printf(",");
    printf("\n");
  }
  printf("};\n\n");

  printf("const fq_elem pt_sec2pub[64][8] = {\n");
  for (i=0; i<64; i++)
  {
    printf("  {\n");
    for (j=0; j<8; j++)
    {
      printf("    UINT64_C(0x%012" PRIxLEAST64 ")", pt_sec2pub[i][j]);
      if (j != 7) printf(",");
      printf("\n");
    }
    printf("  }");
    if (i != 63) printf(",");
    printf("\n");
  }
  printf("};\n\n");

  return 0;
}

