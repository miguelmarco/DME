#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

typedef uint_least64_t fq_elem;
typedef fq_elem fq2_elem[2];
typedef fq_elem fq4_elem[4];

const uint_least64_t min_poly = UINT64_C(0x1000018000003);

int max(int a, int b) { return (a>b) ? a : b; }
int min(int a, int b) { return (a<b) ? a : b; }

fq_elem fq_add(fq_elem a, fq_elem b)
{
  return a ^ b;
}

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

void fq_poly_multiply(fq_elem *a, fq_elem *b, fq_elem *c, int deg_b, int deg_c)
{
  int i, j;
  for (i=0; i<=deg_b+deg_c; i++)
  {
    a[i] = 0;
    for (j=max(0,i-deg_c); j<=min(i,deg_b); j++)
      a[i] = fq_add(a[i], fq_mul(b[j], c[i-j]));
  }
}

fq_elem fq_rnd(void)
{
  fq_elem a;
  a  = (uint_least64_t)(rand() & 0xffffff);
  a <<= 24;
  a += (uint_least64_t)(rand() & 0xffffff);
  return a;
}

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

/* n must be <= 64 */
int fq_matrix_inverse(fq_elem *a, const fq_elem *b, int n)
{
  int i, j, k;
  fq_elem tmp;
  fq_elem c[64][128];

  memset(c, 0, 64*128*sizeof(fq_elem));
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
      c[i][j] = b[i*n+j];
    c[i][i+n] = 1;
  }

  for (j=0; j<n; j++)
  {
    /* Find pivot in j-th column */
    for (i=j; i<n && !c[i][j]; i++);
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

  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      a[i*n+j] = c[i][j+n];

  return 0;
}

fq_elem min_poly_a;
fq_elem min_poly_b;

void fq2_add(fq2_elem a, fq2_elem b, fq2_elem c)
{
  a[0] = fq_add(b[0], c[0]);
  a[1] = fq_add(b[1], c[1]);
}

void fq2_mul(fq2_elem a, fq2_elem b, fq2_elem c)
{
  fq_elem tmp[3];
  fq_poly_multiply(tmp, b, c, 1, 1);
  a[1] = fq_add(tmp[1], fq_mul(tmp[2], min_poly_a));
  a[0] = fq_add(tmp[0], fq_mul(tmp[2], min_poly_b));
}

fq_elem min_poly_c;
fq_elem min_poly_d;
fq_elem min_poly_e;
fq_elem min_poly_f;

void fq4_add(fq4_elem a, fq4_elem b, fq4_elem c)
{
  a[0] = fq_add(b[0], c[0]);
  a[1] = fq_add(b[1], c[1]);
  a[2] = fq_add(b[2], c[2]);
  a[3] = fq_add(b[3], c[3]);
}

void fq4_mul(fq4_elem a, fq4_elem b, fq4_elem c)
{
  fq_elem tmp[7];
  fq_poly_multiply(tmp, b, c, 3, 3);
  tmp[5] = fq_add(tmp[5], fq_mul(tmp[6], min_poly_c));
  tmp[4] = fq_add(tmp[4], fq_mul(tmp[6], min_poly_d));
  tmp[3] = fq_add(tmp[3], fq_mul(tmp[6], min_poly_e));
  tmp[2] = fq_add(tmp[2], fq_mul(tmp[6], min_poly_f));
  tmp[4] = fq_add(tmp[4], fq_mul(tmp[5], min_poly_c));
  tmp[3] = fq_add(tmp[3], fq_mul(tmp[5], min_poly_d));
  tmp[2] = fq_add(tmp[2], fq_mul(tmp[5], min_poly_e));
  tmp[1] = fq_add(tmp[1], fq_mul(tmp[5], min_poly_f));
  a[3] = fq_add(tmp[3], fq_mul(tmp[4], min_poly_c));
  a[2] = fq_add(tmp[2], fq_mul(tmp[4], min_poly_d));
  a[1] = fq_add(tmp[1], fq_mul(tmp[4], min_poly_e));
  a[0] = fq_add(tmp[0], fq_mul(tmp[4], min_poly_f));
}

int main(void)
{
  int i;
  fq_elem m2[3][3], m2i[3][3];
  fq_elem m4[7][7], m4i[7][7];
  fq2_elem x;
  fq4_elem y;

  printf("const uint_least64_t min_poly = UINT64_C(0x%07" PRIxLEAST64 ");\n",
            min_poly);

  while (1)
  {
    min_poly_a = fq_rnd();
    min_poly_b = fq_rnd();
    x[0] = 0;
    x[1] = 1;
    for (i=0; i<48; i++)
      fq2_mul(x, x, x);
    x[1] = fq_add(x[1], 1);
    memset(m2, 0, 9*sizeof(fq_elem));
    m2[0][0] = min_poly_b;
    m2[1][0] = min_poly_a;
    m2[2][0] = 1;
    m2[0][2] = m2[1][1] = x[0];
    m2[1][2] = m2[2][1] = x[1];
    if (!fq_matrix_inverse(&m2i[0][0], &m2[0][0], 3)) break;
  }
  printf("const fq_elem min_poly_a = UINT64_C(0x%06" PRIxLEAST64 ");\n",
            min_poly_a);
  printf("const fq_elem min_poly_b = UINT64_C(0x%06" PRIxLEAST64 ");\n",
            min_poly_b);

  while (1)
  {
    min_poly_c = fq_rnd();
    min_poly_d = fq_rnd();
    min_poly_e = fq_rnd();
    min_poly_f = fq_rnd();
    y[0] = 0;
    y[1] = 1;
    y[2] = 0;
    y[3] = 0;
    for (i=0; i<48; i++)
      fq4_mul(y, y, y);
    y[1] = fq_add(y[1], 1);
    memset(m4, 0, 49*sizeof(fq_elem));
    m4[0][2] = m4[1][1] = m4[2][0] = min_poly_f;
    m4[1][2] = m4[2][1] = m4[3][0] = min_poly_e;
    m4[2][2] = m4[3][1] = m4[4][0] = min_poly_d;
    m4[3][2] = m4[4][1] = m4[5][0] = min_poly_c;
    m4[4][2] = m4[5][1] = m4[6][0] = 1;
    m4[0][6] = m4[1][5] = m4[2][4] = m4[3][3] = y[0];
    m4[1][6] = m4[2][5] = m4[3][4] = m4[4][3] = y[1];
    m4[2][6] = m4[3][5] = m4[4][4] = m4[5][3] = y[2];
    m4[3][6] = m4[4][5] = m4[5][4] = m4[6][3] = y[3];
    if (!fq_matrix_inverse(&m4i[0][0], &m4[0][0], 7)) break;
  }
  printf("const fq_elem min_poly_c = UINT64_C(0x%06" PRIxLEAST64 ");\n",
            min_poly_c);
  printf("const fq_elem min_poly_d = UINT64_C(0x%06" PRIxLEAST64 ");\n",
            min_poly_d);
  printf("const fq_elem min_poly_e = UINT64_C(0x%06" PRIxLEAST64 ");\n",
            min_poly_e);
  printf("const fq_elem min_poly_f = UINT64_C(0x%06" PRIxLEAST64 ");\n",
            min_poly_f);
  printf("\n");

  return 0;
}

