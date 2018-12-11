#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include "dme.h"
#include "hash.h"

/* Constants */
#include "setup.h"

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

/* Arithmetic in Fq2 = Fq[T]/<T^2+min_poly_a*T+min_poly_b> */
/* Elements in Fq2 are represented internally as a pair (r,s) of elements  */
/* in Fq. The pair (r,s) corresponds to the equivalence class of r+s*T mod */
/* the irreducible polynomial T^2+min_poly_a*T+min_poly_b in Fq[T].        */
/* In memory, the "r" component is stored at the locaction with lower      */
/* address, followed by the "s" component. This is consistent with our     */
/* convention for polynomials in the function fq_polynomial_multiply()     */
/* described above. Another reason why this is convenient is that the      */
/* isomorphism (Fq)^2 --> Fq2 that is used all along the algorithms does   */
/* not require to move anything in memory.                                 */

/* Computes b+c in Fq2 and stores the result in a (which is equivalent to  */
/* a coordinate-wise addition of elements in Fq).                          */
void fq2_add(fq2_elem a, const fq2_elem b, const fq2_elem c)
{
  a[0] = fq_add(b[0], c[0]);
  a[1] = fq_add(b[1], c[1]);
}

/* Computes b*c in Fq2 and stores the result in a (which is done by        */
/* multipying the polynomials b[0]+b[1]*T and c[0]+c[1]*T, and reducing the*/
/* resulting polynomial modulo T^2+min_poly_a*T+min_poly_b).               */
void fq2_mul(fq2_elem a, const fq2_elem b, const fq2_elem c)
{
  fq_elem tmp[3];
  fq_poly_multiply(tmp, b, c, 1, 1);
  a[1] = fq_add(tmp[1], fq_mul(tmp[2], min_poly_a));
  a[0] = fq_add(tmp[0], fq_mul(tmp[2], min_poly_b));
}

/* Computes b^n (with 0<=n<2^96) in Fq2 and stores the result in a. */
void fq2_pow(fq2_elem a, const fq2_elem b, const u128 n)
{
  int j;
  a[0] = 1;
  a[1] = 0;
  for (j=31; j>=0; j--)
  {
    fq2_mul(a, a, a);
    if ((n[1] >> j) & 1) fq2_mul(a, a, b);
  }
  for (j=63; j>=0; j--)
  {
    fq2_mul(a, a, a);
    if ((n[0] >> j) & 1) fq2_mul(a, a, b);
  }
}

/* Returns 1 if a=0 in Fq2, and 0 otherwise. */
int fq2_is0(const fq2_elem a)
{
  return (!a[0] && !a[1]);
}

/* Computes a = b^(2^n) in Fq2 */
void fq2_pow_2exp(fq2_elem a, fq2_elem b, unsigned int n)
{
  unsigned int i;
  memcpy(a, b, sizeof(fq2_elem));
  for (i=0; i<n; i++)
    fq2_mul(a, a, a);
}

/* Arithmetic in Fq4 = Fq[S]/<S^4+min_poly_c*S^3+min_poly_d*S^2+min_poly_e*S+min_poly_f> */
void fq4_add(fq4_elem a, const fq4_elem b, const fq4_elem c)
{
  a[0] = fq_add(b[0], c[0]);
  a[1] = fq_add(b[1], c[1]);
  a[2] = fq_add(b[2], c[2]);
  a[3] = fq_add(b[3], c[3]);
}

void fq4_mul(fq4_elem a, const fq4_elem b, const fq4_elem c)
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

void fq4_pow(fq4_elem a, const fq4_elem b, const u192 n)
{
  int i, j;
  a[0] = 1;
  a[1] = 0;
  a[2] = 0;
  a[3] = 0;
  if (n[0] || n[1] || n[2])
  {
    for (i=2; i>=0; i--)
      for (j=63; j>=0; j--)
      {
        fq4_mul(a, a, a);
        if ((n[i] >> j) & 1) fq4_mul(a, a, b);
      }
  }
}

/* Returns 1 if a=0 in Fq4, and 0 otherwise. */
int fq4_is0(const fq4_elem a)
{
  return (!a[0] && !a[1] && !a[2] && !a[3]);
}

/* Computes a = b^(2^n) in Fq4 */
void fq4_pow_2exp(fq4_elem a, fq4_elem b, unsigned int n)
{
  unsigned int i;
  memcpy(a, b, sizeof(fq4_elem));
  for (i=0; i<n; i++)
    fq4_mul(a, a, a);
}

/* API (internal) */
int fq_encrypt_with_skey(fq_elem *ct, const fq_elem *pt, const secret_key *skey)
{
  fq_elem y[8], z[8], w[8], v[8];
  fq2_elem tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  fq4_elem aux1, aux2, aux3, aux4;

  if (fq2_is0(&pt[0]) || fq2_is0(&pt[2]) || fq2_is0(&pt[4]) || fq2_is0(&pt[6]))
    return -1;

  /* L1 */
  fq_matrix_multiply(&y[0], &skey->L1[0][0][0], &pt[0], 2, 2, 1);
  fq_matrix_multiply(&y[2], &skey->L1[1][0][0], &pt[2], 2, 2, 1);
  fq_matrix_multiply(&y[4], &skey->L1[2][0][0], &pt[4], 2, 2, 1);
  fq_matrix_multiply(&y[6], &skey->L1[3][0][0], &pt[6], 2, 2, 1);

  /* F1 */
  fq2_pow(tmp1, &y[0], F1[0][0]);
  fq2_pow(tmp2, &y[2], F1[0][1]);
  fq2_pow(tmp3, &y[4], F1[1][2]);
  fq2_pow(tmp4, &y[6], F1[1][3]);
  fq2_pow(tmp5, &y[0], F1[2][0]);
  fq2_pow(tmp6, &y[4], F1[2][2]);
  fq2_pow(tmp7, &y[2], F1[3][1]);
  fq2_pow(tmp8, &y[6], F1[3][3]);
  fq2_mul(&z[0], tmp1, tmp2);
  fq2_mul(&z[2], tmp3, tmp4);
  fq2_mul(&z[4], tmp5, tmp6);
  fq2_mul(&z[6], tmp7, tmp8);

  /* L2 */
  fq_matrix_multiply(&w[0], &skey->L2[0][0][0], &z[0], 4, 4, 1);
  fq_matrix_multiply(&w[4], &skey->L2[1][0][0], &z[4], 4, 4, 1);

  /* F2 */
  fq4_pow(aux1, &w[0], F2[0][0]);
  fq4_pow(aux2, &w[4], F2[0][1]);
  fq4_pow(aux3, &w[0], F2[1][0]);
  fq4_pow(aux4, &w[4], F2[1][1]);
  fq4_mul(&v[0], aux1, aux2);
  fq4_mul(&v[4], aux3, aux4);

  /* L3 */
  fq_matrix_multiply(&ct[0], &skey->L3[0][0][0], &v[0], 4, 4, 1);
  fq_matrix_multiply(&ct[4], &skey->L3[1][0][0], &v[4], 4, 4, 1);

  return 0;
}

int fq_decrypt(fq_elem *pt, const fq_elem *ct, const secret_key *skey)
{
  fq_elem y[8], z[8], w[8], v[8];
  fq2_elem tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  fq2_elem tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
  fq4_elem aux1, aux2, aux3, aux4;

  /* L3^(-1) */
  fq_matrix_multiply(&v[0], &skey->L3_inverse[0][0][0], &ct[0], 4, 4, 1);
  fq_matrix_multiply(&v[4], &skey->L3_inverse[1][0][0], &ct[4], 4, 4, 1);

  /* F2^(-1) */
  fq4_pow(aux1, &v[0], F2_inverse[0][0]);
  fq4_pow(aux2, &v[4], F2_inverse[0][1]);
  fq4_pow(aux3, &v[0], F2_inverse[1][0]);
  fq4_pow(aux4, &v[4], F2_inverse[1][1]);
  fq4_mul(&w[0], aux1, aux2);
  fq4_mul(&w[4], aux3, aux4);

  /* L2^(-1) */
  fq_matrix_multiply(&z[0], &skey->L2_inverse[0][0][0], &w[0], 4, 4, 1);
  fq_matrix_multiply(&z[4], &skey->L2_inverse[1][0][0], &w[4], 4, 4, 1);

  /* F1^(-1) */
  fq2_pow(tmp1, &z[0], F1_inverse[0][0]);
  fq2_pow(tmp2, &z[2], F1_inverse[0][1]);
  fq2_pow(tmp3, &z[4], F1_inverse[0][2]);
  fq2_pow(tmp4, &z[6], F1_inverse[0][3]);
  fq2_pow(tmp5, &z[0], F1_inverse[1][0]);
  fq2_pow(tmp6, &z[2], F1_inverse[1][1]);
  fq2_pow(tmp7, &z[4], F1_inverse[1][2]);
  fq2_pow(tmp8, &z[6], F1_inverse[1][3]);
  fq2_pow(tmp9, &z[0], F1_inverse[2][0]);
  fq2_pow(tmp10, &z[2], F1_inverse[2][1]);
  fq2_pow(tmp11, &z[4], F1_inverse[2][2]);
  fq2_pow(tmp12, &z[6], F1_inverse[2][3]);
  fq2_pow(tmp13, &z[0], F1_inverse[3][0]);
  fq2_pow(tmp14, &z[2], F1_inverse[3][1]);
  fq2_pow(tmp15, &z[4], F1_inverse[3][2]);
  fq2_pow(tmp16, &z[6], F1_inverse[3][3]);
  fq2_mul(&y[0], tmp1, tmp2);
  fq2_mul(&y[0], &y[0], tmp3);
  fq2_mul(&y[0], &y[0], tmp4);
  fq2_mul(&y[2], tmp5, tmp6);
  fq2_mul(&y[2], &y[2], tmp7);
  fq2_mul(&y[2], &y[2], tmp8);
  fq2_mul(&y[4], tmp9, tmp10);
  fq2_mul(&y[4], &y[4], tmp11);
  fq2_mul(&y[4], &y[4], tmp12);
  fq2_mul(&y[6], tmp13, tmp14);
  fq2_mul(&y[6], &y[6], tmp15);
  fq2_mul(&y[6], &y[6], tmp16);

  /* L1^(-1) */
  fq_matrix_multiply(&pt[0], &skey->L1_inverse[0][0][0], &y[0], 2, 2, 1);
  fq_matrix_multiply(&pt[2], &skey->L1_inverse[1][0][0], &y[2], 2, 2, 1);
  fq_matrix_multiply(&pt[4], &skey->L1_inverse[2][0][0], &y[4], 2, 2, 1);
  fq_matrix_multiply(&pt[6], &skey->L1_inverse[3][0][0], &y[6], 2, 2, 1);

  if (fq2_is0(&pt[0]) || fq2_is0(&pt[2]) || fq2_is0(&pt[4]) || fq2_is0(&pt[6]))
    return -1;

  return 0;
}

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

int fq_encrypt(fq_elem *ct, const fq_elem *pt, const public_key *pkey)
{
  fq_elem vec[2][64];
  if (fq2_is0(&pt[0]) || fq2_is0(&pt[2]) || fq2_is0(&pt[4]) || fq2_is0(&pt[6]))
    return -1;
  compute_monomials(&vec[0][0], pt);
  fq_matrix_multiply(&ct[0], &pkey->coeffs1[0][0], &vec[0][0], 4, 64, 1);
  fq_matrix_multiply(&ct[4], &pkey->coeffs2[0][0], &vec[1][0], 4, 64, 1);
  return 0;
}

int fq_skey_to_pkey(public_key *pkey, const secret_key *skey)
{
  int i, j;
  fq_elem CT1[4][64];
  fq_elem CT2[4][64];
  fq_elem ct[8];
  for (i=0; i<64; i++)
  {
    fq_encrypt_with_skey(ct, &pt_sec2pub[i][0], skey);
    for (j=0; j<4; j++)
    {
      CT1[j][i] = ct[j];
      CT2[j][i] = ct[j+4];
    }
  }
  fq_matrix_multiply(&pkey->coeffs1[0][0], &CT1[0][0], &M1_inverse[0][0], 4, 64, 64);
  fq_matrix_multiply(&pkey->coeffs2[0][0], &CT2[0][0], &M2_inverse[0][0], 4, 64, 64);
  return 0;
}

/* Serialization */
void serialize_fq_elem(unsigned char **p, fq_elem a)
{
  int i;
  for (i=5; i>=0; i--)
    *((*p)++) = (a >> (8*i)) & 0xff;
}

void serialize_skey(unsigned char *p, const secret_key *skey)
{
  int i, j, k;
  for (i=0; i<4; i++)
    for (j=0; j<2; j++)
      for (k=0; k<2; k++)
        serialize_fq_elem(&p, skey->L1[i][j][k]);
  for (i=0; i<2; i++)
    for (j=0; j<4; j++)
      for (k=0; k<4; k++)
        serialize_fq_elem(&p, skey->L2[i][j][k]);
  for (i=0; i<2; i++)
    for (j=0; j<4; j++)
      for (k=0; k<4; k++)
        serialize_fq_elem(&p, skey->L3[i][j][k]);
}

void serialize_pkey(unsigned char *p, const public_key *pkey)
{
  int i, j;
  for (i=0; i<4; i++)
    for (j=0; j<64; j++)
      serialize_fq_elem(&p, pkey->coeffs1[i][j]);
  for (i=0; i<4; i++)
    for (j=0; j<64; j++)
      serialize_fq_elem(&p, pkey->coeffs2[i][j]);
}

void serialize_text(unsigned char *p, const fq_elem *text)
{
  int i;
  for (i=0; i<8; i++)
    serialize_fq_elem(&p, text[i]);
}

fq_elem parse_fq_elem(const unsigned char **p)
{
  int i;
  fq_elem a;
  a = 0;
  for (i=0; i<6; i++)
  {
    a <<= 8;
    a += (unsigned char) *((*p)++);
  }
  return a;
}

int parse_skey(secret_key *skey, const unsigned char *p)
{
  int i, j, k;
  for (i=0; i<4; i++)
  {
    for (j=0; j<2; j++)
      for (k=0; k<2; k++)
        skey->L1[i][j][k] = parse_fq_elem(&p);
    if (fq_matrix_inverse(&skey->L1_inverse[i][0][0], &skey->L1[i][0][0], 2))
      return -1;
  }
  for (i=0; i<2; i++)
  {
    for (j=0; j<4; j++)
      for (k=0; k<4; k++)
        skey->L2[i][j][k] = parse_fq_elem(&p);
    if (fq_matrix_inverse(&skey->L2_inverse[i][0][0], &skey->L2[i][0][0], 4))
      return -1;
  }
  for (i=0; i<2; i++)
  {
    for (j=0; j<4; j++)
      for (k=0; k<4; k++)
        skey->L3[i][j][k] = parse_fq_elem(&p);
    if (fq_matrix_inverse(&skey->L3_inverse[i][0][0], &skey->L3[i][0][0], 4))
      return -1;
  }
  return 0;
}

int parse_pkey(public_key *pkey, const unsigned char *p)
{
  int i, j;
  for (i=0; i<4; i++)
    for (j=0; j<64; j++)
      pkey->coeffs1[i][j] = parse_fq_elem(&p);
  for (i=0; i<4; i++)
    for (j=0; j<64; j++)
      pkey->coeffs2[i][j] = parse_fq_elem(&p);
  return 0;
}

void parse_text(fq_elem *text, const unsigned char *p)
{
  int i;
  for (i=0; i<8; i++)
    text[i] = parse_fq_elem(&p);
}

/* API (external) */
int dme_encrypt(unsigned char *ct_hex, const unsigned char *pt_hex, const
  unsigned char *pkey_hex, const unsigned char *padding)
{
  int i, j;
  struct sha2_256_ctx ct1, ct2;
  unsigned char x1[32], x2[32];
  unsigned char pt_hex_padded[TEXT_BYTES];
  fq_elem ct[8], pt[8];
  public_key pkey;

  sha2_256_init_ctx(&ct1);
  sha2_256_process_bytes((char *)padding, PADDING, &ct1);
  sha2_256_finish_ctx(&ct1, (char *)x1);
  for (i=0; i<PT_BYTES; i++)
    x1[i] ^= pt_hex[i];

  sha2_256_init_ctx(&ct2);
  sha2_256_process_bytes((char *)x1, PT_BYTES, &ct2);
  sha2_256_finish_ctx(&ct2, (char *)x2);
  for (i=0; i<PADDING; i++)
     x2[i] ^= padding[i];

  parse_pkey(&pkey, pkey_hex);
  for (i=0; i<8; i++)
  {
    for (j=0; j<4; j++)
      pt_hex_padded[6*i+j] = x1[4*i+j];
    pt_hex_padded[6*i+4] = x2[2*i];
    pt_hex_padded[6*i+5] = x2[2*i+1];
  }
  parse_text(pt, pt_hex_padded);
  if (fq_encrypt(ct, pt, &pkey)) return -1;
  serialize_text(ct_hex, ct);
  return 0;
}

int dme_decrypt(unsigned char *pt_hex, const unsigned char *ct_hex, const
  unsigned char *skey_hex, unsigned char *padding)
{
  int i, j, e;
  struct sha2_256_ctx ct1, ct2;
  unsigned char x1[32], x2[32], x3[32], x4[32];
  unsigned char pt_hex_padded[TEXT_BYTES];
  fq_elem ct[8], pt[8];
  secret_key skey;
  
  parse_skey(&skey, skey_hex);
  parse_text(ct, ct_hex);
  e = fq_decrypt(pt, ct, &skey);
  serialize_text(pt_hex_padded, pt);
  for (i=0; i<8; i++)
  {
    for (j=0; j<4; j++)
      x1[4*i+j] = pt_hex_padded[6*i+j];
    x2[2*i] = pt_hex_padded[6*i+4];
    x2[2*i+1] = pt_hex_padded[6*i+5];
  }

  sha2_256_init_ctx(&ct1);
  sha2_256_process_bytes((char *)x1, PT_BYTES, &ct1);
  sha2_256_finish_ctx(&ct1, (char *)x3);
  for (i=0; i<PADDING; i++)
    x4[i] = x2[i] ^ x3[i];

  sha2_256_init_ctx(&ct2);
  sha2_256_process_bytes((char *)x4, PADDING, &ct2);
  sha2_256_finish_ctx(&ct2, (char *)x3);
  for (i=0; i<PT_BYTES; i++)
     x1[i] ^= x3[i];
  
  for (i=0; i<8; i++)
  {
    for (j=0; j<4; j++)
      pt_hex[4*i+j] = x1[4*i+j];
    padding[2*i] = x4[2*i];
    padding[2*i+1] = x4[2*i+1];
  }
  return e;
}

int dme_encrypt_with_skey(unsigned char *ct_hex, const unsigned char *pt_hex,
  const unsigned char *skey_hex, const unsigned char *padding)
{
  int i, j;
  struct sha2_256_ctx ct1, ct2;
  unsigned char x1[32], x2[32];
  unsigned char pt_hex_padded[TEXT_BYTES];
  fq_elem ct[8], pt[8];
  secret_key skey;

  sha2_256_init_ctx(&ct1);
  sha2_256_process_bytes((char *)padding, PADDING, &ct1);
  sha2_256_finish_ctx(&ct1, (char *)x1);
  for (i=0; i<PT_BYTES; i++)
    x1[i] ^= pt_hex[i];

  sha2_256_init_ctx(&ct2);
  sha2_256_process_bytes((char *)x1, PT_BYTES, &ct2);
  sha2_256_finish_ctx(&ct2, (char *)x2);
  for (i=0; i<PADDING; i++)
     x2[i] ^= padding[i];

  parse_skey(&skey, skey_hex);
  for (i=0; i<8; i++)
  {
    for (j=0; j<4; j++)
      pt_hex_padded[6*i+j] = x1[4*i+j];
    pt_hex_padded[6*i+4] = x2[2*i];
    pt_hex_padded[6*i+5] = x2[2*i+1];
  }
  parse_text(pt, pt_hex_padded);
  if (fq_encrypt_with_skey(ct, pt, &skey)) return -1;
  serialize_text(ct_hex, ct);
  return 0;
}

int dme_skey_to_pkey(unsigned char *pkey_hex, const unsigned char *skey_hex)
{
  secret_key skey;
  public_key pkey;
  parse_skey(&skey, skey_hex);
  if (fq_skey_to_pkey(&pkey, &skey)) return -1;
  serialize_pkey(pkey_hex, &pkey);
  return 0;
}

#ifdef NIST

int crypto_kem_keypair(unsigned char *pk, unsigned char *sk)
{
  do
  {
    randombytes(sk, CRYPTO_SECRETKEYBYTES);
  }
  while (dme_skey_to_pkey(pk, sk) == -1);
  return 0;
}

int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk)
{
  unsigned char pad[PADDING];
  do
  {
    randombytes(ss, CRYPTO_BYTES);
    randombytes(pad, PADDING);
  }
  while (dme_encrypt(ct, ss, pk, pad) == -1);
  return 0;
}

int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned
  char *sk)
{
  unsigned char pad[PADDING];
  return dme_decrypt(ss, ct, sk, pad);
}

#endif
