#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "bigint.h"

bigint F1[4][4], F1_inverse[4][4];
bigint F2[2][2], F2_inverse[2][2];
bigint EXPS1[64][8], EXPS2[64][8];

int init_constants(int E11, int E12, int E23, int E24, int E31, int E33,
       int E42, int E44, int F11, int F12, int F21, int F22)
{
  int i, j, k;
  bigint ROW1[8][8], ROW2[8][8], tmp1, tmp2, q1m1, q2m1, q4m1;
  if (E11 < 0 || E12 < 0 || E23 < 0 || E24 < 0 || E31 < 0 || E33 < 0) exit(-1);
  if (E42 < 0 || E44 < 0 || F11 < 0 || F12 < 0 || F21 < 0 || F22 < 0) exit(-1);
  /* F1 */
  bigint_powi(F1[0][0], 2, E11 % 96);
  bigint_powi(F1[0][1], 2, E12 % 96);
  bigint_seti(F1[0][2], 0);
  bigint_seti(F1[0][3], 0);
  bigint_seti(F1[1][0], 0);
  bigint_seti(F1[1][1], 0);
  bigint_powi(F1[1][2], 2, E23 % 96);
  bigint_powi(F1[1][3], 2, E24 % 96);
  bigint_powi(F1[2][0], 2, E31 % 96);
  bigint_seti(F1[2][1], 0);
  bigint_powi(F1[2][2], 2, E33 % 96);
  bigint_seti(F1[2][3], 0);
  bigint_seti(F1[3][0], 0);
  bigint_powi(F1[3][1], 2, E42 % 96);
  bigint_seti(F1[3][2], 0);
  bigint_powi(F1[3][3], 2, E44 % 96);
  /* F1_inverse */
  bigint_powi(q2m1, 2, 96);
  bigint_subi(q2m1, q2m1, 1);
  if (bigint_matrix_inverse(&F1_inverse[0][0], &F1[0][0], q2m1, 4))
    return -1;
  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      if (bigint_sgn(F1_inverse[i][j]) < 0)
        bigint_add(F1_inverse[i][j], F1_inverse[i][j], q2m1);
  /* F2 */
  bigint_powi(F2[0][0], 2, F11 % 192);
  bigint_powi(F2[0][1], 2, F12 % 192);
  bigint_powi(F2[1][0], 2, F21 % 192);
  bigint_powi(F2[1][1], 2, F22 % 192);
  /* F2_inverse */
  bigint_powi(q4m1, 2, 192);
  bigint_subi(q4m1, q4m1, 1);
  if (bigint_matrix_inverse(&F2_inverse[0][0], &F2[0][0], q4m1, 2))
    return -1;
  for (i=0; i<2; i++)
    for (j=0; j<2; j++)
      if (bigint_sgn(F2_inverse[i][j]) < 0)
        bigint_add(F2_inverse[i][j], F2_inverse[i][j], q4m1);
  /* ROW1 */
  memset(ROW1, 0, 8*8*sizeof(bigint));
  bigint_set(ROW1[0][0], F1[0][0]);
  bigint_set(ROW1[0][2], F1[0][1]);
  bigint_set(ROW1[1][0], F1[0][0]);
  bigint_set(ROW1[1][3], F1[0][1]);
  bigint_set(ROW1[2][1], F1[0][0]);
  bigint_set(ROW1[2][2], F1[0][1]);
  bigint_set(ROW1[3][1], F1[0][0]);
  bigint_set(ROW1[3][3], F1[0][1]);
  bigint_set(ROW1[4][4], F1[1][2]);
  bigint_set(ROW1[4][6], F1[1][3]);
  bigint_set(ROW1[5][4], F1[1][2]);
  bigint_set(ROW1[5][7], F1[1][3]);
  bigint_set(ROW1[6][5], F1[1][2]);
  bigint_set(ROW1[6][6], F1[1][3]);
  bigint_set(ROW1[7][5], F1[1][2]);
  bigint_set(ROW1[7][7], F1[1][3]);
  /* ROW2 */
  memset(ROW2, 0, 8*8*sizeof(bigint));
  bigint_set(ROW2[0][0], F1[2][0]);
  bigint_set(ROW2[0][4], F1[2][2]);
  bigint_set(ROW2[1][0], F1[2][0]);
  bigint_set(ROW2[1][5], F1[2][2]);
  bigint_set(ROW2[2][1], F1[2][0]);
  bigint_set(ROW2[2][4], F1[2][2]);
  bigint_set(ROW2[3][1], F1[2][0]);
  bigint_set(ROW2[3][5], F1[2][2]);
  bigint_set(ROW2[4][2], F1[3][1]);
  bigint_set(ROW2[4][6], F1[3][3]);
  bigint_set(ROW2[5][2], F1[3][1]);
  bigint_set(ROW2[5][7], F1[3][3]);
  bigint_set(ROW2[6][3], F1[3][1]);
  bigint_set(ROW2[6][6], F1[3][3]);
  bigint_set(ROW2[7][3], F1[3][1]);
  bigint_set(ROW2[7][7], F1[3][3]);
  /* EXPS1 and EXPS2 */
  bigint_powi(q1m1, 2, 48);
  bigint_subi(q1m1, q1m1, 1);
  for (i=0; i<8; i++)
    for (j=0; j<8; j++)
      for (k=0; k<8; k++)
      {
        bigint_mul(tmp1, F2[0][0], ROW1[i][k]);
        bigint_mul(tmp2, F2[0][1], ROW2[j][k]);
        bigint_add(EXPS1[8*i+j][k], tmp1, tmp2);
        if (bigint_sgn(EXPS1[8*i+j][k]))
        {
          bigint_qr(tmp1, EXPS1[8*i+j][k], EXPS1[8*i+j][k], q1m1);
          if (!bigint_sgn(EXPS1[8*i+j][k]))
            bigint_set(EXPS1[8*i+j][k], q1m1);
        }
        bigint_mul(tmp1, F2[1][0], ROW1[i][k]);
        bigint_mul(tmp2, F2[1][1], ROW2[j][k]);
        bigint_add(EXPS2[8*i+j][k], tmp1, tmp2);
        if (bigint_sgn(EXPS2[8*i+j][k]))
        {
          bigint_qr(tmp1, EXPS2[8*i+j][k], EXPS2[8*i+j][k], q1m1);
          if (!bigint_sgn(EXPS2[8*i+j][k]))
            bigint_set(EXPS2[8*i+j][k], q1m1);
        }
      }
  return 0;
}

int num_bits(bigint a)
{
  int i, j, n;
  n = 0;
  for (i=NDIG-1; i>=0; i--)
    for (j=31; j>=0; j--)
      if ((a[i] >> j) & 1) n++;
  return n;
}

int main(void)
{
  int i, j, k, l, ok;
  int E11, E12, E23, E24, E31, E33, E42, E44, F11, F12, F21, F22, val;
  int best_E11, best_E12, best_E23, best_E24, best_E31, best_E33;
  int best_E42, best_E44, best_F11, best_F12, best_F21, best_F22, best;
  bigint base, r;

  best_E11 = best_E12 = best_E23 = best_E24 = best_E31 = best_E33 = 0;
  best_E42 = best_E44 = best_F11 = best_F12 = best_F21 = best_F22 = 0;
  best = 0;
  for (i=0; i<10000; i++)
  {
    E11 = rand() % 96;
    E12 = rand() % 96;
    E23 = rand() % 96;
    E24 = rand() % 96;
    E31 = rand() % 96;
    E33 = rand() % 96;
    E42 = rand() % 96;
    E44 = rand() % 96;
    F11 = rand() % 192;
    F12 = rand() % 192;
    F21 = rand() % 192;
    F22 = rand() % 192;
    val = init_constants(E11, E12, E23, E24, E31, E33, E42, E44,
                         F11, F12, F21, F22);
    if (val == -1) continue;
    ok = 1;
    for (j=0; j<64 && ok; j++)
    {
      for (k=j+1; k<64; k++)
      {
        for (l=0; l<8; l++)
          if (bigint_cmp(EXPS1[j][l], EXPS1[k][l])) break;
        if (l == 8) ok = 0;
      }
    }
    if (!ok) continue;
    ok = 1;
    for (j=0; j<64 && ok; j++)
    {
      for (k=j+1; k<64; k++)
      {
        for (l=0; l<8; l++)
          if (bigint_cmp(EXPS2[j][l], EXPS2[k][l])) break;
        if (l == 8) ok = 0;
      }
    }
    if (!ok) continue;
    val = 0;
    for (j=0; j<4; j++)
      for (k=0; k<4; k++)
        val += num_bits(F1_inverse[j][k]);
    for (j=0; j<2; j++)
      for (k=0; k<2; k++)
        val += num_bits(F2_inverse[j][k]);
    if (val > best)
    {
      best = val;
      best_E11 = E11;
      best_E12 = E12;
      best_E23 = E23;
      best_E24 = E24;
      best_E31 = E31;
      best_E33 = E33;
      best_E42 = E42;
      best_E44 = E44;
      best_F11 = F11;
      best_F12 = F12;
      best_F21 = F21;
      best_F22 = F22;
    }
  }

  printf("const unsigned int E11 = %d;\n", best_E11);
  printf("const unsigned int E12 = %d;\n", best_E12);
  printf("const unsigned int E23 = %d;\n", best_E23);
  printf("const unsigned int E24 = %d;\n", best_E24);
  printf("const unsigned int E31 = %d;\n", best_E31);
  printf("const unsigned int E33 = %d;\n", best_E33);
  printf("const unsigned int E42 = %d;\n", best_E42);
  printf("const unsigned int E44 = %d;\n", best_E44);
  printf("\n");
  
  printf("const unsigned int F11 = %d;\n", best_F11);
  printf("const unsigned int F12 = %d;\n", best_F12);
  printf("const unsigned int F21 = %d;\n", best_F21);
  printf("const unsigned int F22 = %d;\n", best_F22);
  printf("\n");

    init_constants(best_E11, best_E12, best_E23, best_E24, best_E31, best_E33,
                 best_E42, best_E44, best_F11, best_F12, best_F21, best_F22);

  bigint_powi(base, 2, 64);

  printf("const u128 F1[4][4] =\n{\n  ");
  for (i=0; i<4; i++)
  {
    printf("{\n    ");
    for (j=0; j<4; j++)
    {
      printf("{ ");
      printf("UINT64_C(0x");
      bigint_qr(F1[i][j], r, F1[i][j], base);
      bigint_print(r, 16);
      printf("), ");
      printf("UINT64_C(0x");
      bigint_qr(F1[i][j], r, F1[i][j], base);
      bigint_print(r, 16);
      printf(") ");
      printf("}");
      if (j<3) printf(",\n    ");
    }
    printf("\n  }");
    if (i<3) printf(",\n  ");
  }
  printf("\n};\n\n");

  printf("const u128 F1_inverse[4][4] =\n{\n  ");
  for (i=0; i<4; i++)
  {
    printf("{\n    ");
    for (j=0; j<4; j++)
    {
      printf("{ ");
      printf("UINT64_C(0x");
      bigint_qr(F1_inverse[i][j], r, F1_inverse[i][j], base);
      bigint_print(r, 16);
      printf("), ");
      printf("UINT64_C(0x");
      bigint_qr(F1_inverse[i][j], r, F1_inverse[i][j], base);
      bigint_print(r, 16);
      printf(") ");
      printf("}");
      if (j<3) printf(",\n    ");
    }
    printf("\n  }");
    if (i<3) printf(",\n  ");
  }
  printf("\n};\n\n");

  printf("const u192 F2[2][2] =\n{\n  ");
  for (i=0; i<2; i++)
  {
    printf("{\n    ");
    for (j=0; j<2; j++)
    {
      printf("{ ");
      printf("UINT64_C(0x");
      bigint_qr(F2[i][j], r, F2[i][j], base);
      bigint_print(r, 16);
      printf("), ");
      printf("UINT64_C(0x");
      bigint_qr(F2[i][j], r, F2[i][j], base);
      bigint_print(r, 16);
      printf("), ");
      printf("UINT64_C(0x");
      bigint_qr(F2[i][j], r, F2[i][j], base);
      bigint_print(r, 16);
      printf(") ");
      printf("}");
      if (j<1) printf(",\n    ");
    }
    printf("\n  }");
    if (i<1) printf(",\n  ");
  }
  printf("\n};\n\n");

  printf("const u192 F2_inverse[2][2] =\n{\n  ");
  for (i=0; i<2; i++)
  {
    printf("{\n    ");
    for (j=0; j<2; j++)
    {
      printf("{ ");
      printf("UINT64_C(0x");
      bigint_qr(F2_inverse[i][j], r, F2_inverse[i][j], base);
      bigint_print(r, 16);
      printf("), ");
      printf("UINT64_C(0x");
      bigint_qr(F2_inverse[i][j], r, F2_inverse[i][j], base);
      bigint_print(r, 16);
      printf("), ");
      printf("UINT64_C(0x");
      bigint_qr(F2_inverse[i][j], r, F2_inverse[i][j], base);
      bigint_print(r, 16);
      printf(") ");
      printf("}");
      if (j<1) printf(",\n    ");
    }
    printf("\n  }");
    if (i<1) printf(",\n  ");
  }
  printf("\n};\n\n");

  return 0;
}

