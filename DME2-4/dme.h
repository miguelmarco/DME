#ifndef DME_H
#define DME_H

#include <stdint.h>

/* Basic functions */
int max(int a, int b);
int min(int a, int b);

/* Arithmetic in Fq */
typedef uint_least64_t fq_elem;

fq_elem fq_add(fq_elem a, fq_elem b);
fq_elem fq_mul(fq_elem a, fq_elem b);
fq_elem fq_inv(fq_elem a);
fq_elem fq_rnd(void);
fq_elem fq_pow(fq_elem a, uint_least64_t n);
fq_elem fq_pow_2exp(fq_elem a, unsigned int n);

void fq_poly_multiply(fq_elem *a, const fq_elem *b, const fq_elem *c, int deg_b,
  int deg_c);
void fq_matrix_multiply(fq_elem *a, const fq_elem *b, const fq_elem *c, int n,
  int m, int l);
int fq_matrix_inverse(fq_elem *a, const fq_elem *b, int n);

/* Arithmetic in Fq2 */
typedef fq_elem fq2_elem[2];
typedef uint_least64_t u128[2];

void fq2_add(fq2_elem a, const fq2_elem b, const fq2_elem c);
void fq2_mul(fq2_elem a, const fq2_elem b, const fq2_elem c);
void fq2_pow(fq2_elem a, const fq2_elem b, const u128 n);
int fq2_is0(const fq2_elem a);
void fq2_pow_2exp(fq2_elem a, fq2_elem b, unsigned int n);

/* Arithmetic in Fq4 */
typedef fq_elem fq4_elem[4];
typedef uint_least64_t u192[3];

void fq4_add(fq4_elem a, const fq4_elem b, const fq4_elem c);
void fq4_mul(fq4_elem a, const fq4_elem b, const fq4_elem c);
void fq4_pow(fq4_elem a, const fq4_elem b, const u192 n);
int fq4_is0(const fq4_elem a);
void fq4_pow_2exp(fq4_elem a, fq4_elem b, unsigned int n);

/* Key structures */
struct secret_key
{
  fq_elem L1[4][2][2], L1_inverse[4][2][2];
  fq_elem L2[2][4][4], L2_inverse[2][4][4];
  fq_elem L3[2][4][4], L3_inverse[2][4][4];
};

typedef struct secret_key secret_key;

struct public_key
{
  fq_elem coeffs1[4][64];
  fq_elem coeffs2[4][64];
};

typedef struct public_key public_key;

/* API (internal) */
int fq_skey_to_pkey(public_key *pkey, const secret_key *skey);
int fq_encrypt(fq_elem *ct, const fq_elem *pt, const public_key *pkey);
int fq_decrypt(fq_elem *pt, const fq_elem *ct, const secret_key *skey);
int fq_encrypt_with_skey(fq_elem *ct, const fq_elem *pt, const secret_key *skey);

/* Serialization */
#define PKEY_BYTES (6*8*64)
#define SKEY_BYTES (6*(4*2*2+2*4*4+2*4*4))
#define TEXT_BYTES (6*8)
#define PADDING (2*8)
#define PT_BYTES (TEXT_BYTES-PADDING)
#define CT_BYTES TEXT_BYTES
void serialize_skey(unsigned char *p, const secret_key *skey);
void serialize_pkey(unsigned char *p, const public_key *pkey);
void serialize_text(unsigned char *p, const fq_elem *text);
int parse_skey(secret_key *skey, const unsigned char *p);
int parse_pkey(public_key *pkey, const unsigned char *p);
void parse_text(fq_elem *text, const unsigned char *p);

/* API (external) */
int dme_encrypt(unsigned char *ct, const unsigned char *pt, const unsigned char
  *pkey, const unsigned char *padding);
int dme_decrypt(unsigned char *pt, const unsigned char *ct, const unsigned char
  *skey, unsigned char *padding);
int dme_encrypt_with_skey(unsigned char *ct, const unsigned char *pt, const
  unsigned char *skey, const unsigned char *padding);
int dme_skey_to_pkey(unsigned char *pkey, const unsigned char *skey);

/* API (NIST) */
#ifdef NIST
#include "rng.h"
#include "api.h"
#endif

#endif

