#ifndef API_H
#define API_H

#include "dme.h"

#define CRYPTO_SECRETKEYBYTES SKEY_BYTES
#define CRYPTO_PUBLICKEYBYTES PKEY_BYTES
#define CRYPTO_BYTES PT_BYTES
#define CRYPTO_CIPHERTEXTBYTES CT_BYTES

#define CRYPTO_ALGNAME "DME-KEM N=2, M=4, E=48, S=16"

int crypto_kem_keypair(unsigned char *pk, unsigned char *sk);
int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk);
int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *pk);

#endif

