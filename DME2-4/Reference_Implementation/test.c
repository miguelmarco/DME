#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "dme.h"
#include "utils.h"

int main(int argc, char *argv[])
{
  int i, n, num_tests;
  unsigned char pt1[PT_BYTES], pt2[PT_BYTES];
  unsigned char ct1[CT_BYTES], ct2[CT_BYTES];
  unsigned char skey[SKEY_BYTES], pkey[PKEY_BYTES];
  const unsigned char pad[PADDING] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  unsigned char pad_ret[PADDING];
  clock_t t1, t2;
  double s2p_interp, enc, dec, ews;
  double s2p_interp_tot, enc_tot, dec_tot, ews_tot;
  double s2p_interp_sqr, enc_sqr, dec_sqr, ews_sqr;
  FILE *stats;

  num_tests = 10;
  if (argc == 2)
    num_tests = atoi(argv[1]);

  s2p_interp_tot = enc_tot = dec_tot = ews_tot = 0.0;
  s2p_interp_sqr = enc_sqr = dec_sqr = ews_sqr = 0.0;
  printf("--- tests ---\n");
  for (n=0; n<num_tests; n++)
  {
    printf("test #%d\n", n+1);

    /* Create random secret key */
    for (i=0; i<SKEY_BYTES; i++)
      skey[i] = rand() & 0xff;

    /* Convert secret to public key */
    t1 = clock();
    if (dme_skey_to_pkey(pkey, skey))
      continue;
    t2 = clock();
    s2p_interp_tot += s2p_interp = t2-t1;
    s2p_interp_sqr += s2p_interp*s2p_interp;

    printf("secret key = ");
    print_hex(skey, SKEY_BYTES);
    printf("\n");

    printf("public key = ");
    print_hex(pkey, PKEY_BYTES);
    printf("\n");

    /* Create random plaintext */
    for (i=0; i<PT_BYTES; i++)
      pt1[i] = rand() & 0xff;
    printf("plaintext  = ");
    print_hex(pt1, PT_BYTES);
    printf("\n");

    /* Encrypt with secret key */
    t1 = clock();
    dme_encrypt_with_skey(ct1, pt1, skey, pad);
    t2 = clock();
    ews_tot += ews = t2-t1;
    ews_sqr += ews*ews;
    printf("ciphertext = ");
    print_hex(ct1, CT_BYTES);
    printf("\n");

    /* Encrypt with public key */
    t1 = clock();
    if (dme_encrypt(ct2, pt1, pkey, pad))
      error_msg("encrypt()");
    t2 = clock();
    enc_tot += enc = t2-t1;
    enc_sqr += enc*enc;

    /* Make sure both ciphertexts coincide */
    for (i=0; i<CT_BYTES; i++)
      if (ct1[i] != ct2[i])
        error_msg("ciphertexts don't match\n");

    /* Decrypt */
    t1 = clock();
    if (dme_decrypt(pt2, ct2, skey, pad_ret))
      error_msg("decrypt()");
    if (memcmp(pad, pad_ret, PADDING))
      error_msg("decrypt()");
    t2 = clock();
    dec_tot += dec = t2-t1;
    dec_sqr += dec*dec;

    /* Make sure both plaintexts coincide */
    for (i=0; i<PT_BYTES; i++)
      if (pt1[i] != pt2[i])
        error_msg("plaintexts don't match\n");

    printf("skey_to_pkey: %.3f [usec]\n",
      s2p_interp*1e6/CLOCKS_PER_SEC);
    printf("encrypt: %.3f [usec]\n", enc*1e6/CLOCKS_PER_SEC);
    printf("decrypt: %.3f [usec]\n", dec*1e6/CLOCKS_PER_SEC);
    printf("encrypt_with_skey: %.3f [usec]\n", ews*1e6/CLOCKS_PER_SEC);

    printf("\n");
  }

  stats = fopen("stats.txt", "w");
  fprintf(stats, "--- statistics ---\n");
  fprintf(stats,
    "skey_to_pkey: mean = %.3f [usec], std_dev = %.3f [usec]\n",
    s2p_interp_tot*1e6/CLOCKS_PER_SEC/n,
    sqrt(s2p_interp_sqr*n-s2p_interp_tot*s2p_interp_tot)*1e6/CLOCKS_PER_SEC/sqrt(n*(n-1)));
  fprintf(stats,
    "encrypt: mean = %.3f [usec], std_dev = %.3f [usec]\n",
    enc_tot*1e6/CLOCKS_PER_SEC/n,
    sqrt(enc_sqr*n-enc_tot*enc_tot)*1e6/CLOCKS_PER_SEC/sqrt(n*(n-1)));
  fprintf(stats,
    "decrypt: mean = %.3f [usec], std_dev = %.3f [usec]\n",
    dec_tot*1e6/CLOCKS_PER_SEC/n,
    sqrt(dec_sqr*n-dec_tot*dec_tot)*1e6/CLOCKS_PER_SEC/sqrt(n*(n-1)));
  fprintf(stats,
    "encrypt_with_skey: mean = %.3f [usec], std_dev = %.3f [usec]\n",
    ews_tot*1e6/CLOCKS_PER_SEC/n,
    sqrt(ews_sqr*n-ews_tot*ews_tot)*1e6/CLOCKS_PER_SEC/sqrt(n*(n-1)));
  fclose(stats);

  return 0;
}
