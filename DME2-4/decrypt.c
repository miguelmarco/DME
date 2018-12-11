#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dme.h"
#include "utils.h"

int main(int argc, char *argv[])
{
  unsigned char ct[CT_BYTES], pt[PT_BYTES], skey[SKEY_BYTES];
  unsigned char pad[PADDING] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  unsigned char pad_ret[PADDING];

  if (argc < 3 || argc > 4)
  {
    fprintf(stderr, "use: %s cipher_text secret_key [padding]\n", argv[0]);
    exit(-1);
  }

  decode_hex(ct, argv[1], CT_BYTES);
  decode_hex(skey, argv[2], SKEY_BYTES);
  if (argc == 4 && strcmp(argv[3], "-"))
    decode_hex(pad, argv[3], PADDING);
  if (dme_decrypt(pt, ct, skey, pad_ret))
    error_msg("decrypt()");
  if (argc == 4 && memcmp(pad, pad_ret, PADDING))
    error_msg("padding mismatch");
  print_hex(pt, PT_BYTES);
  printf("\n");

  return 0;
}

