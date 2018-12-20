#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dme.h"
#include "utils.h"

int main(int argc, char *argv[])
{
  unsigned char ct[CT_BYTES], pt[PT_BYTES], pkey[PKEY_BYTES];
  unsigned char pad[PADDING] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

  if (argc < 3 || argc > 4)
  {
    fprintf(stderr, "use: %s plain_text public_key [padding]\n", argv[0]);
    exit(-1);
  }

  decode_hex(pt, argv[1], PT_BYTES);
  decode_hex(pkey, argv[2], PKEY_BYTES);
  if (argc == 4)
    decode_hex(pad, argv[3], PADDING);
  if (dme_encrypt(ct, pt, pkey, pad))
    error_msg("encrypt()");
  print_hex(ct, CT_BYTES);
  printf("\n");

  return 0;
}
