#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dme.h"
#include "utils.h"

int main(int argc, char *argv[])
{
  unsigned char pkey[PKEY_BYTES], skey[SKEY_BYTES];

  if (argc != 2)
  {
    fprintf(stderr, "use: %s secret_key\n", argv[0]);
    exit(-1);
  }

  decode_hex(skey, argv[1], SKEY_BYTES);
  if (dme_skey_to_pkey(pkey, skey))
    error_msg("skey_to_pkey()");
  print_hex(pkey, PKEY_BYTES);
  printf("\n");

  return 0;
}

