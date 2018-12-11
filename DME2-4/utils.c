#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "utils.h"

void error_msg(const char *err)
{
  fprintf(stderr, "error: %s!\n", err);
  exit(-1);
}

void print_hex(const unsigned char *hex, int n)
{
  int i;
  for (i=0; i<n; i++)
    printf("%02x", hex[i]);
}

void decode_hex(unsigned char *hex, const char *str, int n)
{
  int i, j, k, c, h;
  for (i=j=0; i<n; i++)
  {
    for (k=h=0; k<2; k++,j++)
    {
      c = tolower(str[j]);
      if (!c)
        error_msg("string too short");
      if (!isxdigit(c))
        error_msg("invalid hex digit");
      h <<= 4;
      h += (c >= '0' && c <= '9')  ? c-'0' : c-'a'+10;
    }
    hex[i] = h;
  }
  if (str[j]) error_msg("string too long");
}

