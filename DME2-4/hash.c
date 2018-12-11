#include <string.h>
#include <stdint.h>
#include "hash.h"

/* utils */
uint32_t lrot32(uint32_t x, uint32_t c)
{
  return (x << c) | (x >> (32-c));
}

uint32_t rrot32(uint32_t x, uint32_t c)
{
  return (x >> c) | (x << (32-c));
}

uint32_t get32be(char *x)
{
  int i;
  uint32_t r = 0;
  for (i=0; i<4; i++)
  {
    r <<= 8;
    r |= (unsigned char)x[i];
  }
  return r;
}

void put32be(char *x, uint32_t r)
{
  int i;
  for (i=0; i<4; i++)
  {
    x[3-i] = r & 0xff;
    r >>= 8;
  }
}

void put64be(char *x, uint64_t r)
{
  int i;
  for (i=0; i<8; i++)
  {
    x[7-i] = r & 0xff;
    r >>= 8;
  }
}

/* sha2_256 */
static const uint32_t sha2_256_k[64] = {
   0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
   0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
   0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
   0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
   0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
   0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
   0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
   0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

void sha2_256_block(char *buf, struct sha2_256_ctx *ctx) /* buf: 64 bytes */
{
  uint32_t a, b, c, d, e, f, g, h, i, s0, s1, t1, t2, maj, ch, w[64];
  a = ctx->h[0];
  b = ctx->h[1];
  c = ctx->h[2];
  d = ctx->h[3];
  e = ctx->h[4];
  f = ctx->h[5];
  g = ctx->h[6];
  h = ctx->h[7];
  for (i=0; i<16; i++)
    w[i] = get32be(&buf[i<<2]);
  for (i=16; i<64; i++)
  {
    s0 = rrot32 (w[i-15],7) ^ rrot32 (w[i-15], 18) ^ (w[i-15] >> 3);
    s1 =  rrot32(w[i-2], 17) ^ rrot32(w[i-2], 19) ^ (w[i-2] >> 10);
    w[i] = w[i-16] + s0 + w[i-7] + s1;
  }
  for (i=0; i<64; i++)
  {
    s0 =  rrot32(a, 2) ^ rrot32(a, 13) ^ rrot32(a, 22);
    maj = (a & b) ^ (a & c) ^ (b & c);
    t2 = s0 + maj;
    s1 =  rrot32(e, 6) ^ rrot32(e, 11) ^ rrot32(e, 25);
    ch = (e & f) ^ ((~e) & g);
    t1 = h + s1 + ch + sha2_256_k[i] + w[i];
    h = g;
    g = f;
    f = e;
    e = d + t1;
    d = c;
    c = b;
    b = a;
    a = t1 + t2;
  }
  ctx->h[0] = ctx->h[0] + a;
  ctx->h[1] = ctx->h[1] + b;
  ctx->h[2] = ctx->h[2] + c;
  ctx->h[3] = ctx->h[3] + d;
  ctx->h[4] = ctx->h[4] + e;
  ctx->h[5] = ctx->h[5] + f;
  ctx->h[6] = ctx->h[6] + g;
  ctx->h[7] = ctx->h[7] + h;
}

void sha2_256_init_ctx(struct sha2_256_ctx *ctx)
{
  ctx->h[0] = 0x6a09e667;
  ctx->h[1] = 0xbb67ae85;
  ctx->h[2] = 0x3c6ef372;
  ctx->h[3] = 0xa54ff53a;
  ctx->h[4] = 0x510e527f;
  ctx->h[5] = 0x9b05688c;
  ctx->h[6] = 0x1f83d9ab;
  ctx->h[7] = 0x5be0cd19;
  ctx->pending = ctx->length  = 0;
}

void sha2_256_process_bytes(char *buf, size_t len, struct sha2_256_ctx *ctx)
{
  ctx->length += len;
  if (ctx->pending)
  {
    if (len < 64-ctx->pending)
    {
      memcpy(ctx->buf + ctx->pending, buf, len);
      ctx->pending += len;
      len = 0;
    }
    else
    {
      memcpy(ctx->buf + ctx->pending, buf, 64-ctx->pending);
      sha2_256_block(ctx->buf, ctx);
      buf += 64-ctx->pending;
      len -= 64-ctx->pending;
      ctx->pending = 0;
    }
  }
  while (len >= 64)
  {
    sha2_256_block(buf, ctx);
    len -= 64;
    buf += 64;
  }
  if (len)
  {
    ctx->pending = len;
    memcpy(ctx->buf, buf, len);
  }
}

void sha2_256_finish_ctx(struct sha2_256_ctx *ctx, char *result) /* result: 32 bytes */
{
  int i;
  ctx->buf[ctx->pending++] = (char) 0x80;
  if (ctx->pending > 56)
  {
    memset(ctx->buf + ctx->pending, 0, 64-ctx->pending);
    sha2_256_block(ctx->buf, ctx);
    ctx->pending = 0;
  }
  memset(ctx->buf + ctx->pending, 0, 56-ctx->pending);
  put64be(&ctx->buf[56], ctx->length<<3);
  sha2_256_block(ctx->buf, ctx);
  for (i=0; i<8; i++)
    put32be(&result[i<<2], ctx->h[i]);
}
