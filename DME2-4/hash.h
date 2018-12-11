#include <stdint.h>

uint32_t lrot32(uint32_t x, uint32_t c);
uint32_t rrot32(uint32_t x, uint32_t c);
uint32_t get32le(char *x);
uint32_t get32be(char *x);
void put32le(char *x, uint32_t r);
void put32be(char *x, uint32_t r);
uint64_t lrot64(uint64_t x, uint64_t c);
uint64_t rrot64(uint64_t x, uint64_t c);
uint64_t get64le(char *x);
uint64_t get64be(char *x);
void put64le(char *x, uint64_t r);
void put64be(char *x, uint64_t r);

struct md4_ctx
{
  uint32_t h[4];
  char buf[64];
  size_t pending;
  uint64_t length;
};

void md4_init_ctx(struct md4_ctx *ctx);
void md4_process_bytes(char *buf, size_t len, struct md4_ctx *ctx);
void md4_finish_ctx(struct md4_ctx *ctx, char *result);

struct md5_ctx
{
  uint32_t h[4];
  char buf[64];
  size_t pending;
  uint64_t length;
};

void md5_init_ctx(struct md5_ctx *ctx);
void md5_process_bytes(char *buf, size_t len, struct md5_ctx *ctx);
void md5_finish_ctx(struct md5_ctx *ctx, char *result);

struct sha1_ctx
{
  uint32_t h[5];
  char buf[64];
  size_t pending;
  uint64_t length;
};

void sha1_init_ctx(struct sha1_ctx *ctx);
void sha1_process_bytes(char *buf, size_t len, struct sha1_ctx *ctx);
void sha1_finish_ctx(struct sha1_ctx *ctx, char *result);

struct sha2_224_ctx
{
  uint32_t h[8];
  char buf[64];
  size_t pending;
  uint64_t length;
};

void sha2_224_init_ctx(struct sha2_224_ctx *ctx);
void sha2_224_process_bytes(char *buf, size_t len, struct sha2_224_ctx *ctx);
void sha2_224_finish_ctx(struct sha2_224_ctx *ctx, char *result);

struct sha2_256_ctx
{
  uint32_t h[8];
  char buf[64];
  size_t pending;
  uint64_t length;
};

void sha2_256_init_ctx(struct sha2_256_ctx *ctx);
void sha2_256_process_bytes(char *buf, size_t len, struct sha2_256_ctx *ctx);
void sha2_256_finish_ctx(struct sha2_256_ctx *ctx, char *result);

struct sha2_384_ctx
{
  uint64_t h[8];
  char buf[128];
  size_t pending;
  uint64_t length[2];
};

void sha2_384_init_ctx(struct sha2_384_ctx *ctx);
void sha2_384_process_bytes(char *buf, size_t len, struct sha2_384_ctx *ctx);
void sha2_384_finish_ctx(struct sha2_384_ctx *ctx, char *result);

struct sha2_512_ctx
{
  uint64_t h[8];
  char buf[128];
  size_t pending;
  uint64_t length[2];
};

void sha2_512_init_ctx(struct sha2_512_ctx *ctx);
void sha2_512_process_bytes(char *buf, size_t len, struct sha2_512_ctx *ctx);
void sha2_512_finish_ctx(struct sha2_512_ctx *ctx, char *result);

