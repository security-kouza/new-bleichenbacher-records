#ifndef HASH_H
#define HASH_H

#define cKeccakB    1600
#define cKeccakR    1024
#define cKeccakR_SizeInBytes    (1344 / 8)
#define crypto_hash_BYTES   64

int hash( unsigned char *out, const unsigned char *in, unsigned long long inlen );

#endif
