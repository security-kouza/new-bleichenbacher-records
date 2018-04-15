#ifndef SIGN_H
#define SIGN_H

#include "hash.h"
#include "ladder.h"
#include "scalar.h"

int keypair(
        unsigned char *pk, 
        unsigned char *sk
        );
int sign(
        unsigned char *sm, unsigned long long *smlen,
        const unsigned char *m, unsigned long long mlen,
        const unsigned char *pk, const unsigned char *sk
        );
int verify(
        unsigned char *m, long long mlen,
        unsigned char *sm, unsigned long long smlen,
        const unsigned char *pk
        );

#endif
