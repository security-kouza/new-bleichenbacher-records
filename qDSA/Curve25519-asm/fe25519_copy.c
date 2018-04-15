#include "fe25519.h"

void fe25519_copy(fe25519 *r, const fe25519 *x)
{
    r->v[0] = x->v[0];
    r->v[1] = x->v[1];
    r->v[2] = x->v[2];
    r->v[3] = x->v[3];
}
void fe25519_cswap(fe25519 *r, fe25519 *x, unsigned char b) {
    unsigned long long mask = -b, t;

#define cswap(a,b) do { t=mask&((a)^(b)); (a)^=t; (b)^=t; } while(0)
    cswap(x->v[0], r->v[0]);
    cswap(x->v[1], r->v[1]);
    cswap(x->v[2], r->v[2]);
    cswap(x->v[3], r->v[3]);
#undef cswap
}

void fe25519_mul121666(fe25519 *r, const fe25519 *x)
{
    const fe25519 c = { .v = {121666,0,0,0} };
    fe25519_mul(r, x, &c);
}
