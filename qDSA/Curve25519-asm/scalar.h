#ifndef GROUP_SCALAR_H
#define GROUP_SCALAR_H

#define GROUP_SCALAR_PACKEDBYTES 32
#include "fe25519.h"
#include "sc25519.h"

#define group_scalar		sc25519
#define group_scalar_get32	sc25519_from32bytes
#define group_scalar_get64	sc25519_from64bytes
#define group_scalar_pack	sc25519_to32bytes
#define group_scalar_mul	sc25519_mul
#define group_scalar_lt		sc25519_lt
#define group_scalar_short  shortsc25519
//#define gs					sc25519
#define gs_get32			sc25519_from32bytes
#define gs_get64			sc25519_from64bytes
#define gs_pack				sc25519_to32bytes
#define gs_mul				sc25519_mul
#define gs_lt				sc25519_lt
//#define gss                 shortsc25519

/* additional functions */

#define group_scalar_sub	sc25519_sub
#define group_scalar_set_pos	sc25519_abs
#define gs_sub					sc25519_sub
#define gs_add					sc25519_add
#define gs_set_pos				sc25519_abs
#define gs_sub_nored			sc25519_sub_nored

void sc25519_sub(sc25519 *r, const sc25519 *x, const sc25519 *y);
void sc25519_abs(sc25519 *r);

#endif
