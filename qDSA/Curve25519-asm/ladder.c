#include "fe25519.h"
#include "ladder.h"

inline void cswap(fe25519 *x, fe25519 *y, int b)
{
    fe25519_cswap(x, y, b);
}

void ecswap(ecp *x, ecp *y, int b)
{
    cswap(&x->X,&y->X,b);
    cswap(&x->Z,&y->Z,b);
}

void xDBLADD(ecp *xp, ecp *xq, const fe25519 *xd)
{
    /*
     * Simultaneous xDBL and xADD operation on
     * the Montgomery curve.
     *
     * Input: 
     *      xp: proj. x-coordinate on Montgomery curve
     *      xq: proj. x-coordinate on Montgomery curve
     *      xd: affine x-coordinate of difference xp-xq
     *
     * Output: 
     *      xp: proj. x-coordinate of 2*xp
     *      xq: proj. x-coordinate of xp+xq
     */

    fe25519 b0,b1;

    fe25519_add(&b0, &xp->X, &xp->Z);
    fe25519_sub(&b1, &xp->X, &xp->Z);
    fe25519_add(&xp->X, &xq->X, &xq->Z);
    fe25519_sub(&xp->Z, &xq->X, &xq->Z);
    fe25519_mul(&xq->X, &xp->Z, &b0);
    fe25519_mul(&xp->Z, &xp->X, &b1);
    fe25519_add(&xp->X, &xp->Z, &xq->X);
    fe25519_sub(&xq->Z, &xq->X, &xp->Z);
    fe25519_square(&xq->X, &xp->X);
    fe25519_square(&xp->X, &xq->Z);
    fe25519_mul(&xq->Z, &xp->X, xd);
    fe25519_square(&xp->X, &b0);
    fe25519_square(&b0, &b1);
    fe25519_sub(&xp->Z, &xp->X, &b0);
    fe25519_mul(&xp->X, &b0, &xp->X);
    fe25519_mul121666(&b1, &xp->Z);
    fe25519_add(&b1, &b1, &b0);
    fe25519_mul(&xp->Z, &b1, &xp->Z);
}

void ladder(ecp *r, ecp *xp, const fe25519 *xpw, const group_scalar *n)
{
    /*
     * Montgomery ladder computing n*xp via repeated
     * differential additions and constant-time
     * conditional swaps.
     *
     * Input: 
     *      xp: proj. x-coordinate on Montgomery curve
     *      xpw: affine x-coordinate of xp
     *      n: Scalar (max 255-bit)
     *
     * Output: 
     *      xr: proj. x-coordinate of n*xq
     *      xp: proj. x-coordinate of (n+1)*xp
     */

    int i,swap; int bit = 0; int prevbit = 0;
    unsigned char nbits[32];

    fe25519_setone(&r->X);
    fe25519_setzero(&r->Z);
    group_scalar_pack(nbits, n);
    for(i=254; i>=0; i--)
    {
        bit = ( nbits[i>>3] >> ( i & 7 ) ) & 1;
        swap = bit ^ prevbit;
        prevbit = bit;

        ecswap(r, xp, swap);
        xDBLADD(r, xp, xpw);
    }
    ecswap(r, xp, bit);
}

void ladder_base(ecp *r, const group_scalar *n)
{
    ecp base;
    fe25519 basex;

    fe25519_setzero(&base.X);
    fe25519_setone(&base.Z);
    base.X.v[0] = 9;
    fe25519_copy(&basex, &base.X);

    ladder(r, &base, &basex, n);
}

/* BEGIN: MODIFIED BY AUTHORS */
void ladder_fault(ecp *r, const group_scalar *n)
{
    ecp base;
    fe25519 basex;

    fe25519_setone(&base.X);
    fe25519_setone(&base.Z);
    fe25519_copy(&basex, &base.X);

    ladder(r, &base, &basex, n);
}
/* END: MODIFIED BY AUTHORS */

void bValues(fe25519 *bZZ, fe25519 *bXZ, fe25519 *bXX, 
        const ecp *xp, const ecp *xq)
{
    /*
     * Three biquadratic forms B_XX, B_XZ and B_ZZ
     * in the coordinates of xp and xq  
     *
     * Input: 
     *      xp: proj. x-coordinate on Montgomery curve
     *      xq: proj. x-coordinate on Montgomery curve
     *
     * Output: 
     *      bZZ: Element B_ZZ of fe25519
     *      bXZ: Element B_XZ of fe25519
     *      bXX: Element B_XX of fe25519
     */

    fe25519 b0,b1,b2;

    fe25519_mul(&b0, &xp->X, &xq->X);
    fe25519_mul(&b1, &xp->Z, &xq->Z);
    fe25519_sub(bZZ, &b0, &b1);
    fe25519_square(bZZ, bZZ);
    fe25519_add(&b0, &b0, &b1);

    fe25519_mul(&b1, &xp->X, &xq->Z);
    fe25519_mul(&b2, &xq->X, &xp->Z);
    fe25519_sub(bXX, &b1, &b2);
    fe25519_square(bXX, bXX);

    fe25519_add(bXZ, &b1, &b2);
    fe25519_mul(bXZ, bXZ, &b0);
    fe25519_mul(&b0, &b1, &b2);
    fe25519_add(&b0, &b0, &b0);
    fe25519_add(&b0, &b0, &b0);
    fe25519_add(&b1, &b0, &b0);
    fe25519_mul121666(&b1, &b1);
    fe25519_sub(&b0, &b1, &b0);
    fe25519_add(bXZ, bXZ, &b0);
    fe25519_add(bXZ, bXZ, bXZ);
}

int check(fe25519 *bZZ, fe25519 *bXZ, fe25519 *bXX,
        const fe25519 *rx)
{
    /*
     * Verify whether B_XXrx^2 - B_XZrx + B_ZZ = 0
     *
     * Input: 
     *      bZZ: Biquadratic form B_ZZ
     *      bXZ: Biquadratic form B_XZ
     *      bXX: Biquadratic form B_XX
     *      rx: affine x-coordinate on Montgomery curve
     *
     * Output: 
     *      1 if B_XXrx^2 - B_XZrx + B_ZZ = 0,
     *      0 otherwise
     */

    fe25519 b0, b1;

    fe25519_square(&b0, rx);
    fe25519_mul(&b0, &b0, bXX);
    fe25519_mul(&b1, rx, bXZ);
    fe25519_sub(&b0, &b0, &b1);
    fe25519_add(&b0, &b0, bZZ);

    return fe25519_iszero(&b0);
}

void compress(fe25519 *r, const ecp *xp)
{
    /*
     * Compress from projective representation (X : Z)
     * to affine x = X*Z^{p-2}, where p = 2^255-19
     *
     * Input: 
     *      xp: proj. x-coordinate (X : Z)
     *
     * Output: 
     *      r: affine x-coordinate x = X*Z^{p-2}
     */

    fe25519_invert(r, &xp->Z);  
    fe25519_mul(r, &xp->X, r);
    fe25519_freeze(r);
}

void decompress(ecp *xp, const fe25519 *r)
{
    /*
     * Decompress from affine representation x
     * to projective (x : 1)
     *
     * Input: 
     *      r: affine x-coordinate x
     *
     * Output: 
     *      xp: proj. x-coordinate (x : 1)
     */

    fe25519_copy(&xp->X, r);
    fe25519_setone(&xp->Z);
}
