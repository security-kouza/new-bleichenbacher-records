#include "dh.h"

int dh_keygen(
        unsigned char *pk,
        unsigned char *sk
        )
{
    /*
     * Generate a public key x-coordinate on the Montgomery curve.
     * The secret key is clamped before usage.
     *
     * Input: 
     *      sk (32 bytes): 32 bytes of randomness
     *
     * Output:
     *      pk (32 bytes): Public key, x-coordinate (no sign bit)
     */

    ecp R;
    fe25519 rx;
    group_scalar d;

    sk[0] &= 248;
    sk[31] &= 127;
    sk[31] |= 64;
    group_scalar_get32(&d, sk);
    ladder_base(&R, &d);

    compress(&rx, &R);
    fe25519_pack(pk, &rx);

    return 0;
}

int dh_exchange(
        unsigned char *ss,
        unsigned char *pk,
        unsigned char *sk
        )
{
    /*
     * Generate a shared secret
     *
     * Input: 
     *      sk (32 bytes): 32 bytes of randomness
     *      pk (32 bytes): Public key, x-coordinate (no sign bit)
     *
     * Output:
     *      ss (32 bytes): Shared secret, x-coordinate (no sign bit)
     */

    ecp SS, R;
    fe25519 rx;
    group_scalar d;

    fe25519_unpack(&rx, pk);
    decompress(&R, &rx);

    group_scalar_get32(&d, sk);
    ladder(&SS, &R, &rx, &d);

    compress(&rx, &SS);
    fe25519_pack(ss, &rx);

    return 0;
}
