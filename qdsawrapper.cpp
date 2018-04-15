//============================================================================
// Name        : qdsa_wrapper.cpp
// Description : Faulty qDSA signature generator based on the C reference implementation
// Standard    : C++11
//============================================================================
#include <iostream>
#include <gmpxx.h>
#include "mocksig.h"
#include "qdsawrapper.h"

mpz_class qdsa::keygen(
        unsigned char *pk,
        unsigned char *sk
        )
{
    /*
     * Generate a 64-byte pseudo-random string and a public
     * key x-coordinate on the Montgomery curve.
     * The secret key is clamped before usage.
     *
     * Input:
     *      sk (64 bytes): 32 bytes of randomness (in lower bytes)
     *
     * Output:
     *      pk (32 bytes): Public key, x-coordinate (no sign bit)
     *      sk (64 bytes): Pseudo-random secret
     *      dgmp (mpz_class): Secret key d'
     */

    ecp R;
    fe25519 rx;
    group_scalar d;

    hash(sk, sk, 32);
    sk[32] &= 248;
    sk[63] &= 127;
    sk[63] |= 64;

    group_scalar_get32(&d, sk+32);
    ladder_base(&R, &d);

    compress(&rx, &R);
    fe25519_pack(pk, &rx);

    mpz_class dgmp;
    mpz_import(dgmp.get_mpz_t(), 32, -1, sizeof(sk[32]), 0, 0, sk+32);

    return dgmp;
}

int qdsa::sign(
		SignatureSimple& sig,
        unsigned char *sm, unsigned long long *smlen,
        const unsigned char *m, unsigned long long mlen,
        const unsigned char *pk, const unsigned char *sk,
        mpz_class n, sc25519 *lim
        )
{
    /*
     * Generate a signature consisting of a 32-byte
     * x-coordinate on the Montgomery curve and an
     * integer modulo the curve order.
     * Append the message to the signature.
     *
     * Input:
     *      m: Message
     *      mlen: Message length in bytes
     *      pk (32 bytes): Public key, x-coordinate (no sign bit)
     *      sk (64 bytes): Pseudo-random secret
     *      n (mpz_class): Order of the base point
     *
     * Output:
     *      sm (64+mlen bytes): Signature + Message
     *      smlen: 64+mlen
     */

	static const group_scalar inv = { .v = {
        0xc20dca53c5b85ef2, 0xfa73b66fa39b5a0, 0, 0xc00000000000000 } };

	// uncomment the following to preprocess 3-bit biased signatures
    //static const group_scalar inv = { .v = {
    //    0x6106E529E2DC2F79, 0x7D39DB37D1CDAD0, 0, 0x600000000000000 } };


    unsigned long long i;
    ecp R;
    fe25519 rx;
    group_scalar r, h, s;

    *smlen = mlen+64;
    for(i=0;i<mlen;i++) { sm[64+i] = m[i]; }
    for(i=0;i<32;i++) { sm[32+i] = sk[i]; } // set d'
    hash(sm, sm+32, mlen+32); // compute r
    group_scalar_get64(&r, sm); // set r

    uint8_t rbytes[32], rr;
    group_scalar_pack(rbytes, &r);
    rr = rbytes[0] & 3;
    if (rr != 0 && rr != 2) return 0; // keep only r mod 4 = 0 or 2;
    //rr = rbytes[0] & 7;
    //if (rr != 0 && rr != 4) return 0; // keep only r mod 8 = 0 or 4;


    ladder_base(&R, &r);
    ladder_fault(&R, &r);
    compress(&rx, &R);

    for(i=0;i<32;i++) { sm[32+i] = pk[i]; }
    fe25519_pack(sm, &rx);
    hash(sm, sm, mlen+64);
    group_scalar_get64(&h, sm);
    group_scalar_get32(&s, sk+32);

    group_scalar_set_pos(&h);

	// preprocess and filter h here
    group_scalar hprime = { .v = { 0, 0, 0, 0 } };
    group_scalar_mul(&hprime, &h, &inv);
    if (sc25519_lt(&hprime,lim)==0) return 0;

    group_scalar_mul(&s, &h, &s);
    group_scalar_sub(&s, &r, &s);

    // Not needed here:
    fe25519_pack(sm, &rx);
    group_scalar_pack(sm+32, &s);

    /* Preprocess signature using leaked nonce
     *
     * Let
     *  h' := h/2^b mod n
     *  s' := (s - (r mod 2^b))/2^b mod n
     *  r' := (r - (r mod 2^b))/2^b mod n
     * Then
     *  r' = s' + h'd mod n
     */
    group_scalar sprime = { .v = { rr, 0, 0, 0 } };
    group_scalar_sub(&sprime, &s, &sprime);
    group_scalar_mul(&sprime, &sprime, &inv);
    
    mpz_class hgmp;
    mpz_import(hgmp.get_mpz_t(), 4, -1, sizeof(hprime.v[0]), 0, 0, hprime.v);

    mpz_class sgmp;
    mpz_import(sgmp.get_mpz_t(), 4, -1, sizeof(sprime.v[0]), 0, 0, sprime.v);

    sig = SignatureSimple(hgmp, sgmp);
    return 1;
}

int qdsa::verify(
        unsigned char *m, long long mlen,
        unsigned char *sm, unsigned long long smlen,
        const unsigned char *pk
        )
{
    /*
     * Verify correctness of a signature with respect
     * to a public key. Return 1 if correct, 0 if
     * incorrect, and return the message.
     *
     * Input:
     *      sm (64+mlen bytes): Signature + Message
     *      smlen: 64+mlen
     *      pk (32 bytes): Public key, x-coordinate (no sign bit)
     *
     * Output:
     *      0 if correct, 1 if incorrect
     *      m: Message
     *      mlen: Message length (bytes)
     */

    unsigned long long i;
    ecp sP, hQ;
    fe25519 rx, pkx;
    fe25519 bZZ, bXZ, bXX;
    group_scalar s, h;

    fe25519_unpack(&rx, sm);
    group_scalar_get32(&s, sm+32);
    for(i=0;i<smlen-64;i++) { m[i] = sm[64+i]; }
    mlen = smlen-64;

    for(i=0;i<32;i++) { sm[32+i] = pk[i]; }
    hash(sm, sm, mlen+64);
    group_scalar_get64(&h, sm);

    fe25519_unpack(&pkx, pk);
    decompress(&sP, &pkx);
    ladder(&hQ, &sP, &pkx, &h);
    ladder_base(&sP, &s);

    bValues(&bZZ, &bXZ, &bXX, &sP, &hQ);
    return check(&bZZ, &bXZ, &bXX, &rx);
}
