#ifndef QDSAWRAPPER_H
#define QDSAWRAPPER_H

extern "C" {

#include "qDSA/Curve25519-asm/hash.h"
#include "qDSA/Curve25519-asm/ladder.h"
#include "qDSA/Curve25519-asm/scalar.h"

}

class FaultException  : public std::exception {
   const char * what () const throw () {
      return "r mod n = ?1";
   }
};

namespace qdsa{
mpz_class keygen(
        unsigned char *pk,
        unsigned char *sk
        );
int sign(
		SignatureSimple& sig,
        unsigned char *sm, unsigned long long *smlen,
        const unsigned char *m, unsigned long long mlen,
        const unsigned char *pk, const unsigned char *sk,
        mpz_class n, sc25519 *lim
        );
int verify(
        unsigned char *m, long long mlen,
        unsigned char *sm, unsigned long long smlen,
        const unsigned char *pk
        );
}
#endif
