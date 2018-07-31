#ifndef MOCKSIG_H
#define MOCKSIG_H

extern "C" {

#include "qDSA/Curve25519-asm/sc25519.h"
#include "qDSA/Curve25519-asm/scalar.h"

}
#include <vector>

struct Domain {
	uint32_t secpar;
	mpz_class n;
};

class SignatureLeak {
public:
    mpz_class h;
    mpz_class s;
    mpz_class rr; // leaked nonce

    SignatureLeak(mpz_class h, mpz_class s, mpz_class rr) : h(h), s(s), rr(rr) {}

    bool operator < (const SignatureLeak& sig) const {
    	return (h < sig.h);
    }

    bool operator > (const SignatureLeak& sig) const {
    	return (h > sig.h);
    }
};

class SignatureSimple {
public:
    mpz_class h;
    mpz_class s;

    SignatureSimple() : h(0), s(0) {}
    SignatureSimple(mpz_class h, mpz_class s) : h(h), s(s) {}

    bool operator < (const SignatureSimple& sig) const {
    	return (h < sig.h);
    }

    bool operator > (const SignatureSimple& sig) const {
    	return (h > sig.h);
    }
};

class SignatureSC25519 {
public:
    sc25519 h;
    sc25519 s;

    SignatureSC25519() : h({.v={0,0,0,0}}), s({.v={0,0,0,0}}) {}
    SignatureSC25519(sc25519 h, sc25519 s) : h(h), s(s) {}

    bool operator < (const SignatureSC25519& sig) const {
    	return (1==sc25519_lt(&h, &sig.h));
    }

    bool operator > (const SignatureSC25519& sig) const {
    	return (0==sc25519_lt(&h, &sig.h));
    }
};

typedef struct {
    uint64_t v[4];
} gs;

typedef struct {
    uint64_t v[2];
} gss;

namespace mock {
Domain setup(uint32_t secpar, mpz_class n);
mpz_class keygen(Domain pp, gmp_randclass& rand);
SignatureLeak sign(Domain pp, mpz_class d, int M, uint32_t leak, gmp_randclass& rand);
SignatureSimple sign_filter(Domain pp, mpz_class d, int M, uint32_t leak, uint32_t filter, gmp_randclass& rand);
}

void sigprint(SignatureSimple sig);
void sigprint(SignatureSC25519 sig);
void sigsave(std::vector<SignatureSimple>& sigs, std::string filename);
void sigsave(std::vector<SignatureSC25519>& sigs, std::string filename);
template <class Iterator>
void sigsave_it(Iterator start, Iterator end, std::string filename) {
    FILE *fp = fopen(filename.c_str(), "wb");

    for (auto it = start; it != end; ++it) {
        mpz_out_raw(fp, it->h.get_mpz_t());
        mpz_out_raw(fp, it->s.get_mpz_t());
    }
    fclose(fp);
}
void sigload(std::vector<SignatureSimple>& sigs, std::string filename);
void sigload(std::vector<SignatureSC25519>& sigs, std::string filename);
std::vector<uint8_t> mpz_to_vector(const mpz_class x);
mpz_class vector_to_mpz(std::vector<uint8_t>& v);
void mpz_to_gs(gs& sc, mpz_class& x);
void mpz_to_gs(sc25519& sc, mpz_class& x);
void mpz_to_gss(gss& sc, mpz_class& x, const int& offset);
void gs_to_mpz(gs& sc, mpz_class& x);
void gs_to_mpz(sc25519& sc, mpz_class& x);
void pack(std::vector<gs>& hlist, std::vector<SignatureSimple>* sigsptr);
void packhalf(std::vector<gss>& hlist, std::vector<SignatureSimple>* sigsptr, const int& offset);
template <class Iterator>
void packhalf_it(Iterator start, Iterator end, std::vector<gss>& hlist, const int& size, const int& offset) {
    hlist.reserve(size);
	for (auto it = start; it != end; ++it) {
	    gss half_h;
		mpz_to_gss(half_h, it->h, offset);
		hlist.emplace_back(half_h);
	}
}
void unpack(std::vector<gs>& hlist, std::vector<SignatureSimple>* sigsptr);
void packedprint(std::vector<gs>& hlist);
void gsprint(gs& sc);
void gsprint(sc25519& sc);
void gss_add(gss *r, const gss *x, const gss *y);
void gss_sub(gss *r, const gss *x, const gss *y);
int gss_lt(const gss *a, const gss *b);
int gss_lteq(const gss *a, const gss *b);
#endif
