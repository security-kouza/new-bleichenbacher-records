//============================================================================
// Name        : mocksig.cpp
// Description : Mock Schnorr signature generator and various utilities
// Standard    : C++11
//============================================================================
#include <inttypes.h>
#include <stdio.h>
#include <iostream>
#include <vector>

#include <gmpxx.h>

#include "mocksig.h"

#define MPZ_LOAD_ERROR 0
#define MPZ_BYTE_SIZE 34

Domain mock::setup(uint32_t secpar, mpz_class n) {
    Domain pp = {
    		secpar,
			n,
    };

    return pp;
}

mpz_class mock::keygen(Domain pp, gmp_randclass& rand) {
    mpz_class d = rand.get_z_range(pp.n);
    return d;
}

SignatureLeak mock::sign(Domain pp, mpz_class d, int M, uint32_t leak, gmp_randclass& rand) {
	mpz_class r = rand.get_z_range(pp.n);
	mpz_class rr;
	mpz_mod_ui(rr.get_mpz_t(), r.get_mpz_t(), 1<<leak);

    mpz_class h = rand.get_z_range(pp.n);
    mpz_class s;
    mpz_class tmp = r - h*d;
	mpz_mod(s.get_mpz_t(), tmp.get_mpz_t(), pp.n.get_mpz_t());

	SignatureLeak sigma(h, s, rr);

    return sigma;
}

SignatureSimple mock::sign_filter(Domain pp, mpz_class d, int M, uint32_t leak, uint32_t filter, gmp_randclass& rand) {
	mpz_class nonce_lim = pp.n >> leak;
	mpz_class hlim = pp.n >> filter;
	mpz_class r,h,s,stmp;
	r = rand.get_z_range(nonce_lim);
	h = rand.get_z_range(hlim);
	stmp = r - h*d;
	mpz_mod(s.get_mpz_t(), stmp.get_mpz_t(), pp.n.get_mpz_t());

	SignatureSimple sigma(h,s);

    return sigma;
}

/* utilities */
void sigprint(SignatureSimple sig) {
	std::cout << "(h=" << sig.h << ", s=" << sig.s << ")" << std::endl;
}

void sigprint(SignatureSC25519 sig) {
	mpz_class h,s;
	gs_to_mpz(sig.h, h);
	gs_to_mpz(sig.s, s);
	gmp_printf("(h=%Zd, s=%Zd)\n",h.get_mpz_t(),s.get_mpz_t());
}

void sigsave(std::vector<SignatureSimple>& sigs, std::string filename) {
    FILE *fp = fopen(filename.c_str(), "wb");

    for (auto& sig: sigs) {
        mpz_out_raw(fp, sig.h.get_mpz_t());
        mpz_out_raw(fp, sig.s.get_mpz_t());
    }
    fclose(fp);
}

void sigsave(std::vector<SignatureSC25519>& sigs, std::string filename) {
    FILE *fp = fopen(filename.c_str(), "wb");

    mpz_class h,s;
    for (auto& sig: sigs) {
        gs_to_mpz(sig.h, h);
        gs_to_mpz(sig.s, s);
        mpz_out_raw(fp, h.get_mpz_t());
        mpz_out_raw(fp, s.get_mpz_t());
    }
    fclose(fp);
}

void sigload(std::vector<SignatureSimple>& sigs, std::string filename) {
    FILE *fp = fopen(filename.c_str(), "rb");
    mpz_class h, s;

    while (1) {
        if (mpz_inp_raw(h.get_mpz_t(), fp) == MPZ_LOAD_ERROR) break;
        if (mpz_inp_raw(s.get_mpz_t(), fp) == MPZ_LOAD_ERROR) break;
        sigs.emplace_back(h, s);
        //sigprint(SignatureSimple(h, s));
    }
    fclose(fp);
};

void sigload(std::vector<SignatureSC25519>& sigs, std::string filename) {
    FILE *fp = fopen(filename.c_str(), "rb");
    mpz_class h, s;

    while (1) {
        if (mpz_inp_raw(h.get_mpz_t(), fp) == MPZ_LOAD_ERROR) break;
        if (mpz_inp_raw(s.get_mpz_t(), fp) == MPZ_LOAD_ERROR) break;
        sc25519 temph, temps;
        mpz_to_gs(temph, h);
        mpz_to_gs(temps, s);

        sigs.emplace_back(temph, temps);
        //sigprint(SignatureSimple(h, s));
    }
    fclose(fp);
};

std::vector<uint8_t> mpz_to_vector(const mpz_class x) {
    size_t size = MPZ_BYTE_SIZE;
    std::vector<uint8_t> v(size);
    mpz_export(&v[0], &size, -1, 1, 0, 0, x.get_mpz_t());
    if (sgn(x) < 0) v.back() = 1; // last byte indicates the sign

    return v;
}

mpz_class vector_to_mpz(std::vector<uint8_t>& v) {
    mpz_class x;
    bool is_negative = false;
    if (v.back() == 1) {
        is_negative = true;
        v.back() = 0;
    }
    mpz_import(x.get_mpz_t(), v.size(), -1, 1, 0, 0, &v[0]);

    if (is_negative) x = -x;
    return x;
}

void mpz_to_gs(gs& sc, mpz_class& x) {
	for (auto& vv : sc.v) vv = 0;
	mpz_export(sc.v, NULL, -1, 8, 0, 0, x.get_mpz_t());
}

void mpz_to_gs(sc25519& sc, mpz_class& x) {
	sc = sc25519({.v = {0,0,0,0}});
	mpz_export(sc.v, NULL, -1, 8, 0, 0, x.get_mpz_t());
}

void mpz_to_gss(gss& sc, mpz_class& x, const int& offset) {
    gs temp;
    mpz_to_gs(temp, x);
	sc.v[0] = temp.v[offset];
	sc.v[1] = temp.v[offset+1];
}

void gs_to_mpz(gs& sc, mpz_class& x) {
	mpz_import(x.get_mpz_t(), 4, -1 , 8, 0, 0, sc.v);
}

void gs_to_mpz(sc25519& sc, mpz_class& x) {
	mpz_import(x.get_mpz_t(), 4, -1 , 8, 0, 0, sc.v);
}

void pack(std::vector<gs>& hlist, std::vector<SignatureSimple>* sigsptr) {
    hlist.reserve(sigsptr->size());
	for (auto& sig: *sigsptr) {
	    gs hscalar;
		mpz_to_gs(hscalar, sig.h);
		hlist.emplace_back(hscalar);
	}
}

void packhalf(std::vector<gss>& hlist, std::vector<SignatureSimple>* sigsptr, const int& offset) {
    hlist.reserve(sigsptr->size());
	for (auto& sig: *sigsptr) {
	    gss half_h;
		mpz_to_gss(half_h, sig.h, offset);
		hlist.emplace_back(half_h);
	}
}

void unpack(std::vector<gs>& hlist, std::vector<SignatureSimple>* sigsptr) {
	mpz_class x;
	for (auto& h: hlist) {
		gs_to_mpz(h, x);
		sigsptr->emplace_back(SignatureSimple(x,0));
	}
}

void packedprint(std::vector<gs>& hlist) {
	for (auto& h: hlist) {
		gsprint(h);
	}
}

void gsprint(gs& sc) {
	mpz_class x;
	gs_to_mpz(sc, x);
	gmp_printf("%Zd\n", x.get_mpz_t());
}

void gsprint(sc25519& sc) {
	mpz_class x;
	gs_to_mpz(sc, x);
	gmp_printf("%Zd\n", x.get_mpz_t());
}

/* arithmetic for short group scalar (=gss) */
void gss_add(gss *r, const gss *x, const gss *y) {
    r->v[1] = x->v[1] + y->v[1];
    if ( (UINT64_MAX - x->v[0]) < y->v[0]) {
        r->v[1] = r->v[1] + 1;
        r->v[0] = y->v[0] - (UINT64_MAX - x->v[0]) - 1;
    } else {
        r->v[0] = x->v[0] + y->v[0];
    }
}

// assume x > y
void gss_sub(gss *r, const gss *x, const gss *y) {
    r->v[1] = x->v[1] - y->v[1];
    if (x->v[0] < y->v[0]) {
        r->v[1] = r->v[1] - 1;
        r->v[0] = UINT64_MAX - (y->v[0] - x->v[0] - 1);
    } else {
        r->v[0] = x->v[0] - y->v[0];
    }
}

int gss_lt(const gss *a, const gss *b) {
    if (a->v[1] < b->v[1]) return 1;
    else if (a->v[1] > b->v[1]) return 0;
    else if (a->v[0] < b->v[0]) return 1;
    else if (a->v[0] > b->v[0]) return 0;
    return 0;
}

int gss_lteq(const gss *a, const gss *b) {
    if (a->v[0] == b->v[0] && a->v[1] == b->v[1]) return 1;
    return gss_lt(a, b);
}
