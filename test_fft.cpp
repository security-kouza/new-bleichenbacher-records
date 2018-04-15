#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <iostream>
#include <vector>

#include <gmpxx.h>

#include "mocksig.h"
#include "fft.h"

int main(void) {
	/* Parameters */
    const uint32_t secpar = 22;
    const mpz_class k = 3;
    const uint32_t leak = 2;
    const int M = 0; // dummy message

    /* Z/nZ */
    mpz_class n;
    mpz_ui_pow_ui(n.get_mpz_t(), 2, secpar);
    n = n - k;
    const uint32_t L = n.get_ui(); // generate as many signatures as the group order

    /* SetUp & KeyGen */
    gmp_randclass rand(gmp_randinit_default);
	rand.seed(12345678);
    Domain pp = mock::setup(secpar, n);
    mpz_class d = mock::keygen(pp, rand);
    gmp_printf("pp = (secpar=%u, n=%Zd), d=%Zd\n", pp.secpar, pp.n.get_mpz_t(), d.get_mpz_t());

    /* Signature Generation */
    printf("generating %u signatures...\n", L);
	mpz_class div, inv, h_tmp, s_tmp;
	div = mpz_class(1) << leak;
	mpz_invert(inv.get_mpz_t(), div.get_mpz_t(), pp.n.get_mpz_t());
    std::vector<SignatureSimple>* sigsptr = new std::vector<SignatureSimple>();
    while (sigsptr->size() < L) {
    	SignatureLeak sig = mock::sign(pp, d, M, leak, rand);
    	sig.h = sig.h*inv;
    	mpz_mod(sig.h.get_mpz_t(), sig.h.get_mpz_t(), pp.n.get_mpz_t());
    	sig.s = (sig.s - sig.rr)*inv;
    	mpz_mod(sig.s.get_mpz_t(), sig.s.get_mpz_t(), pp.n.get_mpz_t());
    	sigsptr->push_back(SignatureSimple(sig.h, sig.s));
    	if (sigsptr->size() % (L/10) == 0) printf("%.0f %% done\n", (float)sigsptr->size()*100/L);
    }

    /* FFT */
    double* W = compute_bias(L, 0, pp, *sigsptr);
    uint32_t peak_at = find_peak(W, L);
    double noise = compute_noise_avg(L, peak_at, W);
    mpf_class estim = (mpf_class)peak_at*n/L;

    gmp_printf("estimated sk w = %.5Ff \n"
    		"correct sk d = %Zd \n", estim.get_mpf_t(), d.get_mpz_t());
    printf("average noise = %lf\n"
    		"correct noise 1/sqrt(L) = %lf\n", noise, 1/sqrt(L));

	return EXIT_SUCCESS;

}
