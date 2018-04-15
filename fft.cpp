//============================================================================
// Name        : fft.cpp
// Description : Bias Computation with FFT
// Standard    : C++11
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <iostream>
#include <vector>

#include <fftw3.h>
#include <gmpxx.h>

#include "mocksig.h"
#include "fft.h"

double find_peak(const double a[], uint32_t L) {
	uint32_t i;
	double peak = a[0];
	uint32_t peak_at = 0;

	for(i = 0; i < L; i++) {
		if (a[i] > peak) {
			peak = a[i];
			peak_at = i;
		}
	}
	printf("find_peak: peak found at w = %u, peak = %lf \n", peak_at, peak);
    return peak_at;
}

double* compute_bias(uint32_t L, uint32_t kb, Domain pp, std::vector<SignatureSimple>& sigs) {
	uint32_t j;
    /* init Z_t */
	fftw_plan p;
    fftw_complex* Z = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L);
    fftw_complex* W = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L);
    double* W_norm = (double *) malloc(sizeof(double) * L);

    p = fftw_plan_dft_1d(L, Z, W, FFTW_BACKWARD, FFTW_ESTIMATE);
    for (j = 0; j < L; j++) {
        Z[j][0] = 0;
        Z[j][1] = 0;
    }
    printf("preparing FFT table...\n");
    mpf_class tmp;
    for (j = 0; j < L; j++) {
    	tmp = (mpf_class)2*sigs[j].s/pp.n;
    	mpz_class idx = sigs[j].h >> kb;
    	Z[idx.get_ui()][0] += cos(M_PI*tmp.get_d());
    	Z[idx.get_ui()][1] += sin(M_PI*tmp.get_d());
    }

    /* doit */
    printf("computing bias...\n");
    fftw_execute(p);

    /* compute norm */
    for(j = 0; j < L; j++){
        W_norm[j] = sqrt(pow(W[j][0], 2.0) + pow(W[j][1], 2.0))/L;
        //if (W_norm[j] > 0.1)
          //printf("Bn(%u) = %lf = |%+lf %+lf*i|\n", j, W_norm[j], W[j][0], W[j][1]);
    }

    fftw_destroy_plan(p);
    fftw_free(Z); fftw_free(W);

    return W_norm;
}

double* compute_bias(uint32_t L, uint32_t kb, Domain pp, std::vector<SignatureSC25519>& sigs) {
	uint32_t j;
    /* init Z_t */
	fftw_plan p;
    fftw_complex* Z = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L);
    fftw_complex* W = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L);
    double* W_norm = (double *) malloc(sizeof(double) * L);

    p = fftw_plan_dft_1d(L, Z, W, FFTW_BACKWARD, FFTW_ESTIMATE);
    for (j = 0; j < L; j++) {
        Z[j][0] = 0;
        Z[j][1] = 0;
    }
    printf("preparing FFT table...\n");
    mpf_class tmp;
    for (j = 0; j < L; j++) {
    	mpz_class h,s;
    	gs_to_mpz(sigs[j].h, h);
    	gs_to_mpz(sigs[j].s, s);
    	tmp = (mpf_class)2*s/pp.n;
    	mpz_class idx = h >> kb;
    	Z[idx.get_ui()][0] += cos(M_PI*tmp.get_d());
    	Z[idx.get_ui()][1] += sin(M_PI*tmp.get_d());
    }

    /* doit */
    printf("computing bias...\n");
    fftw_execute(p);

    /* compute norm */
    for(j = 0; j < L; j++){
        W_norm[j] = sqrt(pow(W[j][0], 2.0) + pow(W[j][1], 2.0))/L;
        if (W_norm[j] > 0.1)
          printf("Bn(%u) = %lf = |%+lf %+lf*i|\n", j, W_norm[j], W[j][0], W[j][1]);
    }

    fftw_destroy_plan(p);
    fftw_free(Z); fftw_free(W);

    return W_norm;
}

double compute_noise_avg(int L, int peak_at, const double W[]) {
	int j;
	double sum = 0;
    for (j = 0; j < L; j++) {
    	if (j != peak_at) sum += W[j];
    }
    return sum/(L-1);
}

