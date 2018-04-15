#ifndef FFT_H
#define FFT_H

double find_peak(const double a[], uint32_t L);

double* compute_bias(uint32_t L, uint32_t kb, Domain pp, std::vector<SignatureSimple>& sigs);
double* compute_bias(uint32_t L, uint32_t kb, Domain pp, std::vector<SignatureSC25519>& sigs);

double compute_noise_avg(int L, int peak_at, const double W[]);

#endif
