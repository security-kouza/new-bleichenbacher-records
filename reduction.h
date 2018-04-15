#ifndef REDUCTION_H
#define REDUCTION_H

class LRComb {

public:
	mpz_class hsum;
	uint32_t idx_L;
	uint32_t idx_R;
	LRComb (mpz_class, uint32_t, uint32_t);
    bool operator < (const LRComb& c) const {
    	return (hsum > c.hsum);
    }
    bool operator > (const LRComb& c) const {
    	return (hsum < c.hsum);
    }
};

class LRCombGSS {

public:
	gss hsum;
	uint32_t idx_L;
	uint32_t idx_R;
	LRCombGSS () : hsum({.v={0,0}}), idx_L(0), idx_R(0) {}
	LRCombGSS (gss hh, uint32_t i, uint32_t j) : hsum(hh), idx_L(i), idx_R(j) {}
    bool operator < (const LRCombGSS& c) const {
    	return (gss_lt(&hsum, &(c.hsum)) == 1);
    }
    bool operator > (const LRCombGSS& c) const {
    	return (gss_lt(&hsum, &(c.hsum)) == 0);
    }
};

class Index {

public:
	uint32_t idx_L1;
	uint32_t idx_R1;
	uint32_t idx_L2;
	uint32_t idx_R2;
	bool flip;
	Index ();
	Index (uint32_t, uint32_t, uint32_t, uint32_t, bool);
};

void idxsave(std::vector<Index>& is, const std::string& filename);
void idxload(std::vector<Index>& is, const std::string& filename);
void restore_from_idx(std::vector<SignatureSimple>& sigs, const std::vector<Index>& idxlist,
			 const std::vector<SignatureSimple>& L1, const std::vector<SignatureSimple>& R1,
			 const std::vector<SignatureSimple>& L2, const std::vector<SignatureSimple>& R2,
			 const mpz_class& threshold_mpz);

void schroeppel_shamir(std::vector<SignatureSimple> * sigs,
		const uint32_t secpar, const uint32_t l, const uint32_t b, const uint32_t filter, const int iota);

void schroeppel_shamir_mpi(std::vector<SignatureSimple>& sigs,
		const uint32_t secpar, const uint32_t l, const double b,
		const uint32_t filter, const uint32_t a, const int log_prec,
		const std::vector<int>& ofst_info,
		const int iota, const std::string out_prefix, const bool istest);

void sort_and_difference(std::vector<SignatureSimple>& sigs,
		const uint32_t secpar, const uint32_t l, const uint32_t b,
		const uint32_t filter, const uint32_t a, const int log_prec,
		const std::vector<int>& ofst_info,
		const int iota, const std::string out_prefix, const bool istest);

void sort_and_difference(std::vector<SignatureSC25519>& sigs,
		const uint32_t secpar, const uint32_t l, const uint32_t b,
		const uint32_t filter, const uint32_t a, const int log_prec,
		const std::vector<int>& ofst_info,
		const int iota, const std::string out_prefix, const bool istest);
#endif
