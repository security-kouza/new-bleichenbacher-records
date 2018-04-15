//============================================================================
// Name        : reduction.cpp
// Description : Reduction phase in Bleichenbacher's attack
// Standard    : C++11
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdint.h>
#include <omp.h>
#include <malloc.h>
#include <iostream>
#include <algorithm>
#include <queue>
#include <bitset>

#include <gmpxx.h>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>

#include "mocksig.h"
#include "reduction.h"

namespace mpi = boost::mpi;

LRComb::LRComb (mpz_class hh, uint32_t i, uint32_t j) {
	hsum = hh;
	idx_L = i;
	idx_R = j;
}

// default constructor
Index::Index () {
	idx_L1 = 0;
	idx_R1 = 0;
	idx_L2 = 0;
	idx_R2 = 0;
	flip = false;
}
Index::Index (uint32_t i1, uint32_t j1, uint32_t i2, uint32_t j2, bool f=false) {
	idx_L1 = i1;
	idx_R1 = j1;
	idx_L2 = i2;
	idx_R2 = j2;
	flip = f;
}

/* serialize user-defined types */
namespace boost {
	namespace serialization {
		template<class Archive>
		void serialize(Archive &ar, gs &scalar, const unsigned int version) {
			ar & scalar.v;
		}

		template<class Archive>
		void serialize(Archive &ar, gss &scalar, const unsigned int version) {
			ar & scalar.v;
		}

		template<class Archive>
		void serialize(Archive &ar, Index &x, const unsigned int version) {
			ar & x.idx_L1;
            ar & x.idx_R1;
            ar & x.idx_L2;
            ar & x.idx_R2;
            ar & x.flip;
		}
	}
}

void idxsave(std::vector<Index>& is, const std::string& filename) {
    FILE *fp = fopen(filename.c_str(), "wb");

    for (auto& i: is) {
    	fwrite(&i.idx_L1, sizeof(uint32_t), 1, fp);
    	fwrite(&i.idx_R1, sizeof(uint32_t), 1, fp);
    	fwrite(&i.idx_L2, sizeof(uint32_t), 1, fp);
    	fwrite(&i.idx_R2, sizeof(uint32_t), 1, fp);
    	fwrite(&i.flip, sizeof(bool), 1, fp);
    }
    fclose(fp);
}

void idxload(std::vector<Index>& is, const std::string& filename) {
    FILE *fp = fopen(filename.c_str(), "rb");
    Index i;
    while (fread(&i.idx_L1, sizeof(uint32_t), 1, fp) &&
    	   fread(&i.idx_R1, sizeof(uint32_t), 1, fp) &&
		   fread(&i.idx_L2, sizeof(uint32_t), 1, fp) &&
		   fread(&i.idx_R2, sizeof(uint32_t), 1, fp) &&
		   fread(&i.flip, sizeof(bool), 1, fp)) {
    	is.emplace_back(i);
    }
    fclose(fp);
}

std::string msb(mpz_class m, uint32_t a, uint32_t bit) {
	std::string bin = m.get_str(2);
	uint32_t padding = bit - bin.length();
	if (padding != 0) bin = std::string(padding, '0') + bin;

	return bin.substr(0, a);
}

/* Contents of sigsptr will be replaced with reduced values */
void schroeppel_shamir(std::vector<SignatureSimple> * sigsptr,
		const uint32_t secpar, const uint32_t l, const uint32_t b, const uint32_t filter, const int iota=1) {
    uint32_t threshold_bit = secpar - filter;
    uint32_t S = sigsptr->size();
    uint32_t keep = 1<<l;
    std::vector<SignatureSimple>* L1;
    std::vector<SignatureSimple>* R1;
    std::vector<SignatureSimple>* L2;
    std::vector<SignatureSimple>* R2;

    printf("schroeppel_shamir: got %u sigs\n", S);

    for (int round = 0; round < iota; round++) {
        S = sigsptr->size();

    	/* split sigs into L1 || R1 || L2 || R2 */
    	uint32_t q1 = S/4;
    	uint32_t q2 = S/2;
    	uint32_t q3 = S*3/4;

        std::cout << "Splitting into 4 lists..." << std::endl;
        L1 = new std::vector<SignatureSimple>(sigsptr->begin(), sigsptr->begin() + q1);
        R1 = new std::vector<SignatureSimple>(sigsptr->begin() + q1, sigsptr->begin() + q2);
        L2 = new std::vector<SignatureSimple>(sigsptr->begin() + q2, sigsptr->begin() + q3);
        R2 = new std::vector<SignatureSimple>(sigsptr->begin() + q3, sigsptr->end());
        delete sigsptr;
        sigsptr = new std::vector<SignatureSimple>();

        threshold_bit -= b;

    	mpz_class threshold;
    	mpz_ui_pow_ui(threshold.get_mpz_t(), 2, threshold_bit);

        /* Init heaps */
    	std::cout << "Sorting R1 and R2..." << std::endl;
        std::sort(R1->begin(), R1->end());
        std::sort(R2->begin(), R2->end());

        std::cout << "Initializing heaps..." << std::endl;
        std::priority_queue<LRComb, std::vector<LRComb>> heap1;
        std::priority_queue<LRComb, std::vector<LRComb>> heap2;
        for (uint32_t i = 0; i != L1->size(); i++) {
        	heap1.push(LRComb((*L1)[i].h + (*R1)[0].h, i, 0));
        	heap2.push(LRComb((*L2)[i].h + (*R2)[0].h, i, 0));
        }

        /* Start looking for collisions */
        std::cout << "Trying to get 2^" << l <<
        		" values less than 2^" << threshold_bit << std::endl;
        uint32_t counter = 0;
        uint32_t num_col = 0;
        mpz_class hsum1, hsum2, h_diff, s_diff;
        uint32_t i1, j1, i2, j2;
        while (num_col < keep && heap1.empty() == false && heap2.empty() == false) {
        	counter++;
        	hsum1 = heap1.top().hsum;
        	i1 = heap1.top().idx_L;
        	j1 = heap1.top().idx_R;

        	hsum2 = heap2.top().hsum;
        	i2 = heap2.top().idx_L;
        	j2 = heap2.top().idx_R;

        	h_diff = hsum1 - hsum2;
        	s_diff = (*L1)[i1].s + (*R1)[j1].s - (*L2)[i2].s - (*R2)[j2].s;

        	if (h_diff > 0) {
        		if (j2 < S/4 - 1) {
        			heap2.pop();
        			heap2.push(LRComb((*L2)[i2].h + (*R2)[j2+1].h, i2, j2 + 1));
        		} else {
        			heap2.pop();
        		}
        	} else {
        		h_diff = -h_diff;
        		s_diff = -s_diff;
        		if (j1 < S/4 -1) {
        			heap1.pop();
        			heap1.push(LRComb((*L1)[i1].h + (*R1)[j1+1].h, i1, j1 + 1));
        		} else {
        			heap1.pop();
        		}
        	}

        	if (h_diff < threshold) {
        	    sigsptr->emplace_back(SignatureSimple(h_diff, s_diff));
        		num_col++;
	        	if (num_col % 100 == 0) printf("%u/%u collisions found, %.2lf %% done \n", num_col, counter, num_col*100.0/keep);
        	}
        }
        if (num_col < keep) puts("WARNING: failed to get expected amount of reduced values!");
        std::cout << "completed after " << counter << " loops" << std::endl;
        puts("-------------------------");
        delete L1;
        delete R1;
        delete L2;
        delete R2;
    }
}

/* helper method that collects linear combinations of two whose (a+1)-MSB is equal to A */
void collect_lc_two_GSS(std::vector<LRCombGSS>& combs, const std::vector<gss>& L, const std::vector<gss>& R,
						const uint32_t& A, const uint32_t& a, const double& current_threshold_bit_f, const int& pckofst, const uint32_t& lim) {
	uint32_t i=0,j=0;
	bool flag_i = false;
	bool flag_j = false;
	uint32_t bad = 0;

	mpz_class A0tmp, A1tmp, Amtmp, tmp;
	gss A0, A1, Amid;
    uint32_t current_threshold_bit_z = current_threshold_bit_f; // rounded
    double extra = pow(2.0, current_threshold_bit_f - current_threshold_bit_z);
    tmp = mpz_class((mpf_class(1) << (current_threshold_bit_z-a))*mpf_class(extra));
    A0tmp = A*tmp;
    A1tmp = (A+1)*tmp-1;
    Amtmp = A*tmp + tmp/2;

    //std::cout << A0tmp << std::endl;
    //std::cout << A1tmp << std::endl;
    //std::cout << Amtmp << std::endl;

	mpz_to_gss(A0, A0tmp, pckofst);
	mpz_to_gss(A1, A1tmp, pckofst);
	mpz_to_gss(Amid, Amtmp, pckofst);

    uint32_t lsize = L.size();
    uint32_t rsize = R.size();

    gss comb;
    gss_add(&comb, &L[i], &R[j]);
	while (combs.size() < lim) {
		if (gss_lt(&comb, &A0)) {
			if (j == rsize-1) break;
			j++;
		}
		else if (gss_lt(&A1, &comb)) {
			if (i == lsize-1) break;
			i++;
		}
		else {
			combs.emplace_back(comb, i, j);
			/* check if indices are at the end */
			if (i == lsize-1 && j == rsize-1) break;
			else if (i == lsize-1) {
			    j++;
			    gss_add(&comb, &L[i], &R[j]);
			    continue;
			}
			else if (j == rsize-1) {
			    i++;
			    gss_add(&comb, &L[i], &R[j]);
			    continue;
			}
			gss peek_i, peek_j;
			gss_add(&peek_i, &L[i+1], &R[j]);
			gss_add(&peek_j, &L[i], &R[j+1]);
			flag_i = (gss_lteq(&A0, &peek_i) && gss_lteq(&peek_i, &A1));
			flag_j = (gss_lteq(&A0, &peek_j) && gss_lteq(&peek_j, &A1));
			if (flag_i ^ flag_j) {
			    if (flag_i) i++;
			    else j++;
			}
			else {
			    if (flag_i && flag_j) bad++;
			    uint32_t c = 1;
                if (gss_lt(&comb, &Amid)) {
                     while (flag_i && i+c < lsize && combs.size() < lim) {
                         gss_add(&peek_i, &L[i+c], &R[j]);
                         flag_i = (gss_lteq(&A0, &peek_i) && gss_lteq(&peek_i, &A1));
                         if (flag_i) combs.emplace_back(peek_i, i+c, j);
                         c++;
                     }
                     j++;
                }
                else {
                     while (flag_j && j+c < rsize && combs.size() < lim) {
                         gss_add(&peek_j, &L[i], &R[j+c]);
                         flag_j = (gss_lteq(&A0, &peek_j) && gss_lteq(&peek_j, &A1));
                         if (flag_j) combs.emplace_back(peek_j, i, j+c);
                         c++;
                     }
                     i++;
                }
			}
		}
		gss_add(&comb, &L[i], &R[j]);
	}
	//std::cout << "bad: " << bad << std::endl;
	//std::cout << i << "," << j << std::endl;
	//printf("A=%u: got %lu collisions\n", A, combs->size());
}

void sort_and_difference(std::vector<SignatureSC25519>& sigs,
		const uint32_t secpar, const uint32_t l, const uint32_t b,
		const uint32_t filter, const uint32_t a, const int log_prec, const std::vector<int>& ofst_info, const int iota=1, const std::string out_prefix="", const bool istest=true){
    uint32_t threshold_bit = secpar - filter;
    uint32_t S;

    for (int round = 0; round < iota; round++){
        S = sigs.size();
        threshold_bit -= b;
    	mpz_class threshold_mpz = mpz_class(1) << threshold_bit;
    	sc25519 threshold_sc;
    	mpz_to_gs(threshold_sc, threshold_mpz);

        std::sort(sigs.begin(), sigs.end());

        sc25519 newh, news;
        uint32_t j = 0;
        for(uint32_t i=0; i < S-1; i++) {
        	sc25519_sub(&newh, &(sigs[i+1].h), &(sigs[i].h));
        	if (sc25519_lt(&newh, &threshold_sc)) {
        		sc25519_sub(&news, &(sigs[i+1].s), &(sigs[i].s));
        		sigs[j] = SignatureSC25519(newh,news);
        		j++;
        	}
        }
        sigs.resize(j);
        printf("sorting done; %lu elements with h < 2^%u obtained\n",sigs.size(),threshold_bit);
        if (out_prefix.length()) {
            // file format: redsigs_sd_round-i.bin
            std::string outsig = out_prefix + "_round-" + std::to_string(round) + ".bin";
            printf("saving signatures to %s... \n", outsig.c_str());
            sigsave(sigs, outsig);
        }
    }
    return;
}

void sort_and_difference(std::vector<SignatureSimple>& sigs,
		const uint32_t secpar, const uint32_t l, const uint32_t b,
		const uint32_t filter, const uint32_t a, const int log_prec, const std::vector<int>& ofst_info, const int iota=1, const std::string out_prefix="", const bool istest=true){
    uint32_t threshold_bit = secpar - filter;
    uint32_t S;

    for (int round = 0; round < iota; round++){
        S = sigs.size();
        threshold_bit -= b;
    	mpz_class threshold_mpz = mpz_class(1) << threshold_bit;
        std::sort(sigs.begin(), sigs.end());

        mpz_class newh, news;
        uint32_t j = 0;
        for(uint32_t i=0; i < S-1; i++) {
        	newh = sigs[i+1].h - sigs[i].h;
        	if (newh < threshold_mpz) {
        		news = sigs[i+1].s - sigs[i].s;
        		sigs[j] = SignatureSimple(newh,news);
        		j++;
        	}
        }
        sigs.resize(j);
        printf("sorting done; %lu elements with h < 2^%u obtained\n",sigs.size(),threshold_bit);
        if (out_prefix.length()) {
            // file format: redsigs_sd_round-i.bin
            std::string outsig = out_prefix + "_round-" + std::to_string(round) + ".bin";
            printf("saving signatures to %s... \n", outsig.c_str());
            sigsave(sigs, outsig);
        }
    }
    return;
}

void schroeppel_shamir_mpi(std::vector<SignatureSimple>& sigs,
		const uint32_t secpar, const uint32_t l, const double b,
		const uint32_t filter, const uint32_t a, const int log_prec, const std::vector<int>& ofst_info, const int iota=1, const std::string out_prefix="", const bool istest=true) {
    mpi::environment env;
    mpi::communicator world;
    const int master = 0;
    const int slave_last = world.size()-1;

    double threshold_bit_f = secpar - filter;
    const uint32_t keep_min = 1<<l;
    const uint32_t keep_max = (1<<l)*11/10;
    uint32_t S, q1, q2, q3;
    const std::string fname_prefix = "ss-mpi";
    if (ofst_info.size()!=iota) {
    	printf("ERROR: invalid offset info\n");
    	return;
    }

    for (int round = 0; round < iota; round++) {
		int pckofst = ofst_info[round];

		/* compute the next threshold value */
		threshold_bit_f -= b;
		uint32_t threshold_bit_z = threshold_bit_f; //rounded
		double extra = pow(2.0, threshold_bit_f-threshold_bit_z);
    	mpz_class threshold_mpz = mpz_class((mpf_class(1) << threshold_bit_z)*mpf_class(extra));
    	gss threshold;
    	mpz_to_gss(threshold, threshold_mpz, pckofst);

    	// file format: ss-mpi_round-i_{L1,R1,L2,R2}.bin
        std::string L1_fname, R1_fname, L2_fname, R2_fname;
        L1_fname = fname_prefix + "_round-" + std::to_string(round) + "_L1.bin";
        R1_fname = fname_prefix + "_round-" + std::to_string(round) + "_R1.bin";
        L2_fname = fname_prefix + "_round-" + std::to_string(round) + "_L2.bin";
        R2_fname = fname_prefix + "_round-" + std::to_string(round) + "_R2.bin";

        std::vector<gss> L1_h, R1_h, L2_h, R2_h;

        if (world.rank() == master) {
            /* split sigs into L1 || R1 || L2 || R2 */
            S = sigs.size();
			q1 = S/4;
			q2 = S/2;
			q3 = S*3/4;

			std::cout << "master: sorting L1 in descending order..." << std::endl;
			std::sort(sigs.begin(), sigs.begin()+q1, std::greater<SignatureSimple>());
			std::cout << "master: sorting R1 in ascending order..." << std::endl;
			std::sort(sigs.begin()+q1, sigs.begin()+q2);
			std::cout << "master: sorting L2 in descending order..." << std::endl;
			std::sort(sigs.begin()+q2, sigs.begin()+q3, std::greater<SignatureSimple>());
			std::cout << "master: sorting R2 in ascending order..." << std::endl;
			std::sort(sigs.begin()+q3, sigs.end());
			std::cout << "master: sorting done" << std::endl;

			/* save split sigs */
			std::cout << "master: saving split lists..." << std::endl;
			sigsave_it(sigs.begin(), sigs.begin()+q1, L1_fname);
			sigsave_it(sigs.begin()+q1, sigs.begin()+q2, R1_fname);
			sigsave_it(sigs.begin()+q2, sigs.begin()+q3, L2_fname);
			sigsave_it(sigs.begin()+q3, sigs.end(), R2_fname);

			/* conversion: mpz_t -> gss vector */
			std::cout << "master: preserving 128-bit of h..." << std::endl;
			packhalf_it(sigs.begin(), sigs.begin()+q1, L1_h, q1, pckofst);
			packhalf_it(sigs.begin()+q1, sigs.begin()+q2, R1_h, q2-q1, pckofst);
			packhalf_it(sigs.begin()+q2, sigs.begin()+q3, L2_h, q3-q2, pckofst);
			packhalf_it(sigs.begin()+q3, sigs.end(), R2_h, S-q3, pckofst);

			std::cout << "master: deleting the original GMP vector" << std::endl;
			std::vector<SignatureSimple>().swap(sigs);
		  	malloc_trim(0); // required to clean up the orphaned memory allocated by GMP
		  	printf("master: trying to get 2^%u values less than 2^%.2f\n"\
    	            "[keep_min, keep_max] = [%u, %u]\n", l, threshold_bit_f, keep_min, keep_max);
    	}
    	world.barrier();

    	/* broadcast sigs */
    	broadcast(world, L1_h, master);
    	broadcast(world, R1_h, master);
    	broadcast(world, L2_h, master);
    	broadcast(world, R2_h, master);
        std::vector<gss>(L1_h).swap(L1_h);
        std::vector<gss>(R1_h).swap(R1_h);
        std::vector<gss>(L2_h).swap(L2_h);
        std::vector<gss>(R2_h).swap(R2_h);
        malloc_trim(0);

    	printf("[%d]/[%d]: filled L1 with %lu signatures \n", world.rank(), world.size(), L1_h.size());
    	printf("[%d]/[%d]: filled R1 with %lu signatures \n", world.rank(), world.size(), R1_h.size());
    	printf("[%d]/[%d]: filled L2 with %lu signatures \n", world.rank(), world.size(), L2_h.size());
    	printf("[%d]/[%d]: filled R2 with %lu signatures \n", world.rank(), world.size(), R2_h.size());

    	/* job scheduling */
    	const uint32_t A_lim = 1<<a;
        uint32_t range;
		range = (A_lim/2 + world.size() - 1)/world.size();

        if (world.rank() == master)
            printf("master: looking for collisions on top %u-bits, range is %u\n", a, range);

        uint32_t A_col_sum = 0;
        size_t log_th = 0;
        float percent = 0;
        std::vector<Index> subresult ;
        const uint32_t keep_min_sub = keep_min/world.size();
        const uint32_t keep_max_sub = keep_max/world.size();
        const uint32_t keep_max_combs = (1<<a)*11/10;
        const uint32_t log_mod = keep_min_sub/log_prec;
        const float percent_delta = 100.0/log_prec;
        subresult.reserve(keep_max_sub);
#pragma omp parallel shared(world,subresult,L1_h,R1_h,L2_h,R2_h,A_col_sum,range,threshold,threshold_bit_f,pckofst,log_th,percent)
{
        std::vector<LRCombGSS> combs1;
        std::vector<LRCombGSS> combs2;
        std::vector<Index> combs_A;
        combs1.reserve(keep_max_combs);
        combs2.reserve(keep_max_combs);
        combs_A.reserve(5000);
		#pragma omp for schedule(dynamic)
        for (uint32_t k = 0; k < range; k++) {
        for (int rev = 0; rev < 2; rev++) {
			combs1.reserve(keep_max_combs);
			combs2.reserve(keep_max_combs);
			uint32_t A = world.size()*k + world.rank();
			if (rev==1)
				A = A_lim-A-1;

            if (A >= A_lim || subresult.size() > keep_max_sub) continue;

            /* Look for collisions on top a-bits */
            collect_lc_two_GSS(combs1, L1_h, R1_h, A, a, threshold_bit_f + b, pckofst, keep_max_combs);
            collect_lc_two_GSS(combs1, L1_h, R1_h, A+A_lim, a, threshold_bit_f + b, pckofst, keep_max_combs);
            if (combs1.size() == 0) continue;
            std::sort(combs1.begin(), combs1.end());
            collect_lc_two_GSS(combs2, L2_h, R2_h, A, a, threshold_bit_f + b, pckofst, keep_max_combs);
            collect_lc_two_GSS(combs2, L2_h, R2_h, A+A_lim, a, threshold_bit_f + b, pckofst, keep_max_combs);
            if (combs2.size() == 0) continue;
            std::sort(combs2.begin(), combs2.end());
            //printf("[%d]/[%d]: A=%u: got (%lu, %lu) collisions \n", world.rank(), world.size(), A, combs1.size(), combs2.size());
            //printf( "[%d]/[%d]@%s: thread %d checks A=%u.\n", world.rank(), world.size(), env.processor_name().c_str(), omp_get_thread_num (), A);

            uint32_t i=0;
            uint32_t j=0;
            bool flip;
            while (i < combs1.size() && j < combs2.size()) {
                flip = (gss_lt(&combs1[i].hsum, &combs2[j].hsum) == 1);
                gss h_diff;
                if (flip) gss_sub(&h_diff, &combs2[j].hsum, &combs1[i].hsum);
                else gss_sub(&h_diff, &combs1[i].hsum, &combs2[j].hsum);

                if (gss_lt(&h_diff, &threshold))
                    combs_A.emplace_back(combs1[i].idx_L,combs1[i].idx_R,combs2[j].idx_L,combs2[j].idx_R,flip);
                if (flip) i++;
                else j++;
            }
            //printf("[%d]/[%d]: A=%u: got %lu subresult\n", world.rank(), world.size(), A, combs_A.size());
#pragma omp critical
{
            subresult.insert(subresult.end(), combs_A.begin(), combs_A.end());
            A_col_sum += combs1.size();
            if (subresult.size() >= log_th) {
            	printf("[%d]/[%d]: %.2f %% done\n", world.rank(), world.size(), percent);
            	log_th += log_mod;
            	percent += percent_delta;
            }
}
			combs1.resize(0);
			combs2.resize(0);
			combs_A.resize(0);
        } // for rev
        } // for A
        std::vector<LRCombGSS>().swap(combs1);
        std::vector<LRCombGSS>().swap(combs2);
        std::vector<Index>().swap(combs_A);
}
        if (subresult.size() > keep_max_sub) subresult.resize(keep_max_sub);
        printf("[%d]/[%d]: mean number of collisions on A: %d\n", world.rank(), world.size(), static_cast<int>(A_col_sum/range));
        printf("[%d]/[%d]: got %lu subresult\n", world.rank(), world.size(), subresult.size());
        if (subresult.size() < keep_min_sub)
        	printf("WARNING:[%d]/[%d]: failed to get expected amount of reduced values!\n", world.rank(), world.size());


        /*
         * save interim report (only indices)
		 * file format: index_round-i_rank-j.bin
		 */
		std::string outidx = "index_round-" + std::to_string(round) + "_rank-" \
							+ std::to_string(world.rank()) + ".bin";
		printf("[%d]/[%d]: saving signatures of h < 2^%.2f to %s... \n", world.rank(), world.size(), threshold_bit_f, outidx.c_str());
		idxsave(subresult, outidx);

        // delete original list
        printf("[%d]/[%d]: deleting the original list \n", world.rank(), world.size());
        std::vector<gss>().swap(L1_h);
        std::vector<gss>().swap(R1_h);
        std::vector<gss>().swap(L2_h);
        std::vector<gss>().swap(R2_h);
        world.barrier();

        /* prepare subres_sizes for gatherv */
        std::vector<int> subres_sizes;
        gather(world, static_cast<int>(subresult.size()), subres_sizes, master);
        // casting is required since MPI gatherv does not allow size specification with unsigned int...

        /* gatherv subresult into result*/
        std::vector<Index> result;
        if (world.rank() == master) {
        	std::cout << "master: gathering results..." << std::endl;
            uint32_t res_size = 0;
            for (auto& subsize : subres_sizes)
                res_size += subsize;
            result.resize(res_size);
            gatherv(world, subresult, &result[0], subres_sizes, master);

        } else {
            gatherv(world, subresult, master);

        }
        std::vector<int>().swap(subres_sizes);
        std::vector<Index>().swap(subresult);

        /* Recover GMP format */
        if (world.rank() == master) {
        	std::cout << "master: computing linear combinations from indices..." << std::endl;
        	std::vector<SignatureSimple> L1, R1, L2, R2;
        	L1.reserve(q1);
        	R1.reserve(q2-q1);
        	L2.reserve(q3-q2);
        	R2.reserve(S-q3);
            sigload(L1, L1_fname);
            sigload(R1, R1_fname);
            sigload(L2, L2_fname);
            sigload(R2, R2_fname);
            sigs.reserve(result.size());
            restore_from_idx(sigs,result,L1,R1,L2,R2,threshold_mpz);


            printf("master: got %lu result \n", sigs.size());
            if (sigs.size() < keep_min) puts("WARNING: failed to get expected amount of reduced values!");
            if (out_prefix.length()) {
            	// file format: redsigs_round-i.bin
				std::string outsig = out_prefix + "_round-" + std::to_string(round) + ".bin";
				printf("master: saving signatures of h < 2^%.2f to %s... \n", threshold_bit_f, outsig.c_str());
				sigsave(sigs, outsig);
			}

            /* cleanups */
            std::vector<Index>().swap(result);
		  	std::vector<SignatureSimple>().swap(L1);
		  	std::vector<SignatureSimple>().swap(R1);
		  	std::vector<SignatureSimple>().swap(L2);
		  	std::vector<SignatureSimple>().swap(R2);
		  	malloc_trim(0);
            puts("-------------------------");
        }
    }// for round
    if (world.rank() == master)
    	printf("master: reduction ended after %d rounds\n", iota);
}

void restore_from_idx(std::vector<SignatureSimple>& sigs, const std::vector<Index>& idxlist,
			 const std::vector<SignatureSimple>& L1, const std::vector<SignatureSimple>& R1,
			 const std::vector<SignatureSimple>& L2, const std::vector<SignatureSimple>& R2,
			 const mpz_class& threshold_mpz) {
	for (auto& r: idxlist) {
		mpz_class h = L1[r.idx_L1].h + R1[r.idx_R1].h - L2[r.idx_L2].h - R2[r.idx_R2].h;
		mpz_class s = L1[r.idx_L1].s + R1[r.idx_R1].s - L2[r.idx_L2].s - R2[r.idx_R2].s;
		if (abs(h) >= threshold_mpz) {
			printf("WARNING: found h of %lu-bit at (%u, %u, %u, %u), skipping\n",
				mpz_sizeinbase(h.get_mpz_t(), 2), r.idx_L1, r.idx_R1, r.idx_L2, r.idx_R2);
			continue;
		}
		if (r.flip) sigs.emplace_back(-h,-s);
		else sigs.emplace_back(h,s);
	}
}
