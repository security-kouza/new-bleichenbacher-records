//============================================================================
// Name        : attack_mpi.cpp
// Author      : Akira Takahashi and Mehdi Tibouchi
// Version     : 0.1
// Copyright   : Public domain
// Description : Parallelized Bleichenbacher's attack
// Standard    : C++11
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <bitset>

#include <gmpxx.h>
#include <boost/mpi.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/serialization/string.hpp>
#include <boost/program_options.hpp>

#include "mocksig.h"
#include "reduction.h"
#include "fft.h"

namespace mpi = boost::mpi;
namespace po = boost::program_options;

int main(int argc, char** argv) {
    /* MPI */
    mpi::environment env;
    mpi::communicator world;
    const int master = 0;

    /* parse command line arguments */
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "help message")
		("verbose", "verbose logging")
        ("test", "run in test mode")
		("red", "perform reduction")
		("sd", "use sort-and-difference")
		("fft", "perform key recovery")
		("known", po::value<int>(), "specify known MSBs of the secret key")
		("gamma", po::value<int>(), "gamma")
        ("in", po::value<std::string>(), "load signature data from a file")
		("out", po::value<std::string>(), "save reduced signature to a file")
		("fcount", po::value<int>(), "number of input files")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return EXIT_SUCCESS;
    }

	int log_prec = 100;
    if (vm.count("verbose")) {
        log_prec *= 100;
    }

    /* Parameters */
    uint32_t secpar;
    mpz_class n;
    uint32_t a;
    uint32_t filter;
    double b;
    uint32_t gamma;
    int iota;
    uint32_t l;
    Domain pp;
    mpz_class d;
    std::vector<SignatureSimple> sigs;
    bool testmode = false;
    std::vector<int> ofst_info;
    uint32_t ub,kb;
    mpz_class d_hi,d_lo;

    if (vm.count("known")) {
        kb = vm["known"].as<int>();
    } else {
        kb = 0;
    }
    if (vm.count("gamma")) {
        gamma = vm["gamma"].as<int>();
    } else {
        gamma = 0;
    }

    /* Set parameters */
    if (vm.count("test")) {
        /* Use signature data from input file */
        secpar = 90;
        n = (mpz_class(1) << secpar) - 33;
        a = 12;
        filter = 10;
        iota = 2;
        l = a+2;
        b = (double)(secpar-filter-kb-l)/iota;
        pp = mock::setup(secpar, n);
        d = mpz_class("924408261565060156037890712");
        ofst_info = {0,0};
        testmode = true;
    } else if (vm.count("sd")){
        /* Use real qDSA */
        secpar = 252;
        n = (mpz_class(1) << 252) + mpz_class("27742317777372353535851937790883648493");
        /* parameters for 3-bit bias */
        a = 28;
        filter = 0;
        iota = 8;
        ofst_info = {};
        l = a+2;
        b = a+2-gamma;
        pp = mock::setup(secpar, n);
        d = mpz_class("5220582922658643192668885191618908575833980181104027493552863441828733052420");
    } else {
        /* Use real qDSA */
        secpar = 252;
        n = (mpz_class(1) << 252) + mpz_class("27742317777372353535851937790883648493");
        /* parameters for 3-bit bias */
        /*
        a = 21;
        filter = 0;
        iota = 4;
        ofst_info = {2,2,1,0};
         */
        /* parameters for 2-bit bias */
        /*
        a = 24;
        filter = 19;
        iota = 3;
        ofst_info = {2,1,0};
        
        l = a+2;
        b = (double)(secpar-filter-kb-l)/iota;
        */
        /* parameters for reduction experiments: test252_a15_b0_f0 */
        a = 15;
        filter = 0;
        iota = 5;
        ofst_info = {2,2,1,0,0};
        l = a+2;
        b = 3*a -1.59;
        pp = mock::setup(secpar, n);
        d = mpz_class("5220582922658643192668885191618908575833980181104027493552863441828733052420");
    }
    gmp_printf("pp=(secpar=%u, n=%Zd),\n        d=%Zd\n", pp.secpar, pp.n.get_mpz_t(), d.get_mpz_t());
    //std::cin.ignore();
    /* Load signature data */
    sigs.reserve(1<<l);
    if (vm.count("in")) {
        if (world.rank() == master) {
            std::string prefix = vm["in"].as<std::string>();
            int fcount = 1;
            if (vm.count("fcount")) {
            	fcount = vm["fcount"].as<int>();
				for (int fc = 0; fc < fcount; fc++) {
					std::string fname = prefix + "_" + std::to_string(fc) + ".bin";
					printf("master: loading %s\n", fname.c_str());
					sigload(sigs, fname);
				}
            } else {
            	std::string fname = prefix + ".bin";
				printf("master: loading %s\n", fname.c_str());
				sigload(sigs, fname);
            }
            /* Preprocess sigs for recovering remaining bits */
			ub = secpar-kb;
			d_hi = (d>>ub)<< ub;
			d_lo = d-d_hi;
			for (SignatureSimple& sig : sigs) {
				mpz_class tmp = sig.s + sig.h*d_hi;
				mpz_mod(sig.s.get_mpz_t(), tmp.get_mpz_t(), n.get_mpz_t());
			}
            printf("master: loaded %lu signatures\n", sigs.size());
        }
    }
    else { return EXIT_SUCCESS;}

    std::string out_prefix = "redsigs"; // set default output file prefix
    if (vm.count("out"))
    	out_prefix = vm["out"].as<std::string>();

    /* REDUCTION */
    if (vm.count("red")) {
        if (vm.count("sd")) {
            sort_and_difference(sigs, secpar, l, b, filter, a, log_prec, ofst_info, iota, "", testmode);
        } else {
            schroeppel_shamir_mpi(sigs, secpar, l, b, filter, a, log_prec, ofst_info, iota, out_prefix, testmode);
        }
    }
    /* FFT */
    if (world.rank() == master && vm.count("fft")) {
        const uint32_t L = sigs.size();
        printf("master: got %u reduced sigs\n", L);

        /* Debug Logging */
        /*
        std::sort(sigs.begin(), sigs.end());
        for (const SignatureSimple& s : sigs) {
            sigprint(s);
        }
        */

        /* FFT */
        double* W = compute_bias(L, kb, pp, sigs);
        uint32_t peak_at = find_peak(W, L);
        double noise = compute_noise_avg(L, peak_at, W);
        mpz_class estim = (mpz_class)peak_at*(mpz_class(1)<<ub)/L + d_hi;
        gmp_printf("sk w = %Zd \n"
                   "sk d = %Zd \n", estim.get_mpz_t(), d.get_mpz_t());
        printf("average noise = %lf\n"
               "estimated noise 1/sqrt(L) = %lf\n", noise, 1/sqrt(L));

		std::string bias_fname("bias.csv");
		printf("saving bias to file %s\n", bias_fname.c_str());
		std::ofstream bias_file(bias_fname.c_str(), std::ofstream::trunc);
		for(uint32_t i=0; i<L; i++) {
			bias_file << i << "," << W[i] << std::endl;
		}
		bias_file.close();

        std::string dbin = d.get_str(2);
        std::string ebin = estim.get_str(2);
        int count_msb = 0;
        int idx = 0;
        while (1) {
        	if (dbin[idx] != ebin[idx])
        		break;
        	count_msb++;
        	idx++;
        }
        std::cout << dbin << std::endl;
        std::cout << ebin << std::endl;
        printf("recovered %d-MSBs of sk d\n", count_msb);
    }

	return EXIT_SUCCESS;

}


