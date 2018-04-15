/*
 * siggen_qdsa.cpp
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <inttypes.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <bitset>

#include <gmpxx.h>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/program_options.hpp>

#include "mocksig.h"
#include "qdsawrapper.h"

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
        ("test", "generate test signatures")
        ("out", po::value<std::string>(), "save signature data to a file with specified prefix")
        ("filter", po::value<int>(), "number of bits filtered")
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
        log_prec *= 10;
    }

    /* Parameters */
    uint32_t secpar;
    mpz_class n,tmplim;
    sc25519 lim;
    uint32_t a;
    uint32_t filter;
    uint32_t leak;
    uint32_t S;
    Domain pp;
    mpz_class d;
    std::vector<SignatureSimple> subresult;

    if (vm.count("filter"))
        filter = vm["filter"].as<int>();
    else
        filter = 0;

    if (vm.count("test")) {
        /* Use mock signature */
        secpar = 90;
        n = (mpz_class(1) << secpar) - 33;
        tmplim = n >> filter;
        a = 12;
        leak = 2;
        S = 1<<(a+2);
        pp = mock::setup(secpar, n);
        d = mpz_class("924408261565060156037890712");
        gmp_printf("pp = (secpar=%lu, n=%Zd), d=%Zd\n", pp.secpar, pp.n.get_mpz_t(), d.get_mpz_t());
        gmp_randclass rand(gmp_randinit_default);
        rand.seed(1234567);

        if (world.rank() == master) {
            printf("master: generating %u signatures with |h| <= %u bits...\n", S, secpar - filter);
            uint32_t numsigs = 0;
            mpz_class div, inv, h_tmp, s_tmp;
            div = mpz_class(1) << leak;
            mpz_invert(inv.get_mpz_t(), div.get_mpz_t(), pp.n.get_mpz_t());
            while (subresult.size() < S) {
                numsigs++;
                SignatureLeak sig = mock::sign(pp, d, 0, leak, rand);
                h_tmp = sig.h*inv;
                mpz_mod(sig.h.get_mpz_t(), h_tmp.get_mpz_t(), pp.n.get_mpz_t());
                if (sig.h >= tmplim) continue;
                s_tmp = (sig.s - sig.rr)*inv;
                mpz_mod(sig.s.get_mpz_t(), s_tmp.get_mpz_t(), pp.n.get_mpz_t());
                subresult.emplace_back(SignatureSimple(sig.h, sig.s));
                if (subresult.size() % (S/log_prec) == 0) printf("%.2f %% done\n", (float)subresult.size()*100/S);
            }
            printf("master: got %u/%u signatures\n", S, numsigs);
        }
    } else {
        /* Use real qDSA */
        secpar = 252;
        n = (mpz_class(1) << 252) + mpz_class("27742317777372353535851937790883648493");
        tmplim = n >> filter;
        mpz_to_gs(lim, tmplim);
        a = 24;
        leak = 2;
        S = 1 << (a+2);
        pp = mock::setup(secpar, n);
        const unsigned long long mlen = sizeof(uint64_t) + sizeof(int); // should be > (a+2+filter)/8
        size_t S_sub;
        if (world.rank() == master) {
            S_sub = (S + world.size() - 1)/world.size();
            printf("master : generating %u signatures...\n", S);
        }
        broadcast(world, S_sub, master);
        const size_t log_mod = S_sub/log_prec;
        const float percent_delta = 100.0/log_prec;
        size_t log_th = 0;
        float percent = 0;

        /* KeyGen */
        uint8_t pk[32];
        uint8_t sk[64] = {0x8D, 0x1E, 0xBD, 0x7D, 0xF8, 0xAE, 0xF4, 0x53, 0x02, 0x49, 0xD3, 0x26, 0x17, 0x58, 0x7A, 0xAC,
                          0x44, 0x09, 0xA7, 0x79, 0x34, 0x5F, 0x1E, 0x6B, 0x20, 0x1C, 0xFC, 0x7C, 0xC5, 0x7E, 0x3C, 0x53,
                          0x9D, 0xF9, 0xD0, 0x95, 0xA7, 0xC4, 0xE9, 0xA9, 0x0D, 0xBC, 0xD0, 0x24, 0x14, 0x4A, 0xD0, 0x58,
                          0x54, 0x78, 0xD1, 0x88, 0xD7, 0xF0, 0xF4, 0xF7, 0x0C, 0xF0, 0x73, 0xD2, 0x6E, 0xAF, 0x25, 0x0B};
        // This leads to d = 5220582922658643192668885191618908575833980181104027493552863441828733052420
        d = qdsa::keygen(pk, sk);
        printf("d is %lu-bit\n", mpz_sizeinbase(d.get_mpz_t(), 2));
        mpz_mod(d.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());
        gmp_printf("pp=(secpar=%u, n=%Zd),\n        d=%Zd\n", pp.secpar, pp.n.get_mpz_t(), d.get_mpz_t());

        /* SigGen */
        printf("[%d]/[%d] : generating %lu signatures with |h| <= %u bits...\n", world.rank(), world.size(), S_sub, secpar - filter);
        uint64_t numsigs = 0;

		union {
        	unsigned char m[mlen]; // message
			struct {
        		uint64_t index;
        		int rank;
			} mstruct;
		};
        unsigned char sm[64+mlen]; // +-R || s || m
        unsigned long long smlen;
	
        mstruct.rank = world.rank();
        while (subresult.size() < S_sub) {
            numsigs++;
            /* generate random message */
            mstruct.index = numsigs;
#if 0
			if(numsigs==17) {
			printf("[%d]/[%d] : m[%u] = ", world.rank(), world.size(), numsigs);
			for(unsigned long long i=0;i<mlen;i++)
				printf("%02x%c", m[i], (i==mlen-1)?'\n':':');
			}
#endif

            /* sign */
            SignatureSimple sig;
            if (qdsa::sign(sig, sm, &smlen, m, mlen, pk, sk, pp.n, &lim) != 1) continue;
            subresult.emplace_back(sig);
            if (subresult.size() >= log_th) {
                printf("[%d]/[%d]: %.2f %% done\n", world.rank(), world.size(), percent);
                log_th += log_mod;
                percent += percent_delta;
            }

            /* verify signature */
#if 0
            int ch = qdsa::verify(m, mlen, sm, smlen, pk);
            if (ch != 1) printf("WARNING: invalid signature\n");
#endif
        }
        printf("[%d]/[%d]: got %lu/%lu signatures\n", world.rank(), world.size(), S_sub, numsigs);
    }

    if (vm.count("out")){
        std::string outsig = vm["out"].as<std::string>() + "_" + std::to_string(world.rank()) + ".bin";
        printf("master: saving %lu signatures in %s... \n", subresult.size(), outsig.c_str());
        sigsave(subresult, outsig);
    }
    return EXIT_SUCCESS;
}
