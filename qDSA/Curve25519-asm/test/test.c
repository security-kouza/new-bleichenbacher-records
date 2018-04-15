#include "../sign.h"
#include "../dh.h"
#include "print.h"

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define mes_len 16
#define loop_min 0
#define loop_len 100000

size_t readhex(unsigned char* out, char *s, size_t maxlen) {
    size_t len;
    memset(out, 0, maxlen);
    for(len=0; len < maxlen && *s; len++) {
       sscanf(s, "%2hhx", &out[len]);
       s+=2;
    }
    
    return len;
}

void printhex(char* prefix, unsigned char* in, size_t len) {
    printf("%s", prefix);
    for(size_t i = 0; i < len; i++)
	printf("%02X%s", in[i], (i+1<len)?"":"\n");
}

int main(int argc, char *argv[])
{
    unsigned long long i;

    /* signatures API */

    int j;
    unsigned char sk[64]; 
    unsigned long long mlen = mes_len;
    unsigned char m[mlen]; 
    unsigned char pk[32];
    unsigned char sm[64+mlen];
    unsigned long long smlen;
    int ch;

    if(argc == 3) {
	readhex(sk, argv[1], 32);
	readhex(m,  argv[2], mlen);
	keypair(pk, sk);
        sign(sm, &smlen, m, mlen, pk, sk);
        ch = verify(m, mlen, sm, smlen, pk);

	printhex("sk = ", sk, 64);
	printhex("pk = ", pk, 32);
	printhex("m  = ", m,  mlen);
	printhex("sm = ", sm, smlen);
	printf("verification: %sok\n", ch?"":"not ");
	return ~ch;
    }

    srand( time(NULL) );
    printf("Checking %d signatures...\n", loop_len);
    for(j=loop_min;j<=loop_min+loop_len-1;j++)
    {
        srand(j);

        for(i=0;i<64;i++) { sk[i] = rand() % 256; }
        for(i=0;i<mlen;i++) { m[i] = rand() % 256; }

        keypair(pk, sk);
        sign(sm, &smlen, m, mlen, pk, sk);
        ch = verify(m, mlen, sm, smlen, pk);

        if ( ch != 1 || j%10000==0 ) 
        {
            printf("============\nSignature %d (%s)\nsk = ", j,
		    ch?"valid":"invalid");
            for(i=0;i<32;i++) { printf("%02X:", sk[i]); }
            printf("\npk = ");
            for(i=0;i<32;i++) { printf("%02X:", pk[i]); }
            printf("\nm  = ");
            for(i=0;i<mlen;i++) { printf("%02X:", m[i]); }
            printf("\nsm = ");
            for(i=0;i<smlen;i++) { printf("%02X:", sm[i]); }
            printf("\n");
        }
    }
    printf("Finished");

    /* end signatures API */

    /* DH */

#if 0
    unsigned char dh_sk[32];
    unsigned char dh_pk[32];
    unsigned char dh_ss[32];
    for(i=0;i<32;i++)
    {
        dh_sk[i] = 0;
        dh_pk[i] = 0;
    }

    dh_keygen(dh_pk, dh_sk);

    printf("\n");
    for(i=0;i<32;i++) { printf("%02X ", dh_pk[i]); }
    printf("\n");

    for(i=0;i<32;i++) { dh_sk[i] = 0; }
    dh_sk[0] = 2;
    dh_exchange(dh_ss, dh_pk, dh_sk);

    printf("\n");
    for(i=0;i<32;i++) { printf("%02X ", dh_ss[i]); }
    printf("\n");
#endif

    /* end DH */

    return 0;
}
