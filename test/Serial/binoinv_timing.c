//
// standard header files
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// https://gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
// variadic macro to print to both file and stdout
#define PRINTF2(fp, ...) {printf(__VA_ARGS__); fprintf(fp, __VA_ARGS__);}

//
// binoinv header files
//

#include "binoinv.h"
#include "binoinvf.h"

//
// linux timing routine
//

#include <sys/time.h>

double elapsed_time(double *et) {
    struct timeval t;
    double old_time = *et;

    gettimeofday( &t, (struct timezone *)0 );
    *et = t.tv_sec + t.tv_usec*1.0e-6;

    return *et - old_time;
}

//
// CDF inverse timing test functions
//

void normcdfinv_test_double_scalar(int M) {
    double x, u;
    int m;

    for (m=0; m<M; m++) {
        // Generate "random" variates
        u = (m+0.5) / M;

#ifdef __INTEL_MKL__
        vdCdfNormInv(1, &u, &x);
#else
        x = norminv(u);
#endif

        // needed to prevent compiler discarding everything
        if (x==-999.0) printf("negative x\n");
    }
}

// Double precision CPU/MIMD-style version
void binoinv_test_double_scalar(int M, double N, double p) {
    double x, u;
    int m;

    for (m=0; m<M; m++) {
        // Generate "random" variates
        u = (m+0.5) / M;

        // extra m*1e-100 is to prevent compiler
        // optimising for fixed N and p, need to take p >~ 1e-80 to retain
        // correctness despite m*1e-100
        x = binoinv(u, N + m*1e-100, p + m*1e-100);

        // needed to prevent compiler discarding everything
        if (x==-999.0) printf("negative x\n");
    }
}

// Double precision GPU/SIMD-style version
void binoinv_test_double_vector(int M, double N, double p) {
    double x, u;
    int m;

    for (m=0; m<M; m++) {
        // Generate "random" variates
        u = (m+0.5) / M;

        // extra m*1e-100 is to prevent compiler
        // optimising for fixed N and p, need to take p >~ 1e-80 to retain
        // correctness despite m*1e-100
        x = binoinv_v(u, N + m*1e-100, p + m*1e-100);

        // needed to prevent compiler discarding everything
        if (x==-999.0) printf("negative x\n");
    }
}

// Single precision version
void normcdfinv_test_float_scalar(int M) {
    float x, u;
    int m;

    for (m=0; m<M; m++) {
        // Generate "random" variates
        u = (m+0.5) / M;

#ifdef __INTEL_MKL__
        vsCdfNormInv(1, &u, &x);
#else
        x = norminvf(u);
#endif

        // needed to prevent compiler discarding everything
        if (x==-999.0) printf("negative x\n");
    }
}

// Single precision version CPU/MIMD-style version
void binoinv_test_float_scalar(int M, float N, float p) {
    float x, u;
    int m;

    for (m=0; m<M; m++) {
        // Generate "random" variates
        u = (m+0.5f) / M;

        // extra n*1e-20 is to prevent compiler
        // optimising for fixed p
        x = binoinvf(u, N, p + m*1.0e-20f);

        // needed to prevent compiler discarding everything
        if (x==-999.0f) printf("negative x\n");
    }
}

// Single precision version GPU/SIMD-style version
void binoinv_test_float_vector(int M, float N, float p) {
    float x, u;
    int m;

    for (m=0; m<M; m++) {
        // Generate "random" variates
        u = (m+0.5f) / M;

        // extra n*1e-20 is to prevent compiler
        // optimising for fixed p
        x = binoinvf_v(u, N, p + m*1.0e-20f);

        // needed to prevent compiler discarding everything
        if (x==-999.0f) printf("negative x\n");
    }
}

//
// main code
//

int main(int argc, char **argv) {
    double timer, elapsed; // timer variable and elapsed time
    float N, p;
    int M, pass, count, Count=10;

    char filename[128];
    FILE *fp;

    if (argc > 3) {
        sprintf(filename, "%s", argv[3]);
    } else {
        sprintf(filename, "binoinv_timing.txt");
    }
    fp = fopen(filename, "w");

#define REPEAT 4

    M = (1<<24);

    // execute code

    for (pass=0; pass<4; pass++) {
        // Option to test only double prec or single prec codes
#ifdef DOUBLE
        if (pass % 2 == 0) continue;
#endif
#ifdef SINGLE
        if (pass % 2 == 1) continue;
#endif

        // Option to only test scalar or vector codes
#ifndef VECTOR
        if (pass > 1) break;
#endif
#ifdef SKIP_SCALAR
        if (pass < 2) continue;
#endif

        // default parameter values
        N = 1.0f;
        p = 0.5f;

        // change N based on input if desired
        if (argc > 1) {
            N = strtod(argv[1], NULL);
        }
        // change p based on input if desired
        if (argc > 2) {
            p = strtod(argv[2], NULL);
        }

        if (pass==0) {
            PRINTF2(fp, "\nscalar single precision algorithm performance tests (CPU)\n");
        } else if (pass==1) {
            PRINTF2(fp, "\nscalar double precision algorithm performance tests (CPU)\n");
        } else if (pass==2) {
            PRINTF2(fp, "\nvector single precision algorithm performance tests (CPU) \n");
        } else {
            PRINTF2(fp, "\nvector double precision algorithm performance tests (CPU)\n");
        }
        PRINTF2(fp, "-------------------------------------------------------\n");
        PRINTF2(fp, "   N       p     execution time    samples/sec \n");

        // binoinv
        for (count=0; count<Count; count++) {
            elapsed_time(&timer);   // initialise timer

            // average over REPEAT runs
            for (int i=0; i<REPEAT; i++) {
                if (pass == 0) {
                    binoinv_test_float_scalar(M, N, p);
                } else if (pass==1) {
                    binoinv_test_double_scalar(M, (double) N, (double) p);
                } else if (pass==2) {
                    binoinv_test_float_vector(M, N, p);
                } else {
                    binoinv_test_double_vector(M, (double) N, (double) p);
                }
            }

            elapsed = elapsed_time(&timer);

            if (count>0) { // skip first one (cache effects?)
                PRINTF2(fp, "%6g  %7.3g     %9.4f      %10.3g \n",
                        N, p, elapsed, REPEAT*M/elapsed);
            }
            // update parameter
#ifdef PLUS
            N += 1.0;
#else
            N *= 2.0;
#endif
        }

        /* norminv */
        /* Run once first (cache effects?) */
        if (pass == 1 || pass == 3) {
            normcdfinv_test_double_scalar(M);
        } else {
            normcdfinv_test_float_scalar(M);
        }
        elapsed_time(&timer);   // initialise timer
        for (int i=0; i<REPEAT; i++) {
            if (pass == 1 || pass == 3) {
                normcdfinv_test_double_scalar(M);
            } else {
                normcdfinv_test_float_scalar(M);
            }
        }
        elapsed = elapsed_time(&timer);
        PRINTF2(fp, "   normcdfinv       %9.4f      %10.3g \n",
                elapsed, REPEAT*M/elapsed);
    }

    fclose(fp);
    return 0;
}
