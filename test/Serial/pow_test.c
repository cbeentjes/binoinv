// Quick test how to best calculate (1-p)^(-N)
// 1) Repeated (1-p)^(-N)
// 2) Exp(-N*log(1-p))
// 3) pow(1-p, -N)

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

// Option 1
inline static double f1(double p, double q, double N) {
    double a = 1.0;
    for (int k = 0; k < N; k++) a *= q;
    a = 1.0/a;
    return a;
}

inline static float f1f(float p, float q, float N) {
    float a = 1.0f;
    for (int k = 0; k < N; k++) a *= q;
    a = 1.0f/a;
    return a;
}

// Option 2
inline static double f2(double p, double q, double N) {
    double a = -N*log1p(-p);
    if (a < 709.782) {
        a = exp(a);
    } else {
        a = exp(a - 128.0);
    }
    return a;
}

inline static float f2f(float p, float q, float N) {
    float a = -N*log1pf(-p);
    if (a < 88.72f) {
        a = expf(a);
    } else {
        a = expf(a - 64.0f);
    }
    return a;
}

// Option 3
inline static double f3(double p, double q, double N) {
    return pow(q, -N + p*1e-100);
}

inline static float f3f(float p, float q, float N) {
    return powf(q, -N + p*1e-20f);
}

int main(int argc, char *argv[]) {
    int r;
    int num_repeats = 50000000;

    double* P;
    double* Q;
    double* F1;
    double* F2;
    double* F3;

    P = (double*) malloc(num_repeats*sizeof(double));
    Q = (double*) malloc(num_repeats*sizeof(double));
    F1 = (double*) malloc(num_repeats*sizeof(double));
    F2 = (double*) malloc(num_repeats*sizeof(double));
    F3 = (double*) malloc(num_repeats*sizeof(double));

    float* Pf;
    float* Qf;
    float* F1f;
    float* F2f;
    float* F3f;

    Pf = (float*) malloc(num_repeats*sizeof(float));
    Qf = (float*) malloc(num_repeats*sizeof(float));
    F1f = (float*) malloc(num_repeats*sizeof(float));
    F2f = (float*) malloc(num_repeats*sizeof(float));
    F3f = (float*) malloc(num_repeats*sizeof(float));

    double N = 3.0;

    // change N based on input if desired
    if (argc > 1) {
        N = strtod(argv[1], NULL);
    }

    float Nf = N;

    float mu_0 = 10.0f;
    // Generate random values for p in [0,0.5]
    srand(time(NULL));   // Initialization, should only be called once.
    for (int k = 0; k < num_repeats; k++) {
        r = rand();
        P[k]  = 0.5*r/RAND_MAX;
        Pf[k] = 0.5f*r/RAND_MAX;
        // Check when NP < mu_0
        P[k]  *= mu_0/N;
        Pf[k] *= mu_0/Nf;

        Q[k] = 1.0 - P[k];
        Qf[k] = 1.0f - P[k];
    }

    // Double precision
    printf("N = %d.\n",(int) N);
    printf("Time elapsed for %d calls.\n", num_repeats);
    clock_t t;
    double time_taken;

    if (N < 50.0) {
        t = clock();
        for (int k = 0; k < num_repeats; k++) {
            F1[k] = f1(P[k],Q[k],N);
        }
        t = clock() - t;
        time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
        printf("f1: %1.2e (s)\n", time_taken);
    }

    t = clock();
    for (int k = 0; k < num_repeats; k++) {
        F2[k] = f2(P[k],Q[k],N);
    }
    t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("f2: %1.2e (s)\n", time_taken);

    t = clock();
    for (int k = 0; k < num_repeats; k++) {
        F3[k] = f3(P[k],Q[k],N);
    }
    t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("f3: %1.2e (s)\n", time_taken);

    // Single precision
    if (Nf < 50.0f) {
        t = clock();
        for (int k = 0; k < num_repeats; k++) {
            F1f[k] = f1f(Pf[k],Qf[k],Nf);
        }
        t = clock() - t;
        time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
        printf("f1f: %1.2e (s)\n", time_taken);
    }

    t = clock();
    for (int k = 0; k < num_repeats; k++) {
        F2f[k] = f2f(Pf[k],Qf[k],Nf);
    }
    t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("f2f: %1.2e (s)\n", time_taken);

    t = clock();
    for (int k = 0; k < num_repeats; k++) {
        F3f[k] = f3f(Pf[k],Qf[k],Nf);
    }
    t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("f3f: %1.2e (s)\n", time_taken);

    double m = 0.0;
    float mf = 0.0;
    for (int k = 0; k < num_repeats; k++) {
        if (fabs(1 - F2[k]/F3[k]) > m) {
            m = fabs(1-F2[k]/F3[k]);
        }
        if (fabsf(1 - F2f[k]/F3f[k]) > mf) {
            mf = fabsf(1-F2f[k]/F3f[k]);
        }
    }

    printf("Max relative error f2 (double): %1.10e\n",m);
    printf("Max relative error f2 (single): %1.10e\n",mf);

    if (Nf < 50.0f) {
        m  = 0.0;
        mf = 0.0f;
        for (int k = 0; k < num_repeats; k++) {
            if (fabs(1 - F1[k]/F3[k]) > m) {
                m = fabs(1-F1[k]/F3[k]);
            }
            if (fabsf(1 - F1f[k]/F3f[k]) > mf) {
                mf = fabsf(1-F1f[k]/F3f[k]);
            }
        }
        printf("Max relative error f1 (double): %1.10e\n",m);
        printf("Max relative error f1 (single): %1.10e\n",mf);
    }

    free(F3f);
    free(F2f);
    free(F1f);
    free(F3);
    free(F2);
    free(F1);
    free(Q);
    free(P);
}
