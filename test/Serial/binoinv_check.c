//
// standard header files
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

// https://gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
// variadic macro to print to both file and stdout
#define PRINTF2(fp, ...) {printf(__VA_ARGS__); fprintf(fp, __VA_ARGS__);}

//
// binoinv header file
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
// function prototype for quad precision evaluation of binomial inverse CDF
//

void binoinv_quad(float, float, float*, float*, double*, double*);
void binocinv_quad(float, float, float*, float*, double*, double*);

//
// functions to calculate the floating point jump points of the approximate
// inverse CDF of the binomial distribution
//


// Double precision CPU/MIMD code version
void binoinv_bisection_scalar(int n, double N, double p,
        double *ulo_d, double *uhi_d) {
    double x, xt, u_lo, u_hi, u_mid;


    u_hi  = 1.0;
    u_lo  = 0.0;
    u_mid = 0.5*(u_hi + u_lo);
    xt    = (double) n;

    while ((u_mid > u_lo) & (u_mid < u_hi)) {
        x = binoinv(u_mid, N, p);

        if (x > xt)
            u_hi = u_mid;
        else
            u_lo = u_mid;

        u_mid = 0.5*(u_hi + u_lo);
    }
    ulo_d[n] = u_lo;
    uhi_d[n] = u_hi;
}

// Double precision GPU/SIMD code version
void binoinv_bisection_vector(int n, double N, double p,
        double *ulo_d, double *uhi_d) {
    double x, xt, u_lo, u_hi, u_mid;

    u_hi  = 1.0;
    u_lo  = 0.0;
    u_mid = 0.5*(u_hi + u_lo);
    xt    = (double) n;

    while ((u_mid > u_lo) & (u_mid < u_hi)) {
        x = binoinv_v(u_mid, N, p);

        if (x > xt)
            u_hi = u_mid;
        else
            u_lo = u_mid;

        u_mid = 0.5*(u_hi + u_lo);
    }
    ulo_d[n] = u_lo;
    uhi_d[n] = u_hi;
}

// Single precision CPU/MIMD code version
void binoinvf_bisection_scalar(int n, float N, float p,
        float *ulo_d, float *uhi_d) {
    float x, xt, u_lo, u_hi, u_mid;

    u_hi  = 1.0f;
    u_lo  = 0.0f;
    u_mid = 0.5f*(u_hi + u_lo);
    xt    = (float) n;

    while ((u_mid > u_lo) & (u_mid < u_hi)) {
        x = binoinvf(u_mid, N, p);

        if (x > xt)
            u_hi = u_mid;
        else
            u_lo = u_mid;

        u_mid = 0.5f*(u_hi + u_lo);
    }
    ulo_d[n] = u_lo;
    uhi_d[n] = u_hi;
}

// Single precision GPU/SIMD code version
void binoinvf_bisection_vector(int n, float N, float p,
        float *ulo_d, float *uhi_d) {
    float x, xt, u_lo, u_hi, u_mid;

    u_hi  = 1.0f;
    u_lo  = 0.0f;
    u_mid = 0.5f*(u_hi + u_lo);
    xt    = (float) n;

    while ((u_mid > u_lo) & (u_mid < u_hi)) {
        x = binoinvf_v(u_mid, N, p);

        if (x > xt)
            u_hi = u_mid;
        else
            u_lo = u_mid;

        u_mid = 0.5f*(u_hi + u_lo);
    }
    ulo_d[n] = u_lo;
    uhi_d[n] = u_hi;
}

//
// functions to calculate the floating point jump points of the approximate
// inverse complementary CDF of the binomial distribution
//

// Double precision CPU/MIMD code version
void binocinv_bisection_scalar(int n, double N, double p,
        double *ulo_d, double *uhi_d) {

    double x, xt, u_lo, u_hi, u_mid;

    u_hi  = 1.0;
    u_lo  = 0.0;
    u_mid = 0.5*(u_hi + u_lo);
    xt    = (double) n;

    while ((u_mid > u_lo) & (u_mid < u_hi)) {
        x = binocinv(u_mid, N, p);

        if (x < xt)
            u_hi = u_mid;
        else
            u_lo = u_mid;

        u_mid = 0.5*(u_hi + u_lo);
    }
    ulo_d[n] = u_lo;
    uhi_d[n] = u_hi;
}

// Single precision CPU/MIMD code version
void binocinvf_bisection_scalar(int n, float N, float p,
        float *ulo_d, float *uhi_d) {

    float x, xt, u_lo, u_hi, u_mid;

    u_hi  = 1.0f;
    u_lo  = 0.0f;
    u_mid = 0.5f*(u_hi + u_lo);
    xt    = (float) n;

    while ((u_mid > u_lo) & (u_mid < u_hi)) {
        x = binocinvf(u_mid, N, p);

        if (x < xt)
            u_hi = u_mid;
        else
            u_lo = u_mid;

        u_mid = 0.5f*(u_hi + u_lo);
    }
    ulo_d[n] = u_lo;
    uhi_d[n] = u_hi;
}

//
// function to compare the binoinv floating point jump points to exact jump
// points calculated using quad precision
//

void binoinv_err_compare(FILE *fp, float N, float p,
        double *Ulo_ex, double *Uhi_ex, double *Ulo_h, double *Uhi_h,
        float  *ulo_ex, float  *uhi_ex, float  *ulo_h, float  *uhi_h) {

    double err1, err1f;
    double timer, elapsed1, elapsed2;

    double NP = N*p;
    double NPQ = NP*(1.0-p);
    double NPQ_root = sqrt(NPQ);

    //////////////////////////////////////////////////
    // compute reference solution in quad precision
    //////////////////////////////////////////////////

    elapsed_time(&timer);
    binoinv_quad(N, p, ulo_ex, uhi_ex, Ulo_ex, Uhi_ex);
    elapsed1 = elapsed_time(&timer);

    //////////////////////////////////////////////////
    //  check double and single precision versions
    //////////////////////////////////////////////////

    err1  = 0.0;
    err1f = 0.0;

    // Let n go from 0 to N when p <= 0.5 and
    // let n go from N to 0 when p >  0.5
    int flag = (p > 0.5);
    int n;
    if (flag)
        n = N;
    else
        n = 0;

    elapsed_time(&timer);
    while (0 <= n && n <= N) {
#ifdef FAST_CHECK // Speed up check by using knowledge of distribution shape
        // Determine good starting point of the summation, as we do not need
        // to evaluate for n so small that it cannot be generated by binoinv(f)
        if (n < NP - 200.0 - 50.0*NPQ_root) {
            Ulo_h[n] = 0.0;
            Uhi_h[n] = DBL_TRUE_MIN;
            ulo_h[n] = 0.0f;
            uhi_h[n] = FLT_TRUE_MIN;
        }
        // Similarly, no need to continue if we are already beyond
        // double precision limit.
        else if (n > NP + 10.0 + 10.0*NPQ_root){
            Ulo_h[n] = 1.0-DBL_EPSILON/2.0;
            Uhi_h[n] = 1.0;
            ulo_h[n] = 1.0f-FLT_EPSILON/2.0f;
            uhi_h[n] = 1.0f;
        }
        else {
            binoinv_bisection_scalar(n,(double) N,(double) p, Ulo_h, Uhi_h);
            binoinvf_bisection_scalar(n, N, p, ulo_h, uhi_h);
        }
#else
        binoinv_bisection_scalar(n,(double) N,(double) p, Ulo_h, Uhi_h);
        binoinvf_bisection_scalar(n, N, p, ulo_h, uhi_h);
#endif /* FAST_CHECK */

        // Double precision
        err1  += 0.5*fabs( (Ulo_ex[n]-Ulo_h[n]) + (Uhi_ex[n]-Uhi_h[n]));

        // check if Linfinity error is <= 1. Note that this is only relevant
        // if ULP < 1
        if (n > 0 && n < 9007199254740992) {
            if ( Uhi_h[n] <= Ulo_ex[n-1] ) {
                PRINTF2(fp,"\n Double precision error \n");
                PRINTF2(fp,"\n error: N = %f, p = %f, n = %d, Uhi_h[n] = %20.16g, Ulo_ex[n-1] = %20.16g \n",
                        N,      p,      n,      Uhi_h[n],           Ulo_ex[n-1]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",Ulo_h[n], Uhi_h[n], Ulo_ex[n-1], Uhi_ex[n-1]);
                exit(1);
            }
        }
        if (n < N && n < 9007199254740992) {
            if ( Uhi_ex[n+1] <= Ulo_h[n] ) {
                PRINTF2(fp,"\n Double precision error \n");
                PRINTF2(fp,"\n error: N = %f, p = %f, n = %d, Uhi_ex[n+1] = %20.16g, Ulo_h[n] = %20.16g \n",
                        N,      p,      n,      Uhi_ex[n+1],           Ulo_h[n]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",Ulo_h[n], Uhi_h[n], Ulo_ex[n+1], Uhi_ex[n+1]);
                exit(1);
            }
        }

        // Single precision
        err1f += 0.5f*fabsf( (ulo_ex[n]-ulo_h[n]) + (uhi_ex[n]-uhi_h[n]));

        // check if Linfinity error is <= 1. Note that this is only relevant
        // if ULP < 1
        if (n > 0 && n < 16777216) {
            if ( uhi_h[n] <= ulo_ex[n-1] ) {
                PRINTF2(fp,"\n Single precision error \n");
                PRINTF2(fp,"\n error: N = %f, p = %f, n = %d, uhi_h[n] = %20.16g, ulo_ex[n-1] = %20.16g \n",
                        N,      p,      n,      uhi_h[n],           ulo_ex[n-1]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",ulo_h[n], uhi_h[n], ulo_ex[n-1], uhi_ex[n-1]);
                exit(1);
            }
        }
        if (n < N && n < 16777216) {
            if ( uhi_ex[n+1] <= ulo_h[n] ) {
                PRINTF2(fp,"\n Single precision error \n");
                PRINTF2(fp,"\n error: N = %f, p = %f, n = %d, uhi_ex[n+1] = %20.16g, ulo_h[n] = %20.16g \n",
                        N,      p,      n,      uhi_ex[n+1],           ulo_h[n]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",ulo_h[n], uhi_h[n], ulo_ex[n+1], uhi_ex[n+1]);
                exit(1);
            }
        }

        // Update n and check if we can break early
        if (flag) {
            if (n < NP - 200.0 - 50.0*NPQ_root) {
                break;
            }
            n -= 1;
        } else {
            if (n > NP + 10.0 + 20.0*NPQ_root) {
                break;
            }
            n += 1;
        }
    }
    elapsed2 = elapsed_time(&timer);
    // print parameters and L1 error
    PRINTF2(fp,"%9.0f %10.4g   %#9.2e     %#9.2e       %.2e    %.2e\n",N,p,err1,err1f,elapsed1,elapsed2);
}

//
// function to compare the binoinv_v floating point jump points to exact jump
// points calculated using quad precision
//

void binoinv_v_err_compare(FILE *fp, float N, float p,
        double *Ulo_ex, double *Uhi_ex, double *Ulo_h, double *Uhi_h,
        float  *ulo_ex, float  *uhi_ex, float  *ulo_h, float  *uhi_h) {

    double err1, err1f;
    double timer, elapsed1, elapsed2;

    double NP = N*p;
    double NPQ = NP*(1.0-p);
    double NPQ_root = sqrt(NPQ);

    //////////////////////////////////////////////////
    // compute reference solution in quad precision
    //////////////////////////////////////////////////

    elapsed_time(&timer);
    binoinv_quad(N, p, ulo_ex, uhi_ex, Ulo_ex, Uhi_ex);
    elapsed1 = elapsed_time(&timer);

    //////////////////////////////////////////////////
    //  check double and single precision versions
    //////////////////////////////////////////////////

    err1  = 0.0;
    err1f = 0.0;

    // Let n go from 0 to N when p <= 0.5 and
    // let n go from N to 0 when p >  0.5
    int flag = (p > 0.5);
    int n;
    if (flag)
        n = N;
    else
        n = 0;

    elapsed_time(&timer);
    while (0 <= n && n <= N) {
#ifdef FAST_CHECK // Speed up check by using knowledge of distribution shape
        // Determine good starting point of the summation, as we do not need
        // to evaluate for n so small that it cannot be generated by binoinv(f)
        if (n < NP - 200.0 - 50.0*NPQ_root) {
            Ulo_h[n] = 0.0;
            Uhi_h[n] = DBL_TRUE_MIN;
            ulo_h[n] = 0.0f;
            uhi_h[n] = FLT_TRUE_MIN;
        }
        // Similarly, no need to continue if we are already beyond
        // double precision limit.
        else if (n > NP + 10.0 + 10.0*NPQ_root){
            Ulo_h[n] = 1.0-DBL_EPSILON/2.0;
            Uhi_h[n] = 1.0;
            ulo_h[n] = 1.0f-FLT_EPSILON/2.0f;
            uhi_h[n] = 1.0f;
        } else {
            binoinv_bisection_vector(n,(double) N,(double) p, Ulo_h, Uhi_h);
            binoinvf_bisection_vector(n, N, p, ulo_h, uhi_h);
        }
#else
        binoinv_bisection_vector(n,(double) N,(double) p, Ulo_h, Uhi_h);
        binoinvf_bisection_vector(n, N, p, ulo_h, uhi_h);
#endif /* FAST_CHECK */

        // Double precision
        err1  += 0.5*fabs( (Ulo_ex[n]-Ulo_h[n]) + (Uhi_ex[n]-Uhi_h[n]));

        // check if Linfinity error is <= 1. Note that this is only relevant
        // if ULP < 1
        if (n > 0 && n < 9007199254740992) {
            if ( Uhi_h[n] <= Ulo_ex[n-1] ) {
                PRINTF2(fp,"\n Double precision error \n");
                PRINTF2(fp,"\n error: N = %f, p = %f, n = %d, Uhi_h[n] = %20.16g, Ulo_ex[n-1] = %20.16g \n",
                        N,      p,      n,      Uhi_h[n],           Ulo_ex[n-1]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",Ulo_h[n], Uhi_h[n], Ulo_ex[n-1], Uhi_ex[n-1]);
                exit(1);
            }
        }
        if (n < N && n < 9007199254740992) {
            if ( Uhi_ex[n+1] <= Ulo_h[n] ) {
                PRINTF2(fp,"\n Double precision error \n");
                PRINTF2(fp,"\n error: N = %f, p = %f, n = %d, Uhi_ex[n+1] = %20.16g, Ulo_h[n] = %20.16g \n",
                        N,      p,      n,      Uhi_ex[n+1],           Ulo_h[n]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",Ulo_h[n], Uhi_h[n], Ulo_ex[n+1], Uhi_ex[n+1]);
                exit(1);
            }
        }

        // Single precision
        err1f += 0.5f*fabsf( (ulo_ex[n]-ulo_h[n]) + (uhi_ex[n]-uhi_h[n]));

        // check if Linfinity error is <= 1. Note that this is only relevant
        // if ULP < 1
        if (n > 0 && n < 16777216) {
            if ( uhi_h[n] <= ulo_ex[n-1] ) {
                PRINTF2(fp,"\n Single precision error \n");
                PRINTF2(fp,"\n error: N = %f, p = %f, n = %d, uhi_h[n] = %20.16g, ulo_ex[n-1] = %20.16g \n",
                        N,      p,      n,      uhi_h[n],           ulo_ex[n-1]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",ulo_h[n], uhi_h[n], ulo_ex[n-1], uhi_ex[n-1]);
                exit(1);
            }
        }
        if (n < N && n < 16777216) {
            if ( uhi_ex[n+1] <= ulo_h[n] ) {
                PRINTF2(fp,"\n Single precision error \n");
                PRINTF2(fp,"\n error: N = %f, p = %f, n = %d, uhi_ex[n+1] = %20.16g, ulo_h[n] = %20.16g \n",
                        N,      p,      n,      uhi_ex[n+1],           ulo_h[n]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",ulo_h[n], uhi_h[n], ulo_ex[n+1], uhi_ex[n+1]);
                exit(1);
            }
        }

        // Update n and check if we can break early
        if (flag) {
            if (n < NP - 200.0 - 50.0*NPQ_root) {
                break;
            }
            n -= 1;
        } else {
            if (n > NP + 10.0 + 20.0*NPQ_root) {
                break;
            }
            n += 1;
        }
    }
    elapsed2 = elapsed_time(&timer);
    // print parameters and L1 error
    PRINTF2(fp,"%9.0f %10.4g   %#9.2e     %#9.2e       %.2e    %.2e\n",N,p,err1,err1f,elapsed1,elapsed2);
}

//
// function to compare the binocinv floating point jump points to exact jump
// points calculated using quad precision
//

void binocinv_err_compare(FILE *fp, float N, float p,
        double *Ulo_ex, double *Uhi_ex, double *Ulo_h, double *Uhi_h,
        float  *ulo_ex, float  *uhi_ex, float  *ulo_h, float  *uhi_h) {

    double err1, err1f;
    double timer, elapsed1, elapsed2; // timer variable and elapsed time

    double NP = N*p;
    double NPQ = NP*(1.0-p);
    double NPQ_root = sqrt(NPQ);

    //////////////////////////////////////////////////
    // compute reference solution in quad precision //
    // ////////////////////////////////////////////////

    elapsed_time(&timer);
    binocinv_quad(N, p, ulo_ex, uhi_ex, Ulo_ex, Uhi_ex);
    elapsed1 = elapsed_time(&timer);

    //////////////////////////////////////////////////
    //  check double and single precision versions
    //////////////////////////////////////////////////

    err1 = 0.0;
    err1f= 0.0f;

    // Let n go from 0 to N when p <= 0.5 and
    // let n go from N to 0 when p >  0.5
    int flag = (p > 0.5);
    int n;
    if (flag)
        n = N;
    else
        n = 0;

    elapsed_time(&timer);
    while (0 <= n && n <= N) {
#ifdef FAST_CHECK // Speed up check by using knowledge of distribution shape
        // Determine good starting point of the summation, as we do not need
        // to evaluate for n so small that it cannot be generated by binoinv(f)
        if (n > NP + 200.0 + 50.0*NPQ_root) {
            Ulo_h[n] = 0.0;
            Uhi_h[n] = DBL_TRUE_MIN;
            ulo_h[n] = 0.0f;
            uhi_h[n] = FLT_TRUE_MIN;
        }
        // Similarly, no need to continue if we are already beyond
        // double precision limit.
        else if (n < NP - 10.0 - 10.0*NPQ_root){
            Ulo_h[n] = 1.0-DBL_EPSILON/2.0;
            Uhi_h[n] = 1.0;
            ulo_h[n] = 1.0f-FLT_EPSILON/2.0f;
            uhi_h[n] = 1.0f;
        } else {
            binocinv_bisection_scalar(n,(double) N,(double) p, Ulo_h, Uhi_h);
            binocinvf_bisection_scalar(n, N, p, ulo_h, uhi_h);
        }
#else
        binocinv_bisection_scalar(n,(double) N,(double) p, Ulo_h, Uhi_h);
        binocinvf_bisection_scalar(n, N, p, ulo_h, uhi_h);
#endif /* FAST_CHECK */

        // Double precision
        err1  += 0.5*fabs( (Ulo_ex[n]-Ulo_h[n]) + (Uhi_ex[n]-Uhi_h[n]));
#ifdef DEBUG2
        printf("n = %d, err = %.8g\n",n,err1);
#endif

        // check if Linfinity error is <= 1. Note that this is only relevant
        // if ULP < 1
        if (n < N && n < 9007199254740992) {
            if ( Uhi_h[n] <= Ulo_ex[n+1] ) {
                PRINTF2(fp,"\n Double precision error \n");
                PRINTF2(fp,"\n error: N = %f, p = %f, n = %d, Uhi_h[n] = %20.16g, Ulo_ex[n-1] = %20.16g \n",
                        N,      p,      n,      Uhi_h[n],           Ulo_ex[n+1]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",Ulo_h[n], Uhi_h[n], Ulo_ex[n-1], Uhi_ex[n-1]);
                exit(1);
            }
        }
        if (n > 0 && n < 9007199254740992) {
            if ( Uhi_ex[n-1] <= Ulo_h[n] ) {
                PRINTF2(fp,"\n Double precision error \n");
                PRINTF2(fp,"\n error: N = %f, p = %f, n = %d, Uhi_ex[n-1] = %20.16g, Ulo_h[n] = %20.16g \n",
                        N,      p,      n,      Uhi_ex[n-1],           Ulo_h[n]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",Ulo_h[n], Uhi_h[n], Ulo_ex[n+1], Uhi_ex[n+1]);
                exit(1);
            }
        }

        // Single precision
        err1f += 0.5f*fabsf( (ulo_ex[n]-ulo_h[n]) + (uhi_ex[n]-uhi_h[n]));
#ifdef DEBUG2
        printf("n = %d, err = %.8g\n",n,err1f);
#endif

        // check if Linfinity error is <= 1. Note that this is only relevant
        // if ULP < 1
        if (n < N && n < 16777216) {
            if ( uhi_h[n] <= ulo_ex[n+1] ) {
                PRINTF2(fp,"\n Single precision error \n");
                PRINTF2(fp,"\n error: N = %f, p = %f, n = %d, uhi_h[n] = %20.16g, ulo_ex[n+1] = %20.16g \n",
                        N,      p,      n,      uhi_h[n],           ulo_ex[n-1]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",ulo_h[n], uhi_h[n], ulo_ex[n-1], uhi_ex[n-1]);
                exit(1);
            }
        }
        if (n > 0 && n < 16777216) {
            if ( uhi_ex[n-1] <= ulo_h[n] ) {
                PRINTF2(fp,"\n Single precision error \n");
                PRINTF2(fp,"\n error: N = %f, p = %f, n = %d, uhi_ex[n-1] = %20.16g, ulo_h[n] = %20.16g \n",
                        N,      p,      n,      uhi_ex[n-1],           ulo_h[n]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",ulo_h[n], uhi_h[n], ulo_ex[n+1], uhi_ex[n+1]);
                exit(1);
            }
        }

        // Update n and check if we can break early
        if (flag) {
            n -= 1;
        } else {
            if (n > NP + 200.0 + 50.0*NPQ_root) {
                break;
            }
            n += 1;
        }
    }
    elapsed2 = elapsed_time(&timer);
    // print parameters and L1 errors
    PRINTF2(fp,"%9.0f %10.4g   %#9.2e     %#9.2e       %.2e    %.2e\n",N,p,err1,err1f,elapsed1,elapsed2);
}


//////////////////////////////////////////////////
// main code
//////////////////////////////////////////////////

int main(int argc, char **argv) {

    // Destination file
    char filename[128];
    FILE *fp;

    if (argc > 3) {
        sprintf(filename, "%s", argv[3]);
    } else {
        sprintf(filename, "binoinv_check.txt");
    }
    fp = fopen(filename, "w");

    // Single precision variables
    float  *ulo_h, *uhi_h, *ulo_ex, *uhi_ex;
    // Double precision variables
    double *Ulo_h, *Uhi_h, *Ulo_ex, *Uhi_ex;


    // Binomial distribution parameters
    float p;
    float N;

    // allocate memory

    int Nmax = (1<<30);

    ulo_ex = (float  *)malloc(Nmax*sizeof(float));
    uhi_ex = (float  *)malloc(Nmax*sizeof(float));
    ulo_h  = (float  *)malloc(Nmax*sizeof(float));
    uhi_h  = (float  *)malloc(Nmax*sizeof(float));
    Ulo_ex = (double *)malloc(Nmax*sizeof(double));
    Uhi_ex = (double *)malloc(Nmax*sizeof(double));
    Ulo_h  = (double *)malloc(Nmax*sizeof(double));
    Uhi_h  = (double *)malloc(Nmax*sizeof(double));

    // number of loop iterations in accuracy tests
#ifndef COUNT_P
#define COUNT_P 16
#endif
#ifndef COUNT_N
#define COUNT_N  4
#endif
#ifndef COUNT_NP
#define COUNT_NP 4
#endif

#ifndef SKIP_SCALAR

    // set values to test
    p = 0.5f;
    N = 1.0f;
    // change N,p based on input if desired
    if (argc > 1) {
        N = strtod(argv[1], NULL);
    }
    if (argc > 2) {
        p = strtod(argv[2], NULL);
    }
    PRINTF2(fp,"\n");
    PRINTF2(fp,"\nscalar CPU algorithm accuracy tests (p fixed)\n");
    PRINTF2(fp,"----------------------------------------------------------------\n");
    PRINTF2(fp,"        N          p    error(FP64)   error(FP32)    T_quad(s)   T_num(s) \n");
    for (int count=0; count<COUNT_P; count++) {
        binoinv_err_compare(fp, N, p, Ulo_ex, Uhi_ex, Ulo_h, Uhi_h, ulo_ex, uhi_ex, ulo_h, uhi_h);
#ifdef PLUS
        N += 1.0f;
        if (N*p*(1.0f-p) > 8e+1f) break;
#else
        N *= 2.0f;
#endif
        if (N > Nmax) break;
    }

    // set values to test
    p = 1.0f;
    N = 1e3f;
    // change N,p based on input if desired
    if (argc > 1) {
        N = strtod(argv[1], NULL);
    }
    if (argc > 2) {
        p = strtod(argv[2], NULL);
    }
    PRINTF2(fp,"\nscalar CPU algorithm accuracy tests (N fixed)\n");
    PRINTF2(fp,"----------------------------------------------------------------\n");
    PRINTF2(fp,"        N          p    error(FP64)   error(FP32)    T_quad(s)   T_num(s) \n");
    for (int count=0; count<COUNT_N; count++) {
        binoinv_err_compare(fp, N, p, Ulo_ex, Uhi_ex, Ulo_h, Uhi_h, ulo_ex, uhi_ex, ulo_h, uhi_h);
        p *= 0.5f;
    }

    // set values to test
    p = 0.5f;
    N = 1e3f;
    // change N,p based on input if desired
    if (argc > 1) {
        N = strtod(argv[1], NULL);
    }
    if (argc > 2) {
        p = strtod(argv[2], NULL);
    }
    PRINTF2(fp,"\nscalar CPU algorithm accuracy tests (Np fixed)\n");
    PRINTF2(fp,"-----------------------------------------------------------------\n");
    PRINTF2(fp,"        N          p    error(FP64)   error(FP32)    T_quad(s)   T_num(s) \n");
    for (int count=0; count<COUNT_NP; count++) {
        binoinv_err_compare(fp, N, p, Ulo_ex, Uhi_ex, Ulo_h, Uhi_h, ulo_ex, uhi_ex, ulo_h, uhi_h);
        N *= 2.0f;
        p *= 0.5f;
        if (N > Nmax) break;
    }

#endif

    //////////////////////////////////////////////////
    // re-do for vector version
    //////////////////////////////////////////////////

#ifdef VECTOR

    // set values to test
    p = 0.5f;
    N = 1.0f;
    // change N,p based on input if desired
    if (argc > 1) {
        N = strtod(argv[1], NULL);
    }
    if (argc > 2) {
        p = strtod(argv[2], NULL);
    }
    PRINTF2(fp,"\n");
    PRINTF2(fp,"\nvector CPU algorithm accuracy tests (p fixed)\n");
    PRINTF2(fp,"----------------------------------------------------------------\n");
    PRINTF2(fp,"        N          p    error(FP64)   error(FP32)    T_quad(s)   T_num(s) \n");
    for (int count=0; count<COUNT_P; count++) {
        binoinv_v_err_compare(fp, N, p, Ulo_ex, Uhi_ex, Ulo_h, Uhi_h, ulo_ex, uhi_ex, ulo_h, uhi_h);
#ifdef PLUS
        N += 1.0f;
        if (N*p*(1.0f-p) > 8e+1f) break;
#else
        N *= 2.0f;
#endif
        if (N > Nmax) break;
    }

    // set values to test
    p = 1.0f;
    N = 1e3f;
    // change N,p based on input if desired
    if (argc > 1) {
        N = strtod(argv[1], NULL);
    }
    if (argc > 2) {
        p = strtod(argv[2], NULL);
    }
    PRINTF2(fp,"\nvector CPU algorithm accuracy tests (N fixed)\n");
    PRINTF2(fp,"----------------------------------------------------------------\n");
    PRINTF2(fp,"        N          p    error(FP64)   error(FP32)    T_quad(s)   T_num(s) \n");
    for (int count=0; count<COUNT_N; count++) {
        binoinv_v_err_compare(fp, N, p, Ulo_ex, Uhi_ex, Ulo_h, Uhi_h, ulo_ex, uhi_ex, ulo_h, uhi_h);
        p *= 0.5f;
    }

    // set values to test
    p = 0.5f;
    N = 1e3f;
    // change N,p based on input if desired
    if (argc > 1) {
        N = strtod(argv[1], NULL);
    }
    if (argc > 2) {
        p = strtod(argv[2], NULL);
    }
    PRINTF2(fp,"\nvector CPU algorithm accuracy tests (Np fixed)\n");
    PRINTF2(fp,"-----------------------------------------------------------------\n");
    PRINTF2(fp,"        N          p    error(FP64)   error(FP32)    T_quad(s)   T_num(s) \n");
    for (int count=0; count<COUNT_NP; count++) {
        binoinv_v_err_compare(fp, N, p, Ulo_ex, Uhi_ex, Ulo_h, Uhi_h, ulo_ex, uhi_ex, ulo_h, uhi_h);
        N *= 2.0f;
        p *= 0.5f;
        if (N > Nmax) break;
    }

#endif /* VECTOR */

    //////////////////////////////////////////////////
    // re-do for complementary version
    //////////////////////////////////////////////////

#ifdef COMPLEMENTARY

    PRINTF2(fp,"\n\n******************************\n");
    PRINTF2(fp,"Complementary version binocinv\n");
    PRINTF2(fp,"******************************\n");

    // set values to test
    p = 0.5f;
    N = 1.0f;
    // change N,p based on input if desired
    if (argc > 1) {
        N = strtod(argv[1], NULL);
    }
    if (argc > 2) {
        p = strtod(argv[2], NULL);
    }
    PRINTF2(fp,"\n");
    PRINTF2(fp,"\nscalar CPU algorithm accuracy tests (p fixed)\n");
    PRINTF2(fp,"----------------------------------------------------------------\n");
    PRINTF2(fp,"        N          p    error(FP64)   error(FP32)    T_quad(s)   T_num(s) \n");
    for (int count=0; count<COUNT_P; count++) {
        binocinv_err_compare(fp, N, p, Ulo_ex, Uhi_ex, Ulo_h, Uhi_h, ulo_ex, uhi_ex, ulo_h, uhi_h);
#ifdef PLUS
        N += 1.0f;
        if (N*p*(1.0f-p) > 8e+1f) break;
#else
        N *= 2.0f;
#endif
        if (N > Nmax) break;
    }

    // set values to test
    p = 1.0f;
    N = 1e3f;
    // change N,p based on input if desired
    if (argc > 1) {
        N = strtod(argv[1], NULL);
    }
    if (argc > 2) {
        p = strtod(argv[2], NULL);
    }

    PRINTF2(fp,"\nscalar CPU algorithm accuracy tests (N fixed)\n");
    PRINTF2(fp,"----------------------------------------------------------------\n");
    PRINTF2(fp,"        N          p    error(FP64)   error(FP32)    T_quad(s)   T_num(s) \n");
    for (int count=0; count<COUNT_N; count++) {
        binocinv_err_compare(fp, N, p, Ulo_ex, Uhi_ex, Ulo_h, Uhi_h, ulo_ex, uhi_ex, ulo_h, uhi_h);
        p *= 0.5;
    }

    // set values to test
    p = 0.5f;
    N = 1e3f;
    // change N,p based on input if desired
    if (argc > 1) {
        N = strtod(argv[1], NULL);
    }
    if (argc > 2) {
        p = strtod(argv[2], NULL);
    }
    PRINTF2(fp,"\nscalar CPU algorithm accuracy tests (Np fixed)\n");
    PRINTF2(fp,"----------------------------------------------------------------\n");
    PRINTF2(fp,"        N          p    error(FP64)   error(FP32)    T_quad(s)   T_num(s) \n");
    for (int count=0; count<COUNT_NP; count++) {
        binocinv_err_compare(fp, N, p, Ulo_ex, Uhi_ex, Ulo_h, Uhi_h, ulo_ex, uhi_ex, ulo_h, uhi_h);
        N *= 2.0f;
        p *= 0.5f;
    }
#endif /* COMPLEMENTARY */

    // free memory

    free(Ulo_h);
    free(Uhi_h);
    free(ulo_h);
    free(uhi_h);

    free(ulo_ex);
    free(uhi_ex);
    free(Ulo_ex);
    free(Uhi_ex);

    fclose(fp);
}
