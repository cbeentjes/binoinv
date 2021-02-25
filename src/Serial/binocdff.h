/////////////////////////////////////////////////////////////////
//                                                             //
// This software was written by Casper Beentjes & Mike Giles   //
// 2021                                                        //
//                                                             //
// It is copyright University of Oxford, and provided under    //
// the terms of the GNU GPLv3 license:                         //
// http://www.gnu.org/licenses/gpl.html                        //
//                                                             //
// Commercial users wanting to use the software under a more   //
// permissive license, such as BSD, should contact the authors://
// giles@maths.ox.ac.uk; beentjes@maths.ox.ac.uk               //
//                                                             //
/////////////////////////////////////////////////////////////////

#ifndef BINOCDFF_H
#define BINOCDFF_H

#include <stdio.h>

// Standard math header file, use Intel math library when compiling
// with the Intel C compiler
#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <math.h>
#endif

// Intel math library contains erfcxf implementation, but the
// standard C math library does not, so include an implementation
#ifndef __MATHIMF_H_INCLUDED
#include "erfcxf.h"
#endif

// Mathematical constants definition
#ifndef M_LN_2PIF
#define M_LN_2PIF	1.8378770f	/* log(2*pi) */
#endif

#ifndef M_SQRT_2PI_INVF
#define M_SQRT_2PI_INVF 0.3989423f /* 1/sqrt(2*pi) */
#endif

#ifndef M_2PIF
#define M_2PIF 6.2831855f /* 2*pi */
#endif

//////////////////////////////////////////////////////////////////////
//                                                                  //
// Binomial distribution (N,p) CDF at k                             //
// i.e. for K~Bin(N,p) calculate binom_cdf(k,N,p) = P(K<=k)         //
//                                                                  //
// Note that P(K<=k) = I_{1-p}(N-k,k+1)                             //
//                                                                  //
// The routine below calculates both                                //
// S0 = P(K<=k)                                                     //
// S1 = 1 - P(K<=k) = P(K>k)                                        //
// and contains an adapted C version of code in                     //
//                                                                  //
// ALGORITHM TOMS708: ACM Transactions on Mathematical Software     //
// (TOMS) (1992) VOL. 18, IS. 3, 360-373.                           //
// DOI: 10.1145/131766.131776                                       //
//                                                                  //
//////////////////////////////////////////////////////////////////////

// prototype declarations for incomplete beta function calculation
static void betaif(float x, float y, float a, float b, float apbm1,
                    float *S0, float *S1);
static float betacff(float x, float y, float a, float b, float lambda);
static float betaexpf(float a, float b, float lambda);
static float brcompf(float a, float b, float lambda);
static float rlog1f(float x);
static float bcorrf(float a, float b);

// prototype declarations for binomial pdf evaluation
static float binom_pmff(float k, float N, float p, float q);

void binom_cdff(float k, float N, float p, float q, float *S0, float *S1) {
    // NOTE: edge case check is not needed for binoinv project as CDF code
    // is never called for such instances.
#ifndef __FAST_CDF
    // Check edge cases
    if (!(N == truncf(N)) || N < 0.0f || isinf(N)) {
        *S0 = FP_NAN, *S1 = FP_NAN;
        return;
    }
    if (!(p >= 0.0f && q >= 0.0f)) {
        *S0 = FP_NAN, *S1 = FP_NAN;
        return;
    }
    if (isnan(k)) {
        *S0 = FP_NAN, *S1 = FP_NAN;
        return;
    }
    if (k <  0.0f) {
        *S0 = 0.0f, *S1 = 1.0f;
        return;
    }
    if (k >= N) {
        *S0 = 1.0f, *S1 = 0.0f;
        return;
    }
    if (p == 0.0f) {
        *S0 = 1.0f, *S1 = 0.0f;
        return;
    }
    if (q == 0.0f) {
        *S0 = ((k==N) ? 1.0f : 0.0f), *S1 = 1.0f - *S0;
        return;
    }
#endif /* __FAST_CDF */

    // If k <= (N+1)p/(2-p) or k >= (N+1)*2p/(1+p) - 1 it is advantageous to
    // evaluate the CDF by direct summation, otherwise resort to evaluation of
    // the incomplete beta function via a continued fraction or asymptotic
    // expansion based on TOMS708 work.
    // The latter approaches only work if k and N-k are big enough (> 20), if
    // k or N-k is small resort to direct summation as well
    if ( ( k*(2.0f - p) > (N + 1.0f)*p ) &&
         ( (k + 1.0f)*(1.0f + p) < (N + 1.0f)*2.0f*p) &&
         ( k > 20.0f && N - k > 20.0f ) ) {
#ifndef __FAST_CDF
        k = truncf(k);
#endif /* __FAST_CDF */
        // Calculate I_{1-p}(N-k,k+1) via incomplete beta function procedure
        betaif(q, p, N-k, k+1.0f, N, S0, S1);
    }
    // Direct summation definition of the CDF; single precision accuracy in 20
    // iterations as p_{k-1}/p_{k} < 1/2 or p_{k+1}/p_{k} < 1/2
    else {
        float S, T, a, n;
        int flag;

        S = 0.0f;

        // Use fast recursive calculation of binomial PMF if N < 2^24 as
        // N - n + 1 != N and n != n - 1 in single precision
        if (N < 16777216.0f) {
            flag = (k > N*p);
            if (flag) {
                k = N - (k + 1.0f);
                a=p; p=q; q=a;  // switch p and q
            }

            // Make sure to convert k to the correct integer value
#ifndef __FAST_CDF
            n = ( flag ? ceilf(k) : truncf(k) );
#else
            n = k;
#endif /* __FAST_CDF */

            T = binom_pmff(n, N, p, q);

            // Perform summation using scaled PMF (T) and CDF (S)
            // to prevent repeated use of div-operation
            S = T;
            a = 1.0f;
            for (int i = 0; i < 26; i++) {
                T *= (n*q);
                S  = (N-n+1.0f)*p*S + T;
                a *= (N-n+1.0f)*p;
                n -= 1.0f;
                // Prevent overflow of a_n, S_n, T_n
                if (a > 1.0e+20f) {
                    a *= 1.0e-10f;
                    S *= 1.0e-10f;
                    T *= 1.0e-10f;
                }
            }
            S /= a;

        }
        // Otherwise need to use direct calculation of binomial PMF because
        // n - 1 == n or N - n + 1 == N - n can hold in single precision
        else {

            // if k > N*0.5 then (N-k) < k so swapping k & (N-k)
            // improves the accuracy of updates of k
            flag = (k > N*0.5f);
            if (flag) {
                k = (N - k);
                a=p; p=q; q=a;  // switch p and q
            }

            // Need to ensure the summation is in the correct direction, i.e.
            // away from the mean of the distribution. If k > N*p sum upwards,
            // otherwise sum downwards
            int sgn = ( k > N*p ? 1: -1);

            // Make sure to convert k to the correct integer value
#ifndef __FAST_CDF
            n = truncf(k);
#else
            n = k;
#endif /* __FAST_CDF */
            // NOTE: the calculated S is the true upper tail if either
            // flag == true or k > N*p, but not both.
            flag = (flag) ^ (sgn > 0);

            // Correct starting point of summation if upper tail is calculated
            n += (flag * sgn);

            for (int i = 0; i < 20; i++) {
                T = binom_pmff(n + sgn*i, N, p, q);
                S += T;
            }
        }

        *S0 = S;
        *S1 = (0.5f - S) + 0.5f;
        if (flag) {
            a=*S0; *S0=*S1; *S1=a;
        }
    }
    return;
} /* binom_cdff */

// Regularized incomplete beta function betai(x,y,a,b,apbm1) = I_{x}(a,b)
// This code is based on BRATIO from TOMS708 but tailored to the specific
// application of the binomial CDF
// NOTE: when betai is called it is assumed x+y=1, 0<x<1 and a,b > 20
static void betaif(float x, float y, float a, float b, float apbm1,
                    float *S0, float *S1) {
    float lambda;

    // NOTE: take a smaller EPS than in TOMS708 to achieve full float
    // precision on the log-scale of the CDF value
#define EPS 1.0e-8f

    // NOTE: lambda = ((a > b) ? (a + b)*y - b : a - (a+b)*x )
    // However, to prevent loss of accuracy rewrite to use (a+b-1)
    // lambda = ((a>b) ? ((a+b-1)*y - b) + y : (a - (a+b-1)*x) - x )
    // because a+b-1 == N is available without cost and
    // cancellation error is reduced when b~(a+b)*y / a~(a+b)*x
    lambda = ((a > b) ? (apbm1*y - b) + y : (a - apbm1*x) - x );
    int flag = (lambda < 0.0f);
    if (flag) {
        *S0=b; b=a; a=*S0; // switch a and b
        *S0=x; x=y; y=*S0; // switch x and y
        lambda = fabsf(lambda);
    }

    // Call asymptotic expansion in central region for large a,b
    if (   (a > 100.0f && b >= a && lambda <= 0.03f*a )
        || (b > 100.0f && a >  b && lambda <= 0.03f*b ) ) {
        *S0 = betaexpf(a,b,lambda);
    // Call continued fraction otherwise
    } else {
        *S0 = betacff(x,y,a,b,lambda);
    }

    *S1 = 0.5f + (0.5f - *S0);
    if (flag) {
        lambda=*S0; *S0=*S1; *S1=lambda;
    }
    return;
} /* betaif */

// betaexpf - asymptotic expansion for the incomplete beta function I_x(a,b)
// This code is equal to BASYM from TOMS708
// NOTE: when betaexpf is called it is assumed in TOMS708 that a,b > 15 and
// lambda = (a+b)y - b with 0<y<1. Here, however, a>20 and b>20 is guaranteed.
static float betaexpf(float a, float b, float lambda) {

#define NUM 20
#define E0  1.1283792f  /* 2/sqrt(pi) */
#define E1  0.35355338f /* 2^(-3/2)   */


    float A0[NUM + 1], B0[NUM + 1], C[NUM + 1], D[NUM + 1];
    float f, u, t, j0, j1, r0, r1, w0, w, h, h2, hn, sum, s;
    float t0, t1, r, z0, z, z2, znm1, zn, bsum, dsum;

    float eps = 100.0f * EPS;

    // NOTE: due to the fact that the betai routine is only called when
    // (N+1)p * 1/(2-p) + 1 < k + 1 < (N+1)p * 2/(1+p) the value of
    // ea = -lambda/a and eb = lambda/b satisfy the following bounds
    // -max(p,1-p)/2 <  ea <= 0
    //        0      <= eb < max(p,1-p)
    // so in general -1/2 < e < 1 for either of the two and rlog1 is tailored
    // to this range of input.
    f = a*rlog1f(-lambda/a) + b*rlog1f(lambda/b);
    t = expf(-f);
    if (t == 0.0f) return 0.0f;

    z0 = sqrtf(f);
    z  = 0.5f*(z0/E1);
    z2 = f + f;

    if (a < b) {
        h = a/b;
        w0 = 1.0f / sqrtf(a*(1.0f + h));
        r1 = (b - a) / b;
    } else {
        h = b/a;
        w0 = 1.0f / sqrtf(b*(1.0f + h));
        r1 = (b - a) / a;
    }
    r0 = 1.0f / (1.0f + h);

    A0[0] = (2.0f/3.0f)*r1;
    C[0]  = -0.5f*A0[0];
    D[0]  = -C[0];
    j0 = (0.5f/E0)*erfcxf(z0);
    j1 = E1;
    sum = j0 + D[0]*w0*j1;

    s  = 1.0f;
    h2 = h*h;
    hn = 1.0f;
    w  = w0;
    znm1 = z;
    zn = z2;

    int mmj;
    for (int n=2; n <= NUM; n+=2) {
        hn *= h2;
        A0[n-1] = r0*2.0f*(h*hn + 1.0f)/(n + 2.0f);
        s += hn;
        A0[n]   = r1*2.0f*s/(n + 3.0f);

        for (int i=n; i<=n+1; i++) {
            r = -0.5f*(i + 1);
            B0[0] = r*A0[0];
            for (int m=2; m<=i; m++) {
                bsum = 0.0f;
                for (int j=1; j<=m-1; j++) {
                    mmj = m - j;
                    bsum += (j*r - mmj)*A0[j-1]*B0[mmj-1];
                }
                B0[m-1] = r*A0[m-1] + bsum/m;
            }
            C[i-1] = B0[i-1]/(i+1);

            dsum = 0.0f;
            for (int j=1; j<=i-1; j++) {
                dsum += D[i-j-1]*C[j-1];
            }
            D[i-1] = -(dsum + C[i-1]);
        }

        j0 = E1*znm1 + (n-1)*j0;
        j1 = E1*zn   + n*j1;
        znm1 *= z2;
        zn   *= z2;

        w *= w0;
        t0 = D[n-1]*w*j0;
        w *= w0;
        t1 = D[n]*w*j1;
        sum += (t0 + t1);
        if (fabsf(t0) + fabsf(t1) <= eps*sum) {
            break;
        }
    }

    u = expf(-bcorrf(a,b));
    return E0*t*u*sum;
} /* betaexpf */

// betacff - continued fraction expansion for the inc. beta function I_x(a,b)
// This code is equal to BFRAC from TOMS708
// NOTE: when betacff is called it is assumed in TOMS708 that a,b>1 and
// lambda = (a+b)y - b with 0<y<1. Here, however, a>20 and b>20 is guaranteed.
static float betacff(float x, float y, float a, float b, float lambda) {

#define MAXIT 1000

    float n, c, c0, c1, yp1, p, s, t, e, w, bt;
    float an0, an1, bn0, bn1, alpha, beta, r0, r1;

    float eps = 15.0f*EPS;

    bt = brcompf(a,b,lambda); // limiting step for accuracy
    if (bt <= 0.0f) return 0.0f;

	c  = 1.0f + lambda;
    c0 = b/a;
    c1 = 1.0f + 1.0f/a;
    yp1 = y + 1.0f;

    n = 0.0f;
    p = 1.0f;
    s = a + 1.0f;

    an0 = 0.0f;
    bn0 = 1.0f;
    an1 = 1.0f;
    bn1 = c/c1;
    r1  = c1/c;

	while (n <= MAXIT) {
        n += 1.0f;
        t = n/a;
        w = n*(b-n)*x;
        e = a/s;
        alpha = (p*(p + c0)*e*e)*(w*x);
        e  = (1.0f + t)/(c1 + t + t);
        beta = n + w/s + e*(c + n*yp1);
        p  = 1.0f + t;
        s += 2.0f;

        t = alpha*an0 + beta*an1;
        an0 = an1; an1 = t;
        t = alpha*bn0 + beta*bn1;
        bn0 = bn1; bn1 = t;

        r0 = r1;
        r1 = an1/bn1;

        // NOTE: if n == b then r1 == r0 holds
        if (fabsf(r1 - r0) < eps*r1) {
            return bt*r1;
        }

        an0 = an0/bn1;
        bn0 = bn0/bn1;
        an1 = r1;
        bn1 = 1.0f;
	}
	if (n > MAXIT) printf("Warning: reached max iterations in continued "
                          "fraction expansion for I_x(a,b).\n");
	return bt*r1;
} /* betacff */

// brcompf(a,b,lambda) = x^a*y^b / Beta(a,b)
// This code is based on BRCOMP from TOMS708 but tailored to the specific
// application of the binomial CDF
// NOTE: when brcompf is called it is assumed in TOMS708 that a,b>8 and
// lambda = (a+b)*y - b with 0<y<1. Here a>20 and b>20 is guaranteed.
static float brcompf(float a, float b, float lambda) {
    float h, x0, e, z, u, v;
    // This step is still used to circumvent overflow issues when calculating
    // a*b/(a+b) if a & b are both very large
    if (a > b) {
        h  = a/b;
        x0 = h  /(1.0f + h);
    }
    else {
        h  = b/a;
        x0 = 1.0f/(1.0f + h);
    }

    // NOTE: due to the fact that the betaif routine is only called when
    // (N+1)p * 1/(2-p) + 1 < k + 1 < (N+1)p * 2/(1+p) the value of
    // ea = -lambda/a and eb = lambda/b satisfy the following bounds
    // -max(p,1-p)/2 <  ea <= 0
    //        0      <= eb < max(p,1-p)
    // so in general -1/2 < e < 1 for either of the two and rlog1 is tailored
    // to this range of input.
    e = -lambda/a;
    u = rlog1f(e);

    e = lambda/b;
    v = rlog1f(e);

    z = expf(-(a*u + b*v)); // this step is limiting for accuracy

    return M_SQRT_2PI_INVF*sqrtf(b*x0)*z*expf(-bcorrf(a,b));
} /* brcompf */

// rlog1f(x) = x - log(1+x) for -0.618<x<1.618 using rational minimax approx.
// This code is inspired by RLOG1 from TOMS708 but is significantly more
// accurate and is tailored to the range of inputs that occur for the binomial
// CDF evaluation
static float rlog1f(float x) {

#define P0  0.19999999f
#define P1 -0.11665086f
#define Q1 -1.2975472f
#define Q2  0.37142160f

    float r, t, w;

    // rlog1f(x) = f(r) = 2*r^2*(1/(1-r) + r*f1(r))
    // where f1(r) = (log((1+r)/(1-r)) - 2*r)/(2*r^3)
    // which has series expansion sum_{n=0} r^(2n)/(2n + 3)
    // Calculate f1(r) = 1/3 + r^2*P(r^2)/Q(r^2) where P,Q follow
    // from rational minimax approximation for 0 <= r^2 <= 0.2
    r = x/(x + 2.0f);
    t = r*r;
    w = (1.0f/3.0f) + t*(P1*t + P0) /
        ((Q2*t + Q1)*t + 1.0f);
    return t*((x + 2.0f) - 2.0f*r*w);
} /* rlog1f */

// bcorrf(a, b) = D(a) + D(b) - D(a+b) using polynomial minimax approx.
// where D(a) = log(Gamma(a)) - (a-0.5)log(a) + a - 0.5log(2*pi)
// This code is inspired by BCORR from TOMS708 but is slightly faster because
// it is tailored to the range of inputs that occur for the binomial
// CDF evaluation. This guarantees that a,b > 20. The minimax approximation
// for D(a),D(b) and D(a+b) is derived for 20<a,b<80.
static float bcorrf(float a, float b) {

#define C0  0.083333333f  /* 1/12 */
#define C1 -0.0027767253f /* ~-1/360 */

    float h, c, S3, t, x, w;

    if (a > b) {
        h=b; b=a; a=h; //Switch a and b
    }

    h  = a/b;
    c  = h   /(1.0f + h);
    x  = 1.0f/(1.0f + h);

    S3 = 1.0f + (x + x*x);

    h = 1.0f/b;
    t = h*h;

    w = C1*S3*t + C0;
    w *= (c/b);

    h = 1.0f/a;
    t = h*h;
    return (C1*t + C0 )/a + w;
} /* bcorrf */

#undef MAXIT
#undef EPS
#undef NUM
#undef E0
#undef E1
#undef C0
#undef C1
#undef P0
#undef P1
#undef Q1
#undef Q2

//////////////////////////////////////////////////////////////////////
//                                                                  //
// Binomial distribution (N,p) p.m.f. at k                          //
// i.e. for K~Bin(N,p) calculate binom_pmf(k,N,p) = P(K=k)          //
//                                                                  //
// The routine below is an adapted C version of the code in         //
//                                                                  //
// Fast and Accurate Computation of Binomial Probabilities          //
// by Catherine Loader (2000)                                       //
//                                                                  //
//////////////////////////////////////////////////////////////////////


// Lookup table for Stirling approximation error for N<=15
static float sfef[16] = {
    0.0f        ,  0.081061467f,
    0.041340696f,  0.027677926f,
    0.020790672f,  0.016644691f,
    0.013876129f,  0.011896710f,
    0.010411265f,  0.0092554622f,
    0.0083305634f, 0.0075736755f,
    0.0069428401f, 0.0064089942f,
    0.0059513701f, 0.0055547336f
};

// stirlerrf(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
static float stirlerrf(float N) {

#define S0 0.083333333f     /* 1/12 */
#define S1 0.0027777778f    /* 1/360 */

    if (N < 16.0f)  return sfef[(int) N];
    if (N > 500.0f) return S0/N;
    return (S0-S1/(N*N))/N;
} /* stirlerrf */

// bd0f(x,np) = x log(x/np) + np - x, deviance term in binomial pmf
// NOTE: this code is based on bd0 in Loader(2000), but adapted to improve
// precision for small values of P(K=k)
static float bd0f(float k, float NP) {

#define P0  0.19999997f
#define P1 -0.11809043f
#define Q1 -1.3047534f
#define Q2  0.37667380f

    float t, s, w, v;

    // Rational minimax expansion for small x/NP
    // NOTE: this approximation is valid for NP/3 < k < 3*NP
    if (fabsf(k-NP) < 0.5f*(k+NP)) {
        v = (k-NP)/(k+NP); // |v| < 0.5
        s = (k-NP)*v;
        t = v*v;
        w = (1.0f/3.0f) + t*(P1*t + P0) /
                          ((Q2*t + Q1)*t + 1.0f);
        s += (2.0f*t*v*w)*k;
        return s;
    }
    // Direct evaluation (limiting step for accuracy)
    return k*logf(k/NP)+NP-k;
} /* bd0f */

// binom_pmff - Evaluate the binomial p.m.f. P(K=k)
static float binom_pmff(float k, float N, float p, float q) {
    float lc;
    // Check edge cases
    if (!(N == truncf(N)) || N < 0.0f || isinf(N)) return FP_NAN;
    if (!(p >= 0.0f && q >= 0.0f)) return FP_NAN;
    if (isnan(k))  return FP_NAN;
    if (p == 0.0f) return ((k == 0.0f) ? 1.0f : 0.0f);
    if (q == 0.0f) return ((k == N)    ? 1.0f : 0.0f);
    if (k == 0.0f) return expf(N*log1pf(-p));
    if (k == N)    return expf(N*logf(p));
    if (k < 0.0f || k > N || k != truncf(k) ) return 0.0f;
    // NOTE: calculating bd0 is the limiting step for accuracy for small
    // values of P(K=k), i.e. large values of log(P(K=k)), when
    // |v| = |k - NP|/(k + NP) is larger than 0.5
    lc = stirlerrf(N) - stirlerrf(k) - stirlerrf(N-k)
            - bd0f(k,N*p) - bd0f(N-k,N*q);

    // NOTE: this is not needed for binoinv project as code is never called
    // for such edge cases.
#ifndef __FAST_CDF
    float lf;
    // NOTE: sqrt(N/(2*PI*k*(N-k))) can cause overflow/underflow
    // errors, so if N >> 1 use logarithm transformation.
    // Here N > 2^63 means N >> 1 as k*(N-k) <= N^2/4
    if ( N > 9.2233720e18f ) {
        if (k < 0.6f*N) {
            lf = M_LN_2PIF + logf(k) + log1pf(-k/N);
        } else {
            lf = M_LN_2PIF + logf(k) + logf(N-k) - logf(N);
        }
        return expf(lc - 0.5f*lf);
    }
#endif /* __FAST_CDF */
    return expf(lc)*sqrtf(N/(M_2PIF*k*(N-k)));

} /* binom_pmff */

#undef S0
#undef S1
#undef P0
#undef P1
#undef Q1
#undef Q2

#endif /* BINOCDFF_H */
