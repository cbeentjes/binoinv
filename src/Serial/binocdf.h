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

#ifndef BINOCDF_H
#define BINOCDF_H

#include <stdio.h>

// Standard math header file, use Intel math library when compiling
// with the Intel C compiler
#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <math.h>
#endif

// Intel math library contains erfcx implementation, but the
// standard C math library does not, so include an implementation
#ifndef __MATHIMF_H_INCLUDED
#include "erfcx.h"
#endif

// Mathematical constants definition
#ifndef M_LN_2PI
#define M_LN_2PI	1.837877066409345483560659472811	/* log(2*pi) */
#endif

#ifndef M_SQRT_2PI_INV
#define M_SQRT_2PI_INV 0.39894228040143267794 /* 1/sqrt(2*pi) */
#endif

#ifndef M_2PI
#define M_2PI 6.2831853071795865 /* 2*pi */
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
static void betai(double x, double y, double a, double b, double apbm1,
                    double *S0, double *S1);
static double betacf(double x, double y, double a, double b, double lambda);
static double betaexp(double a, double b, double lambda);
static double brcomp(double a, double b, double lambda);
static double rlog1(double x);
static double bcorr(double a, double b);

// prototype declarations for binomial pdf evaluation
static double binom_pmf(double k, double N, double p, double q);

void binom_cdf(double k, double N, double p, double q, double *S0, double *S1) {
    // NOTE: edge case check is not needed for binoinv project as CDF code
    // is never called for such instances.
#ifndef __FAST_CDF
    // Check edge cases
    if (!(N == trunc(N)) || N < 0.0 || isinf(N)) {
        *S0 = FP_NAN, *S1 = FP_NAN;
        return;
    }
    if (!(p >= 0.0 && q >= 0.0)) {
        *S0 = FP_NAN, *S1 = FP_NAN;
        return;
    }
    if (isnan(k)) {
        *S0 = FP_NAN, *S1 = FP_NAN;
        return;
    }
    if (k <  0.0) {
        *S0 = 0.0, *S1 = 1.0;
        return;
    }
    if (k >= N) {
        *S0 = 1.0, *S1 = 0.0;
        return;
    }
    if (p == 0.0) {
        *S0 = 1.0, *S1 = 0.0;
        return;
    }
    if (q == 0.0) {
        *S0 = ((k==N) ? 1.0 : 0.0), *S1 = 1.0 - *S0;
        return;
    }
#endif /* __FAST_CDF */

    // If k <= (N+1)p/(2-p) or k >= (N+1)*2p/(1+p) - 1 it is advantageous to
    // evaluate the CDF by direct summation, otherwise resort to evaluation of
    // the incomplete beta function via a continued fraction or asymptotic
    // expansion based on TOMS708 work.
    // The latter approaches only work if k and N-k are big enough (> 20), if
    // k or N-k is small resort to direct summation as well
    if ( ( k*(2.0 - p) > (N + 1.0)*p ) &&
         ( (k + 1.0)*(1.0 + p) < (N + 1.0)*2.0*p) &&
         ( k > 20.0 && N - k > 20.0 ) ) {
#ifndef __FAST_CDF
        k = trunc(k);
#endif /* __FAST_CDF */
        // Calculate I_{1-p}(N-k,k+1) via incomplete beta function procedure
        betai(q, p, N-k, k+1.0, N, S0, S1);
    }
    // Direct summation definition of the CDF; double precision accuracy in 50
    // iterations as p_{k-1}/p_{k} < 1/2 or p_{k+1}/p_{k} < 1/2
    else {
        double S, T, a, n;
        int flag;

        S = 0.0;

        // Use fast recursive calculation of binomial PMF if N < 2^53 as
        // N - n + 1 != N and n != n - 1 in double precision
        if (N < 9007199254740992.0) {
            flag = (k > N*p);
            if (flag) {
                k = N - (k + 1.0);
                a=p; p=q; q=a;  // switch p and q
            }

            // Make sure to convert k to the correct integer value
#ifndef __FAST_CDF
            n = ( flag ? ceil(k) : trunc(k) );
#else
            n = k;
#endif /* __FAST_CDF */

            T = binom_pmf(n, N, p, q);

            // Perform summation using scaled PMF (T) and CDF (S)
            // to prevent repeated use of div-operation
            S = T;
            a = 1.0;
            for (int i = 0; i < 50; i++) {
                T *= (n*q);
                S  = (N-n+1.0)*p*S + T;
                a *= (N-n+1.0)*p;
                n -= 1.0;
                // Prevent overflow of a_n, S_n, T_n
                if (a > 1.0e270) {
                    T *= 1.0e-220;
                    a *= 1.0e-220;
                    S *= 1.0e-220;
                }
            }
            S /= a;

        }
        // Otherwise need to use direct calculation of binomial PMF because
        // n - 1 == n or N - n + 1 == N - n can hold in double precision
        else {

            // if k > N*0.5 then (N-k) < k so swapping k & (N-k)
            // improves the accuracy of updates of k
            flag = (k > N*0.5);
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
            n = trunc(k);
#else
            n = k;
#endif /* __FAST_CDF */

            // NOTE: the calculated S is the true upper tail if either
            // flag == true or k > N*p, but not both.
            flag = (flag) ^ (sgn > 0);

            // Correct starting point of summation if upper tail is calculated
            n += (flag * sgn);

            for (int i = 0; i < 50; i++) {
                T = binom_pmf(n + sgn*i, N, p, q);
                S += T;
            }
        }

        *S0 = S;
        *S1 = (0.5 - S) + 0.5;
        if (flag) {
            a=*S0; *S0=*S1; *S1=a;
        }
    }
    return;
} /* binom_cdf */

// Regularized incomplete beta function betai(x,y,a,b,apbm1) = I_{x}(a,b)
// This code is based on BRATIO from TOMS708 but tailored to the specific
// application of the binomial CDF
// NOTE: when betai is called it is assumed x+y=1, 0<x<1 and a,b > 20
static void betai(double x, double y, double a, double b, double apbm1,
                    double *S0, double *S1) {
    double lambda;

    // NOTE: take a smaller EPS than in TOMS708 to achieve full double
    // precision on the log-scale of the CDF value
#define EPS 1.0e-17

    // NOTE: lambda = ((a > b) ? (a + b)*y - b : a - (a+b)*x )
    // However, to prevent loss of accuracy rewrite to use (a+b-1)
    // lambda = ((a>b) ? ((a+b-1)*y - b) + y : (a - (a+b-1)*x) - x )
    // because a+b-1 == N is available without cost and
    // cancellation error is reduced when b~(a+b)*y / a~(a+b)*x
    lambda = ((a > b) ? (apbm1*y - b) + y : (a - apbm1*x) - x );
    int flag = (lambda < 0.0);
    if (flag) {
        *S0=b; b=a; a=*S0; // switch a and b
        *S0=x; x=y; y=*S0; // switch x and y
        lambda = fabs(lambda);
    }

    // Call asymptotic expansion in central region for large a,b
    if (   (a > 100.0 && b >= a && lambda <= 0.03*a )
        || (b > 100.0 && a >  b && lambda <= 0.03*b ) ) {
        *S0 = betaexp(a,b,lambda);
    // Call continued fraction otherwise
    } else {
        *S0 = betacf(x,y,a,b,lambda);
    }

    *S1 = 0.5 + (0.5 - *S0);
    if (flag) {
        lambda=*S0; *S0=*S1; *S1=lambda;
    }
    return;
} /* betai */

// betaexp - asymptotic expansion for the incomplete beta function I_x(a,b)
// This code is equal to BASYM from TOMS708
// NOTE: when betaexp is called it is assumed in TOMS708 that a,b > 15 and
// lambda = (a+b)y - b with 0<y<1. Here, however, a>20 and b>20 is guaranteed.
static double betaexp(double a, double b, double lambda) {

#define NUM 20
#define E0  1.12837916709551  /* 2/sqrt(pi) */
#define E1  0.353553390593274 /* 2^(-3/2)   */


    double A0[NUM + 1], B0[NUM + 1], C[NUM + 1], D[NUM + 1];
    double f, u, t, j0, j1, r0, r1, w0, w, h, h2, hn, sum, s;
    double t0, t1, r, z0, z, z2, znm1, zn, bsum, dsum;

    double eps = 100.0 * EPS;

    // NOTE: due to the fact that the betai routine is only called when
    // (N+1)p * 1/(2-p) + 1 < k + 1 < (N+1)p * 2/(1+p) the value of
    // ea = -lambda/a and eb = lambda/b satisfy the following bounds
    // -max(p,1-p)/2 <  ea <= 0
    //        0      <= eb < max(p,1-p)
    // so in general -1/2 < e < 1 for either of the two and rlog1 is tailored
    // to this range of input.
    f = a*rlog1(-lambda/a) + b*rlog1(lambda/b);
    t = exp(-f);
    if (t == 0.0) return 0.0;

    z0 = sqrt(f);
    z  = 0.5*(z0/E1);
    z2 = f + f;

    if (a < b) {
        h = a/b;
        w0 = 1.0 / sqrt(a*(1.0 + h));
        r1 = (b - a) / b;
    } else {
        h = b/a;
        w0 = 1.0 / sqrt(b*(1.0 + h));
        r1 = (b - a) / a;
    }
    r0 = 1.0 / (1.0 + h);

    A0[0] = (2.0/3.0)*r1;
    C[0]  = -0.5*A0[0];
    D[0]  = -C[0];
    j0 = (0.5/E0)*erfcx(z0);
    j1 = E1;
    sum = j0 + D[0]*w0*j1;

    s  = 1.0;
    h2 = h*h;
    hn = 1.0;
    w  = w0;
    znm1 = z;
    zn = z2;

    int mmj;
    for (int n=2; n <= NUM; n+=2) {
        hn *= h2;
        A0[n-1] = r0*2.0*(h*hn + 1.0)/(n + 2.0);
        s += hn;
        A0[n]   = r1*2.0*s/(n + 3.0);

        for (int i=n; i<=n+1; i++) {
            r = -0.5*(i + 1);
            B0[0] = r*A0[0];
            for (int m=2; m<=i; m++) {
                bsum = 0.0;
                for (int j=1; j<=m-1; j++) {
                    mmj = m - j;
                    bsum += (j*r - mmj)*A0[j-1]*B0[mmj-1];
                }
                B0[m-1] = r*A0[m-1] + bsum/m;
            }
            C[i-1] = B0[i-1]/(i+1);

            dsum = 0.0;
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
        if (fabs(t0) + fabs(t1) <= eps*sum) {
            break;
        }
    }

    u = exp(-bcorr(a,b));
    return E0*t*u*sum;
} /* betaexp */

// betacf - continued fraction expansion for the inc. beta function I_x(a,b)
// This code is equal to BFRAC from TOMS708
// NOTE: when betacf is called it is assumed in TOMS708 that a,b>1 and
// lambda = (a+b)y - b with 0<y<1. Here, however, a>20 and b>20 is guaranteed.
static double betacf(double x, double y, double a, double b, double lambda) {

#define MAXIT 1000

    double n, c, c0, c1, yp1, p, s, t, e, w, bt;
    double an0, an1, bn0, bn1, alpha, beta, r0, r1;

    double eps = 15.0*EPS;

    bt = brcomp(a,b,lambda); // limiting step for accuracy
    if (bt <= 0.0) return 0.0;

	c  = 1.0 + lambda;
    c0 = b/a;
    c1 = 1.0 + 1.0/a;
    yp1 = y + 1.0;

    n = 0.0;
    p = 1.0;
    s = a + 1.0;

    an0 = 0.0;
    bn0 = 1.0;
    an1 = 1.0;
    bn1 = c/c1;
    r1  = c1/c;

	while (n <= MAXIT) {
        n += 1.0;
        t = n/a;
        w = n*(b-n)*x;
        e = a/s;
        alpha = (p*(p + c0)*e*e)*(w*x);
        e  = (1.0 + t)/(c1 + t + t);
        beta = n + w/s + e*(c + n*yp1);
        p  = 1.0 + t;
        s += 2.0;

        t = alpha*an0 + beta*an1;
        an0 = an1; an1 = t;
        t = alpha*bn0 + beta*bn1;
        bn0 = bn1; bn1 = t;

        r0 = r1;
        r1 = an1/bn1;

        // NOTE: if n == b then r1 == r0 holds
        if (fabs(r1 - r0) < eps*r1) {
            return bt*r1;
        }

        an0 = an0/bn1;
        bn0 = bn0/bn1;
        an1 = r1;
        bn1 = 1.0;
	}
	if (n > MAXIT) printf("Warning: reached max iterations in continued "
                          "fraction expansion for I_x(a,b).\n");
	return bt*r1;
} /* betacf */

// brcomp(a,b,lambda) = x^a*y^b / Beta(a,b)
// This code is based on BRCOMP from TOMS708 but tailored to the specific
// application of the binomial CDF
// NOTE: when brcomp is called it is assumed in TOMS708 that a,b>8 and
// lambda = (a+b)*y - b with 0<y<1. Here a>20 and b>20 is guaranteed.
static double brcomp(double a, double b, double lambda) {
    double h, x0, e, z, u, v;
    // This step is still used to circumvent overflow issues when calculating
    // a*b/(a+b) if a & b are both very large
    if (a > b) {
        h  = a/b;
        x0 = h  /(1.0 + h);
    }
    else {
        h  = b/a;
        x0 = 1.0/(1.0 + h);
    }

    // NOTE: due to the fact that the betai routine is only called when
    // (N+1)p * 1/(2-p) + 1 < k + 1 < (N+1)p * 2/(1+p) the value of
    // ea = -lambda/a and eb = lambda/b satisfy the following bounds
    // -max(p,1-p)/2 <  ea <= 0
    //        0      <= eb < max(p,1-p)
    // so in general -1/2 < e < 1 for either of the two and rlog1 is tailored
    // to this range of input.
    e = -lambda/a;
    u = rlog1(e);

    e = lambda/b;
    v = rlog1(e);

    z = exp(-(a*u + b*v)); // this step is limiting for accuracy

    return M_SQRT_2PI_INV*sqrt(b*x0)*z*exp(-bcorr(a,b));
} /* brcomp */

// rlog1(x) = x - log(1+x) for -0.618<x<1.618 using rational minimax approx.
// This code is inspired by RLOG1 from TOMS708 but is significantly more
// accurate and is tailored to the range of inputs that occur for the binomial
// CDF evaluation
static double rlog1(double x) {

#define P0  0.2000000000000000
#define P1 -0.3636535967811319
#define P2  0.2015244511825799
#define P3 -0.03274937605228191
#define P4  0.00004542775258423288
#define Q1 -2.532553698191349
#define Q2  2.261033627633934
#define Q3 -0.8263420776370621
#define Q4  0.1008870710149048

    double r, t, w;

    // rlog1(x) = f(r) = 2*r^2*(1/(1-r) + r*f1(r))
    // where f1(r) = (log((1+r)/(1-r)) - 2*r)/(2*r^3)
    // which has series expansion sum_{n=0} r^(2n)/(2n + 3)
    // Calculate f1(r) = 1/3 + r^2*P(r^2)/Q(r^2) where P,Q follow
    // from rational minimax approximation for 0 <= r^2 <= 0.2
    r = x/(x + 2.0);
    t = r*r;
    w = (1.0/3.0) + t*((((P4*t + P3)*t + P2)*t + P1)*t + P0 ) /
                      ((((Q4*t + Q3)*t + Q2)*t + Q1)*t + 1.0);
    return t*((x + 2.0) - 2.0*r*w);
} /* rlog1 */

// bcorr(a, b) = D(a) + D(b) - D(a+b) using polynomial minimax approx.
// where D(a) = log(Gamma(a)) - (a-0.5)log(a) + a - 0.5log(2*pi)
// This code is inspired by BCORR from TOMS708 but is slightly faster because
// it is tailored to the range of inputs that occur for the binomial
// CDF evaluation. This guarantees that a,b > 20. The minimax approximation
// for D(a),D(b) and D(a+b) is derived for 20<a,b<80.
static double bcorr(double a, double b) {

#define C0  0.083333333333333333    /* 1/12 */
#define C1 -0.0027777777777760393   /* ~-1/360 */
#define C2  0.00079365078287927436  /* ~1/1260 */
#define C3 -0.00059522072260950274  /* ~-1/1680 */
#define C4  0.00083168805235713158  /* ~1/1188 */

    double h, c, S3, S5, S7, S9, t, x, x2, w;

    if (a > b) {
        h=b; b=a; a=h; //Switch a and b
    }

    h  = a/b;
    c  = h  /(1.0 + h);
    x  = 1.0/(1.0 + h);
    x2 = x*x;

    S3 = 1.0 + (x + x2);
    S5 = 1.0 + (x + x2*S3);
    S7 = 1.0 + (x + x2*S5);
    S9 = 1.0 + (x + x2*S7);

    h = 1.0/b;
    t = h*h;

    w = (((C4*S9*t + C3*S7)*t + C2*S5)*t + C1*S3)*t + C0;
    w *= (c/b);

    h = 1.0/a;
    t = h*h;
    return ((((C4*t + C3)*t + C2)*t + C1)*t + C0 )/a + w;

} /* bcorr */

#undef MAXIT
#undef EPS
#undef NUM
#undef E0
#undef E1
#undef C0
#undef C1
#undef C2
#undef C3
#undef C4
#undef C5
#undef P0
#undef P1
#undef P2
#undef P3
#undef P4
#undef Q1
#undef Q2
#undef Q3
#undef Q4

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
static double sfe[16] = {
    0.0                          , 0.081061466795327258219670264,
    0.041340695955409294093822081, 0.0276779256849983391487892927,
    0.020790672103765093111522771, 0.0166446911898211921631948653,
    0.013876128823070747998745727, 0.0118967099458917700950557241,
    0.010411265261972096497478567, 0.0092554621827127329177286366,
    0.008330563433362871256469318, 0.0075736754879518407949720242,
    0.006942840107209529865664152, 0.0064089941880042070684396310,
    0.005951370112758847735624416, 0.0055547335519628013710386899
};

// stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
static double stirlerr(double N) {

#define S0 0.083333333333333333333          /* 1/12 */
#define S1 0.00277777777777777777778        /* 1/360 */
#define S2 0.00079365079365079365079365     /* 1/1260 */
#define S3 0.000595238095238095238095238    /* 1/1680 */
#define S4 0.0008417508417508417508417508   /* 1/1188 */

    if (N < 16.0)  return sfe[(int) N];
    double nn = N*N;
    if (N > 500.0) return (S0-S1/nn)/N;
    if (N > 80.0)  return (S0-(S1-S2/nn)/nn)/N;
    if (N > 35.0)  return (S0-(S1-(S2-S3/nn)/nn)/nn)/N;
    return (S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/N;
} /* stirlerr */

// bd0(x,np) = x log(x/np) + np - x, deviance term in binomial pmf
// NOTE: this code is based on bd0 in Loader(2000), but adapted to improve
// precision for small values of P(K=k)
static double bd0(double k, double NP) {

#define P0  0.2000000000000000
#define P1 -0.3669101082562791
#define P2  0.2056846557845579
#define P3 -0.03394955490410171
#define P4  0.00005007384892489343
#define Q1 -2.548836255566920
#define Q2  2.293465048754533
#define Q3 -0.8464624069944372
#define Q4  0.1046656636182324

    double t, s, w, v;

    // Rational minimax expansion for small x/NP
    // NOTE: this approximation is valid for NP/3 < k < 3*NP
    if (fabs(k-NP) < 0.5*(k+NP)) {
        v = (k-NP)/(k+NP); // |v| < 0.5
        s = (k-NP)*v;
        t = v*v;
        w = (1.0/3.0) + t*((((P4*t + P3)*t + P2)*t + P1)*t + P0) /
                          ((((Q4*t + Q3)*t + Q2)*t + Q1)*t + 1.0);
        s += (2.0*t*v*w)*k;
        return s;
    }
    // Direct evaluation (limiting step for accuracy)
    return k*log(k/NP)+NP-k;
} /* bd0 */

// binom_pmf - Evaluate the binomial p.m.f. P(K=k)
static double binom_pmf(double k, double N, double p, double q) {
    double lc;
    // Check edge cases
    if (!(N == trunc(N)) || N < 0.0 || isinf(N)) return FP_NAN;
    if (!(p >= 0.0 && q >= 0.0)) return FP_NAN;
    if (isnan(k)) return FP_NAN;
    if (p == 0.0) return ((k == 0.0) ? 1.0 : 0.0);
    if (q == 0.0) return ((k == N)   ? 1.0 : 0.0);
    if (k == 0.0) return exp(N*log1p(-p));
    if (k == N)   return exp(N*log(p));
    if (k < 0.0 || k > N || k != trunc(k) ) return 0.0;
    // NOTE: calculating bd0 is the limiting step for accuracy for small
    // values of P(K=k), i.e. large values of log(P(K=k)), when
    // |v| = |k - NP|/(k + NP) is larger than 0.5
    lc = stirlerr(N) - stirlerr(k) - stirlerr(N-k)
            - bd0(k,N*p) - bd0(N-k,N*q);

    // NOTE: this is not needed for binoinv project as code is never called
    // for such edge cases.
#ifndef __FAST_CDF
    double lf;
    // NOTE: sqrt(N/(2*PI*k*(N-k))) can cause overflow/underflow
    // errors, so if N >> 1 use logarithm transformation.
    // Here N > 2^511 means N >> 1 as k*(N-k) <= N^2/4
    if ( N > 6.703903964971299e+153 ) {
        if (k < 0.6*N) {
            lf = M_LN_2PI + log(k) + log1p(-k/N);
        } else {
            lf = M_LN_2PI + log(k) + log(N-k) - log(N);
        }
        return exp(lc - 0.5*lf);
    }
#endif /* __FAST_CDF */
    return exp(lc)*sqrt(N/(M_2PI*k*(N-k)));

} /* binom_pmf */

#undef S0
#undef S1
#undef S2
#undef S3
#undef S4
#undef P0
#undef P1
#undef P2
#undef P3
#undef P4
#undef Q1
#undef Q2
#undef Q3
#undef Q4

#endif /* BINOCDF_H */
