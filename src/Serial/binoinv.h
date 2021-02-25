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

#ifndef BINOINV_H
#define BINOINV_H

// Standard Math Library or Intel Math library
#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <math.h>
#endif

// declare prototype for inverse Normal CDF function
// defined at the bottom of this header file

static double normcdfinv_as241(double);
#define norminv(u) normcdfinv_as241(u);

// helper header file containing binomial CDF evaluation
// use __FAST_CDF macro to ignore edge case checks in binomial CDF routines
#ifndef BINOCDF_H
#define __FAST_CDF
#include "binocdf.h"
#undef __FAST_CDF
#endif

// As described in the TOMS paper, there are two versions;
// the first  (binoinv) is optimised for MIMD execution, whereas
// the second (binoinv_v)  is designed for SIMD/vector execution

// binoinv(u, N, p):
//
// This double precision function computes the inverse
// of the binomial CDF
//
// u   = CDF value in range [0,1]
// N   = Number of trials
// p   = Success probability for each trial
//
//  max |error| no more than 1 for N*p < ~9e+15
//  ave |error| < 1.29e-15*sqrt(max(1 ,N*p*(1-p)))  for 0 <= N*p*(1-p) <  7
//              < 6.04e-16*sqrt(max(32,N*p*(1-p)))  for 7 <= N*p*(1-p) < ~9e+15
//
//  For N*p >= ~9e+15, the errors will be about 1 ulp.
//

static double binoinv_core(double U, double V, double N, double p, double q);
static double binoinv_summation(double U, double V, double N, double p, double q);

double binoinv(double U, double N, double p) {
    return binoinv_core(U, 1.0-U, N, p, 1.0-p);
}

double binocinv(double V, double N, double p) {
    return binoinv_core(1.0-V, V, N, p, 1.0-p);
}

// Forcing inline here improves performance typically by ~1% to 5%
inline __attribute__((always_inline)) static double binoinv_core(double U, double V, double N, double p, double q) {

    double X=0.0, NP, NQ, NPQ_root, NPQ_rooti, W, Wsqr, S, Del;

// handle exceptions -- constants defined in <math.h>

    // handles NAN inputs as well
    if (!(N == trunc(N)) || N < 0.0 || isinf(N)) return FP_NAN;
    if (!(U >= 0.0 && V >= 0.0)) return FP_NAN;
    if (!(p >= 0.0 && q >= 0.0)) return FP_NAN;
    if (U == 0.0 || p == 0.0) return 0.0;
    if (V == 0.0 || q == 0.0) return N;

    // To maintain accuracy even for large N work with min(p,q)
    int flag = (p > 0.5);
    // Swap p & q and U & V.
    if (flag) {
        W=p; p=q; q=W;
        W=U; U=V; V=W;
    }

    NP = N*p; NQ = N*q;

#define mu_0 10.0

// If Np or Nq is small ( <= mu_0), i.e. p is close to 0 or 1 we perform a
// direct summation. The number of loop iterations is then bounded
// by ~ 30 + Np + 17*sqrt(Np) for Np small,
// by ~180 + Nq + 46*sqrt(Nq) for Nq small.
//
// For large Np and Nq (i.e. > mu_0) we attempt to use asymptotic
// approximations rather than direct summation.

    if (NP > mu_0 && NQ > mu_0) {

        W = norminv(fmin(U,V));
        if (U > V) W = -W;

        NPQ_root  = sqrt(NP*q);
        NPQ_rooti = 1.0/NPQ_root;

// use polynomial approximations in central region

        if (fabs(W) < 3.0) {

            Wsqr = W*W;
            // Q_N2
            S = (NPQ_root*W) +
                ( 0.5 + (q - p)*(Wsqr - 1.0)*(1.0/6.0) ) +
                ( (7.0*p*q - 1.0) - (p*q + 0.5)*Wsqr )*W*NPQ_rooti*(1.0/36.0);

            // Get Del from Q_N3 difference with Q_N2 (for p != 1/2)
            // and from Q_N4 difference with Q_N2 (for p == 1/2)
            Del = 0.0025;
            Del = 0.0050 + Del*Wsqr;
            Del = 0.0100 + Del*Wsqr;
            Del = Del * (0.25*NPQ_rooti + fabs(q - p));
            Del = Del * (1.0 + p) * (1.0 + q);
            Del = Del * (NPQ_rooti*NPQ_rooti);

            // Sum from smallest to largest to minimise rounding error;
            S += Del;
        }

// otherwise use Temme uniform approximation

        else {
            double nu, eta0, eta1, eta, sqpq, pr, a0, a1, a2, a3, a4, a5,
                   xi, xi_min_p, xi_eta;
            sqpq = sqrt(p*q);
            nu   = N + 1.0;

            eta0 = -W / sqrt(nu);

            // References refer to 2020 paper by Gil, Segura & Temme
            // "ASYMPTOTIC INVERSION OF THE BINOMIAL AND NEGATIVE BINOMIAL
            // CUMULATIVE DISTRIBUTION FUNCTIONS"
            // DOI: 10.1553/etna_vol52s270

            // use asymptotic expansions to solve for xi when eta near 0
            if (fabs(eta0) < 1.25e-3*sqpq) {
                eta0 = eta0 / sqpq;
                pr = p - 0.5;

                // Taylor approximation to xi - p

                // coeffs on page 273 below (2.7)
                a2 = pr*(1.0/3.0);
                a3 =    (-0.75   + pr*pr)*(1.0/36.0);
                a4 = pr*( 2.25   - pr*pr)*(1.0/270.0);
                a5 =    (-2.4375 + pr*pr*(-13.5 + pr*pr))*(1.0/4320.0);
                // expansion on page 273 in (2.7)
                xi_min_p = a5;
                xi_min_p = a4   + eta0*xi_min_p;
                xi_min_p = a3   + eta0*xi_min_p;
                xi_min_p = a2   + eta0*xi_min_p;
                xi_min_p = 1.0  + eta0*xi_min_p;
                xi_min_p *= -(p - p*p)*eta0;

                xi = p + xi_min_p;

                // Taylor approximation to eta1 * xi_eta
                // NOTE 1: use a direct expansion of the product, rather than
                // using separate expansions for eta1 and xi_eta.
                // NOTE 2: only compute xi_eta*eta1 to third order because
                // it is not multiplied by nu in S = nu*xi + eta1*xi_eta
                a0 = pr*(-2.0/3.0);
                a1 =    (  3.75  - 11.0*pr*pr)*(1.0/36.0);
                a2 = pr*(-15.75  +  7.0*pr*pr)*(1.0/810.0);
                a3 =    ( 30.375 + pr*pr*(-2.25 + pr*pr*71.0))*(1.0/9720.0);

                eta1 = a3;
                eta1 = a2 + eta0*eta1;
                eta1 = a1 + eta0*eta1;
                eta1 = a0 + eta0*eta1;

                xi_eta = 1.0; // xi_eta is already incorporated in eta1
            }
            // use Newton iteration to solve for xi in (2.3) on page 272
            // when eta not near 0
            else {
                xi = p;
                eta = 0.0;
                xi_eta = -sqpq;

                for (int k = 0; k < 5; k++) {
                    xi -= (eta - eta0)*xi_eta;
                    xi = fmax(fmin(xi, 1.0 - 1.0/nu), 1.0/nu);   // clamp at two ends
                    eta = sqrt(-2.0*(xi*log1p((p - xi)/xi) + (1.0-xi)*log1p(-(p - xi)/(1.0 - xi))));
                    if (p < xi) eta = -eta;
                    xi_eta = - eta / log1p((p - xi)/(q*xi));   // log1p used for accuracy
                }

                eta1 = log(sqrt(xi - xi*xi)*eta0/(p - xi)) / eta0;   // see (3.6)

                xi_min_p = xi - p;
            }

            Del = 2.25e-2/(nu*(xi - xi*xi));

            // Sum from smallest to largest to minimise rounding error;
            S = ((p + Del) + eta1*xi_eta) + xi_min_p*nu;
        }


        // NOTE: use of round down mode in final sum is important for
        // large values of NP to ensure fast execution
        // If NP >> 1 then Del can be smaller than ULP and the standard
        // IEEE_754 rounding mode (round to nearest, ties to even)
        // causes issues when the asymptotic approximation returns an
        // answer very close to an integer value. Notably it results in a
        // significant slowdown as the CDF must be checked increasingly
        // often. Round down mode for the final step of the summation in
        // the asymptotic approximation (somewhat) mitigates this.
        // Only apply the round down transformation if needed as this
        // operation is (relatively) costly on CPUs and not often needed.

        // Round down to nearest integer and check accuracy

        if (NP < 1e+15 || NP > 1e+16) {
            X = trunc(S+NP);
        } else {
            X = S;
            X += NP;
            if ((X - NP) >= S) {       // If X += NP was rounded up
                X = nextafter(X,0.0);  // use round down mode instead
            }
            X = trunc(X);
        }


        if (10.0 <= X && X <= N - 10.0) {
            // If NO_CDF is defined we skip this step, this increases speed
            // but can significantly impact accuracy
#ifdef NO_CDF
            return (flag) ? N - X : X;
#endif

            // Accurate enough, no need to check the CDF
            if (S - 2.0*Del >= X - NP || X > 9.007199254740992e+15) {
                return (flag) ? N - X : X;
            }
            // Not accurate enough, need to check the CDF
            //
            // Calculate both C0 = F(X-1) and C1 = 1 - F(X-1) and
            // check both U against C0 and V against C1 (if necessary) to
            // prevent loss of accuracy in the tails
            double C0, C1;
            X -= 1.0;
            binom_cdf(X, N, p, q, &C0, &C1);
            if ( ( U <= 0.5 && C0 < U ) || (V < 0.5 && V < C1 ) ) {
                X += 1.0;
            }
            return (flag) ? N - X : X;
        }
        // At this point it is guaranteed that either X < 10 or X > N - 10.
        // Decide if p & q and U & V are swapped depending on whether direct
        // summation will be used when in the right tail of the distribution.
        if (U > 0.5) {
            flag = !flag;
            W=p; p=q; q=W;
            W=U; U=V; V=W;
        }
    }

// Use direct summation routine either if Np or N(1-p) is small,
// or U is too close to 0 or 1 for Temme approximation to be valid
    X = binoinv_summation(U, V, N, p, q);
    return (flag) ? N - X : X;
}

// Fast direct summation routine
// ICC automatically inlines at -O3, gcc does not and inline makes
// execution slightly faster for small values of N
inline static double binoinv_summation(double U, double V, double N,
        double p, double q) {

    double n, m, a, S, T, d;

    n = 0.0;
    m = N + 1.0;

    T = 1.0;
    // Calculate a0 = 1/P(n==0) = (1-p)^(-N)
    // NOTE 1: Calculation of P(n==0) = (1-p)^N via iterative multiplication
    // for moderate N is fastest and safe accuracy-wise. For Np or Nq small
    // the standard pow-function is safe to use, even for N large, maintaining
    // accuracy whilst being faster than iterative calculation.
    // NOTE 2: For larger N when both Np and Nq are large the exponential of
    // logarithm approach is needed due to the possibility of overflow of a0.
    // Use of exponential of logarithm does mean accuracy is lost due to error
    // magnification, but this approach is only used in the tails of the
    // distribution and the overall error impact is thus small.
    //
    // The value of N for which iterative multiplaction is faster depends on
    // compiler.
#ifdef __INTEL_COMPILER
#define N_SWITCH 20.0
#else
#define N_SWITCH 27.0
#endif
    // Direct iterative calculation for small N
    if (N <= N_SWITCH) {
        a = 1.0;
        for (int k = 0; k < N; k++) a *= q;
        a = 1.0/a;
    } else {
        // NP or NQ small enough (note p,q swap), use pow-function
        if (N*p < 200.0) {
            a = pow(q, -N);
        // NP and NQ large, need to prevent overflow
        } else {
            a = -N*log1p(-p);
            // Include a factor to prevent (1-p)^(-N) from overflow
            // Let a0 = exp(-128)*(1-p)^(-N) & T = exp(-128)
            a = exp(a - 128.0);
            T = 2.572209372642415e-56; /* exp(-128.0) */
        }
    }

    d = 0.0;
    if (U > 0.5) d = 1.0e-13;

    S = T - a*(U-d);

    // bottom-up summation
    // NOTE: For N > 2^53 there is an extra slight loss of accuracy due to
    // incremental updates m -= 1.0 being lost as N - 1 == N. Effect on
    // ultimate outcome, however, is minimal.
    while (S < 0.0 && m > 1.0) {
        n += 1.0;
        m -= 1.0;
        T  = m*p*T;
        S  = n*q*S + T;
        a  = n*q*a;
    }

    // top-down summation if needed
    if (S < 2.0*d*a) {
        // Important to first multiply V*a to keep accuracy when V << 1
        while (m > 1.0 && T*1.0e17 > (V*a)) {
            n += 1.0;
            m -= 1.0;
            T  = m*p*T;
            a  = n*q*a;
            // Scale down a_n & T_n if a_n is close to overflow
            if (a > 1.0e270) {
                // NOTE: T > 1e-17*2^-1074*a > 1e-71 so no need
                // need to check if T is large enough to rescale
                T *= 1.0e-220;
                a *= 1.0e-220;
            }
        }
        S = T - a*V;

        while (S < 0.0) {
            T  = n*q*T;
            S  = m*p*S + T;
            n -= 1.0;
            m += 1.0;
        }
    }
    return n;
}

#undef mu_0
#undef N_SWITCH

// binoinv_v(u, N, p):
//
// This double precision function computes the inverse
// of the binomial CDF using an approach more suitable for SIMD execution
//
// u   = CDF value in range [0,1]
// N   = Number of trials
// p   = Success probability for each trial
//
//  max |error| no more than 1 for N*p < ~9e+15
//  ave |error| < 1.29e-15*sqrt(max(1 ,N*p*(1-p)))  for 0 <= N*p*(1-p) <  7
//              < 6.04e-16*sqrt(max(32,N*p*(1-p)))  for 7 <= N*p*(1-p) < ~9e+15
//
//  For N*p >= ~9e+15, the errors will be about 1 ulp.
//

double binoinv_v(double U, double N, double p) {

    double V = 1.0 - U;
    double q = 1.0 - p;

    double X=0.0, NP, NQ, W, S, Del;

// handle exceptions -- constants defined in <math.h>

    // handles NAN inputs as well
    if (!(N == trunc(N)) || N < 0.0 || isinf(N)) return FP_NAN;
    if (!(U >= 0.0 && V >= 0.0)) return FP_NAN;
    if (!(p >= 0.0 && q >= 0.0)) return FP_NAN;
    if (U == 0.0 || p == 0.0) return 0.0;
    if (V == 0.0 || q == 0.0) return N;

    // To maintain accuracy for large N work with min(p,q)
    int flag = (p > 0.5);
    W = U; // copy U to W to be used in argument for norminv
    // Swap p & q and U & V.
    if (flag) {
        W=p; p=q; q=W;
        W=U; U=V; V=W;
    }

    NP = N*p; NQ = N*q;

#define mu_0 10.0

// If Np or Nq is small ( <= mu_0), i.e. p is close to 0 or 1 we perform a
// direct summation. The number of loop iterations is then bounded
// by ~ 30 + Np + 17*sqrt(Np) for Np small,
// by ~180 + Nq + 46*sqrt(Nq) for Nq small.
//
// For large Np and Nq (i.e. > mu_0) we attempt to use asymptotic
// approximations rather than direct summation.

    if (NP > mu_0 && NQ > mu_0) {

        W = norminv(W);
        if (flag) W = -W;

// use Temme uniform approximation

        double nu, eta0, eta1, eta, sqpq, pr, a0, a1, a2, a3, a4, a5, a6, a7,
               a8, a9, xi, xi_min_p, xi_eta;
        sqpq = sqrt(p*q);
        nu   = N + 1.0;

        eta0 = -W / sqrt(nu);

        // References refer to 2020 paper by Gil, Segura & Temme
        // "ASYMPTOTIC INVERSION OF THE BINOMIAL AND NEGATIVE BINOMIAL
        // CUMULATIVE DISTRIBUTION FUNCTIONS"
        // DOI: 10.1553/etna_vol52s270

        Del = 0.0;

        // use asymptotic expansions to solve for xi when eta near 0
        if (fabs(eta0) < 5.48e-1*sqpq) {
            eta0 = eta0 / sqpq;
            pr = p - 0.5;

            // Taylor approximation to xi - p

            // first coeffs from page 273 below (2.7)
            a2 = pr*(1.0/3.0);
            a3 =    (-0.75   + pr*pr)*(1.0/36.0);
            a4 = pr*( 2.25   - pr*pr)*(1.0/270.0);
            a5 =    (-2.4375 + pr*pr*(-13.5 + pr*pr))*(1.0/4320.0);
            a6 = pr*(14.0625 + pr*pr*( 19.5 + pr*pr))*(1.0/17010.0);
            a7 =   -(13095.0 + pr*pr*(  270756.0 + pr*pr*(143568.0
                            + 8896.0*pr*pr)))*(1.0/348364800.0);
            a8 = pr*( 1287.0 + pr*pr*(    7692.0 + pr*pr*(1872.0
                            + 64.0*pr*pr)))*(1.0/13063680.0);
            a9 =   -(2049867.0 + pr*pr*(90331632.0 + pr*pr*(234206496.0
                            + pr*pr*(28628736.0 + pr*pr*146176.0))))*(1.0/601974374400.0);
            // expansion on page 273 in (2.7)
            xi_min_p = a9;
            xi_min_p = a8   + eta0*xi_min_p;
            xi_min_p = a7   + eta0*xi_min_p;
            xi_min_p = a6   + eta0*xi_min_p;
            xi_min_p = a5   + eta0*xi_min_p;
            xi_min_p = a4   + eta0*xi_min_p;
            xi_min_p = a3   + eta0*xi_min_p;
            xi_min_p = a2   + eta0*xi_min_p;
            xi_min_p = 1.0  + eta0*xi_min_p;
            xi_min_p *= -(p - p*p)*eta0;

            xi = p + xi_min_p;

            // Taylor approximation to eta1 * xi_eta
            // NOTE 1: use a direct expansion of the product, rather than
            // using separate expansions for eta1 and xi_eta.
            // NOTE 2: only compute xi_eta*eta1 to fifth order because
            // it is not multiplied by nu in S = nu*xi + eta1*xi_eta
            a0 = pr*(-2.0/3.0);
            a1 =    (  3.75  - 11.0*pr*pr)*(1.0/36.0);
            a2 = pr*(-15.75  +  7.0*pr*pr)*(1.0/810.0);
            a3 =    ( 30.375 + pr*pr*(-2.25 + pr*pr*71.0))*(1.0/9720.0);
            a4 = pr*(-58.375 + pr*pr*(43.0  + pr*pr*38.0))*(1.0/15120.0);
            a5 =    (104574.375 + pr*pr*(1104205.5 + pr*pr*(-822006 + pr*pr*147688)))*(1/391910400);
            xi_eta = 1.0;

            eta1 = a5;
            eta1 = a4 + eta0*eta1;
            eta1 = a3 + eta0*eta1;
            eta1 = a2 + eta0*eta1;
            eta1 = a1 + eta0*eta1;
            eta1 = a0 + eta0*eta1;

            if (fabs(eta0) > 2.74e-0) {
                Del = 4.54e-4;
            }
        }
        // use Newton iteration to solve for xi in (2.3) on page 272
        // when eta not near 0
        else {
            xi = p;
            eta = 0.0;
            xi_eta = -sqpq;

            for (int k = 0; k < 5; k++) {
              xi -= (eta - eta0)*xi_eta;
              xi = fmax(fmin(xi, 1.0 - 1.0/nu), 1.0/nu);   // clamp at two ends
              eta = sqrt(-2.0*(xi*log1p((p - xi)/xi) + (1.0-xi)*log1p(-(p - xi)/(1.0 - xi))));
              if (p < xi) eta = -eta;
              xi_eta = - eta / log1p((p - xi)/(q*xi));   // log1p used for accuracy
            }

            eta1 = log(sqrt(xi - xi*xi)*eta0/(p - xi)) / eta0;   // see (3.6)

            xi_min_p = xi - p;
        }

        Del = fmax(Del, 2.25e-2*(1.0/(nu*(xi - xi*xi))));

        // Sum from smallest to largest to minimise rounding error;
        S = ((p + Del) + eta1*xi_eta) + xi_min_p*nu;

        // NOTE: use of round down mode in final sum is important for
        // large values of NP to ensure fast execution
        // If NP >> 1 then Del can be smaller than ULP and the standard
        // IEEE_754 rounding mode (round to nearest, ties to even)
        // causes issues when the asymptotic approximation returns an
        // answer very close to an integer value. Notably it results in a
        // significant slowdown as the CDF must be checked increasingly
        // often. Round down mode for the final step of the summation in
        // the asymptotic approximation (somewhat) mitigates this.
        // On GPU/SIMD devices simply always use this instruction.

        // Round down to nearest integer and check accuracy if 10 <= X <= N-10

        // This emulates the use of round-down mode, rather than relying on
        // the specific instruction
        X = S;
        X += NP;
        if ((X - NP) >= S) {       // If X += NP was rounded up
            X = nextafter(X,0.0);  // use round down mode instead
        }
        X = trunc(X);

        if (10.0 <= X && X <= N - 10.0) {
            // If NO_CDF is defined we skip this step, this increases speed
            // but can significantly impact accuracy
#ifdef NO_CDF
            return (flag) ? N - X : X;
#endif

            // Accurate enough, no need to check the CDF
            if (S - 2.0*Del >= X - NP || X > 9.007199254740992e+15 ) {
                return (flag) ? N - X : X;
            }
            // Not accurate enough, need to check the CDF
            //
            // Calculate both C0 = F(X-1) and C1 = 1 - F(X-1) and
            // check both U against C0 and V against C1 (if necessary) to
            // prevent loss of accuracy in the tails
            double C0, C1;
            X -= 1.0;
            binom_cdf(X, N, p, q, &C0, &C1);
            if ( ( U <= 0.5 && C0 < U ) || (V < 0.5 && V < C1 ) ) {
                X += 1.0;
            }
            return (flag) ? N - X : X;
        }
        // At this point know that either X < 10 or X > N - 10.
        // Decide if p & q and U & V are swapped depending on whether direct
        // summation will be used when in the right tail of the distribution.
        if (U > 0.5) {
            flag = !flag;
            W=p; p=q; q=W;
            W=U; U=V; V=W;
        }
    }

// Use direct summation routine either if Np or N(1-p) is small,
// or U is too close to 0 or 1 for Temme approximation to be valid
    X = binoinv_summation(U, V, N, p, q);
    return (flag) ? N - X : X;
}

#undef mu_0

//////////////////////////////////////////////////////////////////////
//                                                                  //
// The routine below is a C version of the code in                  //
//                                                                  //
// ALGORITHM AS241: APPLIED STATS (1988) VOL. 37, NO. 3, 477-44.    //
// http://lib.stat.cmu.edu/apstat/241                               //
//                                                                  //
// The relative error is less than 1e-15, and the accuracy is       //
// verified in the accompanying MATLAB code as241_test.m            //
//                                                                  //
//////////////////////////////////////////////////////////////////////

static double normcdfinv_as241(double p) {

  double q, r, num, den, res;

  q = p - 0.5;
  if (fabs(q) <= 0.425) {
    r = 0.180625 - q*q;

    num =         2.5090809287301226727e+3;
    num = r*num + 3.3430575583588128105e+4;
    num = r*num + 6.7265770927008700853e+4;
    num = r*num + 4.5921953931549871457e+4;
    num = r*num + 1.3731693765509461125e+4;
    num = r*num + 1.9715909503065514427e+3;
    num = r*num + 1.3314166789178437745e+2;
    num = r*num + 3.3871328727963666080e0;

    den =         5.2264952788528545610e+3;
    den = r*den + 2.8729085735721942674e+4;
    den = r*den + 3.9307895800092710610e+4;
    den = r*den + 2.1213794301586595867e+4;
    den = r*den + 5.3941960214247511077e+3;
    den = r*den + 6.8718700749205790830e+2;
    den = r*den + 4.2313330701600911252e+1;
    den = r*den + 1.0000000000e+00;

    res = q * num / den;

    return res;
  }

  else {

    if (q < 0.0)
      r = p;
    else
      r = 1.0 - p;

    r = sqrt(-log(r));

    if (r <= 5.0) {
      r = r - 1.6;

      num =         7.74545014278341407640e-4;
      num = r*num + 2.27238449892691845833e-2;
      num = r*num + 2.41780725177450611770e-1;
      num = r*num + 1.27045825245236838258e0;
      num = r*num + 3.64784832476320460504e0;
      num = r*num + 5.76949722146069140550e0;
      num = r*num + 4.63033784615654529590e0;
      num = r*num + 1.42343711074968357734e0;

      den =         1.05075007164441684324e-9;
      den = r*den + 5.47593808499534494600e-4;
      den = r*den + 1.51986665636164571966e-2;
      den = r*den + 1.48103976427480074590e-1;
      den = r*den + 6.89767334985100004550e-1;
      den = r*den + 1.67638483018380384940e0;
      den = r*den + 2.05319162663775882187e0;
      den = r*den + 1.0000000000e+00;

      res = num / den;
    }

    else {
      r = r - 5.0;

      num =         2.01033439929228813265e-7;
      num = r*num + 2.71155556874348757815e-5;
      num = r*num + 1.24266094738807843860e-3;
      num = r*num + 2.65321895265761230930e-2;
      num = r*num + 2.96560571828504891230e-1;
      num = r*num + 1.78482653991729133580e0;
      num = r*num + 5.46378491116411436990e0;
      num = r*num + 6.65790464350110377720e0;

      den =         2.04426310338993978564e-15;
      den = r*den + 1.42151175831644588870e-7;
      den = r*den + 1.84631831751005468180e-5;
      den = r*den + 7.86869131145613259100e-4;
      den = r*den + 1.48753612908506148525e-2;
      den = r*den + 1.36929880922735805310e-1;
      den = r*den + 5.99832206555887937690e-1;
      den = r*den + 1.0000000000e+00;

      res = num / den;
    }

    if (q < 0.0)
      res = - res;

    return res;
  }
}

#endif /* BINOINV_H */
