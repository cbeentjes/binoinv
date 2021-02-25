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

#ifndef BINOINVF_H
#define BINOINVF_H

// Standard Math Library or Intel Math library
#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <math.h>
#endif

// declare prototype for inverse Normal CDF function
// defined at the bottom of this header file

static float normcdfinvf_as241(float);
#define norminvf(u) normcdfinvf_as241(u);

// helper header file containing binomial CDF evaluation
// use __FAST_CDF macro to ignore edge case checks in binomial CDF routines
#ifndef BINOCDFF_H
#define __FAST_CDF
#include "binocdff.h"
#undef __FAST_CDF
#endif

// As described in the TOMS paper, there are two versions;
// the first  (binoinvf) is optimised for MIMD execution, whereas
// the second (binoinvf_v)  is designed for SIMD/vector execution

// binoinvf(u, N, p):
//
// This single precision function computes the inverse
// of the binomial CDF
//
// u   = CDF value in range [0,1]
// N   = Number of trials
// p   = Success probability for each trial
//
//  max |error| no more than 1 for N*p < ~1e+07
//  ave |error| < 5.14e-07*sqrt(max(0.5,N*p*(1-p))) for 0 <= N*p*(1-p) <  8
//              < 2.57e-07*sqrt(max(32 ,N*p*(1-p))) for 8 <= N*p*(1-p) < ~1e+07
//
//  For N*p >= ~1e+07, the errors will be about 1 ulp.
//

static float binoinvf_core(float U, float V, float N, float p, float q);
static float binoinvf_summation(float U, float V, float N, float p, float q);

float binoinvf(float U, float N, float p) {
    return binoinvf_core(U, 1.0f-U, N, p, 1.0f-p);
}

float binocinvf(float U, float N, float p) {
    return binoinvf_core(1.0f-U, U, N, p, 1.0f-p);
}

// Forcing inline here improves performance typically by ~1% to 5%
inline __attribute__((always_inline)) float binoinvf_core(float U, float V, float N, float p, float q) {

    float X=0.0f, NP, NQ, NPQ_root, NPQ_rooti, W, Wsqr, S, Del;

    // handle exceptions -- constants defined in <math.h>

    // handles NAN inputs as well
    if (!(N == truncf(N)) || N < 0.0f || isinf(N)) return FP_NAN;
    if (!(U >= 0.0f && V >= 0.0f)) return FP_NAN;
    if (!(p >= 0.0f && q >= 0.0f)) return FP_NAN;
    if (U == 0.0f || p == 0.0f) return 0.0f;
    if (V == 0.0f || q == 0.0f) return N;

    // To maintain accuracy even for large N work with min(p,q)
    int flag = (p > 0.5f);
    // Swap p & q and U & V.
    if (flag) {
        W=p; p=q; q=W;
        W=U; U=V; V=W;
    }

    NP = N*p; NQ = N*q;

#define mu_0 10.0f

    // If Np or Nq is small ( <= mu_0), i.e. p is close to 0 or 1 we perform a
    // direct summation. The number of loop iterations is then bounded
    // by ~ 16 + Np + 11*sqrt(Np) for Np small,
    // by ~ 30 + Nq + 20*sqrt(Nq) for Nq small.
    //
    // For large Np and Nq (i.e. > mu_0) we attempt to use asymptotic
    // approximations rather than direct summation.

    if (NP > mu_0 && NQ > mu_0) {

        W = norminvf(fminf(U,V));
        if (U > V) W = -W;

        NPQ_root  = sqrtf(NP*q);
        NPQ_rooti = 1.0f/NPQ_root;

        // use polynomial approximations in central region

        if (fabsf(W) < 3.0f) {

            Wsqr = W*W;
            // Q_N2
            S = (NPQ_root*W) +
                ( 0.5f + (q - p)*(Wsqr - 1.0f)*(1.0f/6.0f) ) +
                ( (7.0f*p*q - 1.0f) - (p*q + 0.5f)*Wsqr )*W*NPQ_rooti*(1.0f/36.0f);

            // Get Del from Q_N3 difference with Q_N2 (for p != 1/2)
            // and from Q_N4 difference with Q_N2 (for p == 1/2)
            Del = 0.0025f;
            Del = 0.0050f + Del*Wsqr;
            Del = 0.0100f + Del*Wsqr;
            Del = Del * (0.25f*NPQ_rooti + fabsf(q - p));
            Del = Del * (1.0f + p) * (1.0f + q);
            Del = Del * (NPQ_rooti*NPQ_rooti);

            // Sum from smallest to largest to minimise rounding error;
            S += Del;
        }

        // otherwise use Temme uniform approximation

        else {
            float nu, eta0, eta1, eta, sqpq, pr, a0, a1, a2, a3,
                  xi, xi_min_p, xi_eta;
            sqpq = sqrtf(p*q);
            nu   = N + 1.0f;

            eta0 = -W / sqrtf(nu);

            // References refer to 2020 paper by Gil, Segura & Temme
            // "ASYMPTOTIC INVERSION OF THE BINOMIAL AND NEGATIVE BINOMIAL
            // CUMULATIVE DISTRIBUTION FUNCTIONS"
            // DOI: 10.1553/etna_vol52s270

            // use asymptotic expansions to solve for xi when eta near 0
            if (fabsf(eta0) < 6e-2f*sqpq) {
                eta0 = eta0 / sqpq;
                pr = p - 0.5f;

                // Taylor approximation to xi - p

                // coeffs on page 273 below (2.7)
                a2 = pr*(1.0f/3.0f);
                a3 =    (-0.75f   + pr*pr)*(1.0f/36.0f);
                // expansion on page 273 in (2.7)
                xi_min_p = a3;
                xi_min_p = a2   + eta0*xi_min_p;
                xi_min_p = 1.0f + eta0*xi_min_p;
                xi_min_p *= -(p - p*p)*eta0;

                xi = p + xi_min_p;

                // Taylor approximation to eta1 * xi_eta
                // NOTE 1: use a direct expansion of the product, rather than
                // using separate expansions for eta1 and xi_eta.
                // NOTE 2: only compute xi_eta*eta1 to first order because
                // it is not multiplied by nu in S = nu*xi + eta1*xi_eta
                a0 = pr*(-2.0f/3.0f);
                a1 =    ( 3.75f   - 11.0f*pr*pr)*(1.0f/36.0f);

                eta1 = a1;
                eta1 = a0 + eta0*eta1;

                xi_eta = 1.0f; // xi_eta is already incorporated in eta1
            }
            // use Newton iteration to solve for xi in (2.3) on page 272
            // when eta not near 0
            else {
                xi = p;
                eta = 0.0f;
                xi_eta = -sqpq;

                for (int k = 0; k < 5; k++) {
                    xi -= (eta - eta0)*xi_eta;
                    xi = fmaxf(fminf(xi, 1.0f - 1.0f/nu), 1.0f/nu);   // clamp at two ends
                    eta = sqrtf(-2.0f*(xi*log1pf((p - xi)/xi) + (1.0f-xi)*log1pf(-(p - xi)/(1.0f - xi))));
                    if (p < xi) eta = -eta;
                    xi_eta = - eta / log1pf((p - xi)/(q*xi));   // log1p used for accuracy
                }

                eta1 = logf(sqrtf(xi - xi*xi)*eta0/(p - xi)) / eta0;   // see (3.6)

                xi_min_p = xi - p;
            }

            Del = 2.25e-2f/(nu*(xi - xi*xi));

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
        // the asymptotic approximation (somewhat) mitigates this, whilst
        // the loss of accuracy is limited.
        // Only apply the round down transformation if needed as this
        // operation is (relatively) costly on CPUs and not often needed.

        // Round down to nearest integer and check accuracy

        if (NP < 5e+4f || NP > 2e+7f) {
            X = truncf(S+NP);
        } else {
            X = S;
            X += NP;
            if ((X - NP) >= S) {         // If X += NP was rounded up
                X = nextafterf(X,0.0f);  // use round down mode instead
            }
            X = truncf(X);
        }

        if (10.0f <= X && X <= N - 10.0f) {
            // If NO_CDF is defined we skip this step, this increases speed
            // but can significantly impact accuracy
#ifdef NO_CDF
            return (flag) ? N - X : X;
#endif

            // Accurate enough, no need to check the CDF
            if (S - 2.0f*Del >= X - NP || X > 1.6777216e+7f) {
                return (flag) ? N - X : X;
            }
            // Not accurate enough, need to check the CDF
            //
            // Calculate both C0 = F(X-1) and C1 = 1 - F(X-1) and
            // check both U against C0 and V against C1 (if necessary) to
            // prevent loss of accuracy in the tails
            float C0, C1;
            X -= 1.0f;
            binom_cdff(X, N, p, q, &C0, &C1);
            if ( ( U <= 0.5f && C0 < U ) || (V < 0.5f && V < C1 ) ) {
                X += 1.0f;
            }
            return (flag) ? N - X : X;
        }

        // At this point it is guaranteed that either X < 10 or X > N - 10.
        // Decide if p & q and U & V are swapped depending on whether direct
        // summation will be used when in the right tail of the distribution.
        if (U > 0.5f) {
            flag = !flag;
            W=p; p=q; q=W;
            W=U; U=V; V=W;
        }
    }

    // Use direct summation routine either if Np or N(1-p) is small,
    // or U is too close to 0 or 1 for Temme approximation to be valid
    X = binoinvf_summation(U, V, N, p, q);
    return (flag) ? N - X : X;
}

// Fast direct summation routine
// ICC automatically inlines at -O3, gcc does not and inline makes
// execution slightly faster for small values of N
// NOTE: this summation routine implicitly relies on the fact that both
// NP < 14 and NQ < 14 hold. Otherwise more checks for overflow/underflow
// will need to be carried out.
inline float binoinvf_summation(float U, float V, float N,
        float p, float q) {

    float n, m, a, S, T, d;

    n = 0.0f;
    m = N + 1.0f;

    T = 1.0f;
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
#define N_SWITCH 16.0f
#else
#define N_SWITCH 17.0f
#endif
    // Direct iterative calculation for small N
    if (N <= N_SWITCH) {
        a = 1.0f;
        for (int k = 0; k < N; k++) a *= q;
        a = 1.0f/a;
    } else {
        // NP or NQ small enough (note p,q swap), use pow-function
        if (N*p < 40.0f) {
            a = powf(q, -N);
        // NP and NQ large, need to prevent overflow
        } else {
            a = -N*log1pf(-p);
            // Include a factor to prevent (1-p)^(-N) from overflow
            // Let a0 = exp(-64)*(1-p)^(-N) & T = exp(-64)
            a = expf(a - 64.0f);
            T = 1.60381083e-28f; /* exp(-64) */
        }
    }

    d = 0.0f;
    if (U > 0.5f) d = 1.0e-5f;

    S = T - a*(U-d);

    // bottom-up summation
    // NOTE: For N > 2^24 there is an extra slight loss of accuracy due to
    // incremental updates m -= 1.0f being lost as N - 1 == N. Effect on
    // ultimate outcome, however, is minimal.
    while (S < 0.0f && m > 1.0f) {
        n += 1.0f;
        m -= 1.0f;
        T  = m*p*T;
        S  = n*q*S + T;
        a  = n*q*a;
    }

    // top-down summation if needed
    if (S < 2.0f*d*a) {
        // Important to first multiply V*a to keep accuracy
        while (m > 1.0f && T*1.0e8f > (V*a)) {
            n += 1.0f;
            m -= 1.0f;
            T  = m*p*T;
            a  = n*q*a;
            // Scale down a_n & T_n if a_n is close to overflow
            if (a > 1.0e+25f) {
                // NOTE: T > 1e-8*2^-149*a > 1e-28 so no need
                // need to check if T is large enough to rescale
                T *= 1.0e-10f;
                a *= 1.0e-10f;
            }
        }
        // Scale down a_n & T_n if a_n or T_n can reach overflow in final
        // iterations.
        if (a > 1.0e+14f && T > 1.0f ) {
            T *= 1.0e-14f;
            a *= 1.0e-14f;
        }

        S = T - a*V;

        while (S < 0.0f) {
            T  = n*q*T;
            S  = m*p*S + T;
            n -= 1.0f;
            m += 1.0f;
        }
    }
    return n;
}

#undef mu_0
#undef N_SWITCH

// binoinvf_v(u, N, p):
//
// This single precision function computes the inverse
// of the binomial CDF using an approach more suitable for SIMD execution
//
// u   = CDF value in range [0,1]
// N   = Number of trials
// p   = Success probability for each trial
//
//  max |error| no more than 1 for N*p < ~1e+07
//  ave |error| < 5.14e-07*sqrt(max(0.5,N*p*(1-p))) for 0 <= N*p*(1-p) <  8
//              < 2.57e-07*sqrt(max(32 ,N*p*(1-p))) for 8 <= N*p*(1-p) < ~1e+07
//
//  For N*p >= ~1e+07, the errors will be about 1 ulp.
//

float binoinvf_v(float U, float N, float p) {

    float V = 1.0f - U;
    float q = 1.0f - p;

    float X=0.0f, NP, NQ, W, S, Del;

// handle exceptions -- constants defined in <math.h>

    // handles NAN inputs as well
    if (!(N == truncf(N)) || N < 0.0f || isinf(N)) return FP_NAN;
    if (!(U >= 0.0f && V >= 0.0f)) return FP_NAN;
    if (!(p >= 0.0f && q >= 0.0f)) return FP_NAN;
    if (U == 0.0f || p == 0.0f) return 0.0f;
    if (V == 0.0f || q == 0.0f) return N;

    // To maintain accuracy for large N work with min(p,q)
    int flag = (p > 0.5f);
    W = U; // copy U to W to be used in argument for norminvf
    // Swap p & q and U & V.
    if (flag) {
        W=p; p=q; q=W;
        W=U; U=V; V=W;
    }

    NP = N*p; NQ = N*q;

#define mu_0 10.0f

// If Np or Nq is small ( <= mu_0), i.e. p is close to 0 or 1 we perform a
// direct summation. The number of loop iterations is then bounded
// by ~ 30 + Np + 17*sqrt(Np) for Np small,
// by ~180 + Nq + 46*sqrt(Nq) for Nq small.
//
// For large Np and Nq (i.e. > mu_0) we attempt to use asymptotic
// approximations rather than direct summation.

    if (NP > mu_0 && NQ > mu_0) {

        W = norminvf(W);
        if (flag) W = -W;

// use Temme uniform approximation

        float nu, eta0, eta1, eta, sqpq, pr, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9,
               xi, xi_min_p, xi_eta;
        sqpq = sqrtf(p*q);
        nu   = N + 1.0f;

        eta0 = -W / sqrtf(nu);

        // References refer to 2020 paper by Gil, Segura & Temme
        // "ASYMPTOTIC INVERSION OF THE BINOMIAL AND NEGATIVE BINOMIAL
        // CUMULATIVE DISTRIBUTION FUNCTIONS"
        // DOI: 10.1553/etna_vol52s270

        // use asymptotic expansions to solve for xi when eta near 0
        if (fabsf(eta0) < 5.8e-1f*sqpq) {
            eta0 = eta0 / sqpq;
            pr = p - 0.5f;

            // Taylor approximation to xi - p

            // first coeffs from page 273 below (2.7)
            a2 = pr*(1.0f/3.0f);
            a3 =    (pr*pr - 0.75f)*(1.0f/36.0f);
            a4 = pr*(2.25f - pr*pr)*(1.0f/270.0f);
            a5 =    (-2.4375f + pr*pr*(-13.5f + pr*pr))*(1.0f/4320.0f);
            a6 = pr*(14.0625f + pr*pr*(19.5f + pr*pr))*(1.0f/17010.0f);
            a7 =   -(  13095.0f + pr*pr*(  270756.0f + pr*pr*(143568.0f
                            + 8896.0f*pr*pr)))*(1.0f/348364800.0f);
            a8 = pr*(   1287.0f + pr*pr*(    7692.0f + pr*pr*(1872.0f
                            + 64.0f*pr*pr)))*(1.0f/13063680.0f);
            a9 =   -(2049867.0f + pr*pr*(90331632.0f + pr*pr*(234206496.0f
                            + pr*pr*(28628736.0f + pr*pr*146176.0f))))*(1.0f/601974374400.0f);
            // expansion on page 273 in (2.7)
            xi_min_p = a9;
            xi_min_p = a8   + eta0*xi_min_p;
            xi_min_p = a7   + eta0*xi_min_p;
            xi_min_p = a6   + eta0*xi_min_p;
            xi_min_p = a5   + eta0*xi_min_p;
            xi_min_p = a4   + eta0*xi_min_p;
            xi_min_p = a3   + eta0*xi_min_p;
            xi_min_p = a2   + eta0*xi_min_p;
            xi_min_p = 1.0f + eta0*xi_min_p;
            xi_min_p *= -(p - p*p)*eta0;

            xi = p + xi_min_p;

            // Taylor approximation to eta1 * xi_eta
            // NOTE 1: use a direct expansion of the product, rather than
            // using separate expansions for eta1 and xi_eta.
            // NOTE 2: only compute xi_eta*eta1 to fifth order because
            // it is not multiplied by nu in S = nu*xi + eta1*xi_eta
            a0 = pr*(-2.0f/3.0f);
            a1 =    (  3.75f  - 11.0f*pr*pr)*(1.0f/36.0f);
            a2 = pr*(-15.75f  +  7.0f*pr*pr)*(1.0f/810.0f);
            a3 =    ( 30.375f + pr*pr*(-2.25f + pr*pr*71.0f))*(1.0f/9720.0f);
            a4 = pr*(-58.375f + pr*pr*(43.0f  + pr*pr*38.0f))*(1.0f/15120.0f);
            a5 =    (104574.375f + pr*pr*(1104205.5f + pr*pr*(-822006.0f
                            + pr*pr*147688.0f)))*(1.0f/391910400.0f);

            eta1 = a5;
            eta1 = a4 + eta0*eta1;
            eta1 = a3 + eta0*eta1;
            eta1 = a2 + eta0*eta1;
            eta1 = a1 + eta0*eta1;
            eta1 = a0 + eta0*eta1;

            xi_eta = 1.0f; // xi_eta is already incorporated in eta1
        }
        // use Newton iteration to solve for xi in (2.3) on page 272
        // when eta not near 0
        else {
            xi = p;
            eta = 0.0f;
            xi_eta = -sqpq;

            for (int k = 0; k < 5; k++) {
              xi -= (eta - eta0)*xi_eta;
              xi = fmaxf(fminf(xi, 1.0f - 1.0f/nu), 1.0f/nu);   // clamp at two ends
              eta = sqrtf(-2.0f*(xi*log1pf((p - xi)/xi) + (1.0f-xi)*log1pf(-(p - xi)/(1.0f - xi))));
              if (p < xi) eta = -eta;
              xi_eta = - eta / log1pf((p - xi)/(q*xi));   // log1p used for accuracy
            }

            eta1 = logf(sqrtf(xi - xi*xi)*eta0/(p - xi)) / eta0;   // see (3.6)

            xi_min_p = xi - p;
        }

        Del = 2.25e-2f*(1.0f/(nu*(xi - xi*xi)));

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
        if ((X - NP) >= S) {         // If X += NP was rounded up
            X = nextafterf(X,0.0f);  // use round down mode instead
        }
        X = truncf(X);

        if (10.0f <= X && X <= N - 10.0f) {
            // If NO_CDF is defined we skip this step, this increases speed
            // but can significantly impact accuracy
#ifdef NO_CDF
            return (flag) ? N - X : X;
#endif

            // Accurate enough, no need to check the CDF
            if (S - 2.0f*Del >= X - NP || X > 9.007199254740992e+15f ) {
                return (flag) ? N - X : X;
            }
            // Not accurate enough, need to check the CDF
            //
            // Calculate both C0 = F(X-1) and C1 = 1 - F(X-1) and
            // check both U against C0 and V against C1 (if necessary) to
            // prevent loss of accuracy in the tails
            float C0, C1;
            X -= 1.0f;
            binom_cdff(X, N, p, q, &C0, &C1);
            if ( ( U <= 0.5f && C0 < U ) || (V < 0.5f && V < C1 ) ) {
                X += 1.0f;
            }
            return (flag) ? N - X : X;
        }
        // At this point know that either X < 10 or X > N - 10.
        // Decide if p & q and U & V are swapped depending on whether direct
        // summation will be used when in the right tail of the distribution.
        if (U > 0.5f) {
            flag = !flag;
            W=p; p=q; q=W;
            W=U; U=V; V=W;
        }
    }

// Use direct summation routine either if Np or N(1-p) is small,
// or U is too close to 0 or 1 for Temme approximation to be valid
    X = binoinvf_summation(U, V, N, p, q);
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
// The relative error is less than 5 ULP in single precision, and   //
// the accuracy is verified in the accompanying MATLAB code         //
// as241_test.m                                                     //
//                                                                  //
//////////////////////////////////////////////////////////////////////

static float normcdfinvf_as241(float p) {

    float q, r, num, den, res;

    q = p - 0.5f;
    if (fabsf(q) <= 0.425f) {
        r = 0.180625f - q*q;

        num =         5.9109374720e+1f;
        num = r*num + 1.5929113202e+2f;
        num = r*num + 5.0434271938e+1f;
        num = r*num + 3.3871327179e+0f;

        den =         6.7187563600e+1f;
        den = r*den + 7.8757757664e+1f;
        den = r*den + 1.7895169469e+1f;
        den = r*den + 1.0000000000e+0f;

        res = q * num / den;

        return res;
    }

    else {

        if (q < 0.0f)
            r = p;
        else
            r = 1.0f - p;

        r = sqrtf(-logf(r));

        if (r <= 5.0f) {
            r = r - 1.6f;

            num =         1.7023821103e-1f;
            num = r*num + 1.3067284816e+0f;
            num = r*num + 2.7568153900e+0f;
            num = r*num + 1.4234372777e+0f;

            den =         1.2021132975e-1f;
            den = r*den + 7.3700164250e-1f;
            den = r*den + 1.0000000000e+0f;

            res = num / den;
        }

        else {
            r = r - 5.0f;

            num =         1.7337203997e-2f;
            num = r*num + 4.2868294337e-1f;
            num = r*num + 3.0812263860e+0f;
            num = r*num + 6.6579051150e+0f;

            den =         1.2258502635e-2f;
            den = r*den + 2.4197894225e-1f;
            den = r*den + 1.0000000000e+0f;

            res = num / den;
        }

        if (q < 0.0f)
            res = - res;

        return res;
    }
}

#endif /* BINOINVF_H */
