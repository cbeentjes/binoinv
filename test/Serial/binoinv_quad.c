#include <stdlib.h>

#include <math.h>
#include <quadmath.h>

#include <float.h>

//////////////////////////////////////////////////
// compute binomial p.m.f in quad precision
//////////////////////////////////////////////////

__float128 binom_pmf_quad(int n, float N, float p) {
    __float128 q_N, q_p, q_n;

    q_p = (__float128) p;
    q_N = (__float128) N;
    q_n = (__float128) n;

    if (p==0.0f) return( (n==0) ? 1.0q : 0.0q);
    if (p==1.0f) return( (n==N) ? 1.0q : 0.0q);
    if (n==0)   return expq(q_N*log1pq(-q_p));
    if (n==N)   return expq(q_N*logq(q_p));

    return expq(  q_n*logq(q_p) + (q_N-q_n)*log1pq(-q_p)
            + lgammaq(q_N+1.0q)
            - lgammaq(q_n+1.0q) - lgammaq(q_N-q_n+1.0q) );
}


//////////////////////////////////////////////////
// compute reference solution in quad precision
//////////////////////////////////////////////////

void binoinv_quad(float N, float p, float *u_lo, float *u_hi,
        double *U_lo, double *U_hi) {
    __float128 q_s;

    q_s = 0.0q;

    double NP = N*p;
    double NPQ = N*p*(1.0f-p);
    double NPQ_root = sqrt(NPQ);

    for (int n=0; n<=N; n++) {
#ifdef FAST_CHECK // Speed up check by using knowledge of distribution shape

        // For large Np(1-p) use normal approximation to determine good
        // starting point of the summation, as we do not need to evaluate
        // the PMF when it is smaller than quad precision.
        if (NPQ > 25000. && n < NP - 152.0*NPQ_root) {
            u_lo[n] = 0.0f;
            u_hi[n] = FLT_TRUE_MIN;
            U_lo[n] = 0.0;
            U_hi[n] = DBL_TRUE_MIN;
            continue;
        }

        // Similarly, no need to continue summing if we are already beyond
        // double precision.
        if (q_s > 1.0 - DBL_EPSILON/2.0) {
            u_hi[n] = 1.0f;
            u_lo[n] = 1.0f-FLT_EPSILON/2.0f;
            U_hi[n] = 1.0;
            U_lo[n] = 1.0-DBL_EPSILON/2.0;
            continue;
        }
        // No need to continue too far into the right tail of the distribution
        if (n > 10.0 + NP + 20.0*NPQ_root) {
            break;
        }
#endif /* FAST_CHECK */

        q_s += binom_pmf_quad(n, N, p);

        //
        // single precision interval bisection
        //
#ifdef FAST_CHECK
        if (q_s < FLT_TRUE_MIN ) {
            u_hi[n] = FLT_TRUE_MIN;
            u_lo[n] = 0.0f;
        } else if (q_s > 1.0f - FLT_EPSILON/2.0f) {
            u_hi[n] = 1.0f;
            u_lo[n] = 1.0f-FLT_EPSILON/2.0f;
        } else {
#endif
            u_hi[n] = 1.0f;
            u_lo[n] = 0.0f;
            float u_mid = 0.5f*(u_hi[n] + u_lo[n]);

            while ( (u_mid>u_lo[n]) & (u_mid<u_hi[n]) ) {
                if ((__float128) u_mid > q_s)
                    u_hi[n] = u_mid;
                else
                    u_lo[n] = u_mid;

                u_mid = 0.5f*(u_hi[n] + u_lo[n]);
            }
#ifdef FAST_CHECK
        }
#endif

        //
        // double precision interval bisection
        //
#ifdef FAST_CHECK
        if (q_s < DBL_TRUE_MIN ) {
            U_hi[n] = DBL_TRUE_MIN;
            U_lo[n] = 0.0;
        } else if (q_s > 1.0 - DBL_EPSILON/2.0) {
            U_hi[n] = 1.0;
            U_lo[n] = 1.0-DBL_EPSILON/2.0;
        } else {
#endif
            U_hi[n] = 1.0;
            U_lo[n] = 0.0;
            double U_mid = 0.5*(U_hi[n] + U_lo[n]);

            while ( (U_mid>U_lo[n]) & (U_mid<U_hi[n]) ) {
                if ((__float128) U_mid > q_s)
                    U_hi[n] = U_mid;
                else
                    U_lo[n] = U_mid;

                U_mid = 0.5*(U_hi[n] + U_lo[n]);
            }
#ifdef FAST_CHECK
        }
#endif
    }
}


void binocinv_quad(float N, float p, float *u_lo, float *u_hi,
        double *U_lo, double *U_hi) {
    __float128 q_s;

    q_s = 0.0q;

    double NP = N*p;
    double NPQ = N*p*(1.0-p);
    double NPQ_root = sqrt(NPQ);

    for (int n=N; n>=0; n--) {
#ifdef FAST_CHECK // Speed up check by using knowledge of distribution shape

        // For large Np(1-p) use normal approximation to determine good
        // starting point of the summation, as we do not need to evaluate
        // the PMF when it is smaller than quad precision.
        if (NPQ > 25000. && n > NP + 152.0*NPQ_root) {
            u_lo[n] = 0.0f;
            u_hi[n] = FLT_TRUE_MIN;
            U_lo[n] = 0.0;
            U_hi[n] = DBL_TRUE_MIN;
            continue;
        }

        // Similarly, no need to continue summing if we are already beyond
        // double precision.
        if (q_s > 1.0 - DBL_EPSILON/2.0) {
            u_hi[n] = 1.0f;
            u_lo[n] = 1.0f-FLT_EPSILON/2.0f;
            U_hi[n] = 1.0;
            U_lo[n] = 1.0-DBL_EPSILON/2.0;
            continue;
        }
#endif /* FAST_CHECK */

        q_s += binom_pmf_quad(n, N, p);

        //
        // single precision interval bisection
        //
#ifdef FAST_CHECK
        if (q_s < FLT_TRUE_MIN ) {
            u_hi[n] = FLT_TRUE_MIN;
            u_lo[n] = 0.0f;
        } else if (q_s > 1.0f - FLT_EPSILON/2.0f) {
            u_hi[n] = 1.0f;
            u_lo[n] = 1.0f-FLT_EPSILON/2.0f;
        } else {
#endif
            u_hi[n] = 1.0f;
            u_lo[n] = 0.0f;
            float u_mid = 0.5f*(u_hi[n] + u_lo[n]);

            while ( (u_mid>u_lo[n]) & (u_mid<u_hi[n]) ) {
                if ((__float128) u_mid > q_s)
                    u_hi[n] = u_mid;
                else
                    u_lo[n] = u_mid;

                u_mid = 0.5f*(u_hi[n] + u_lo[n]);
            }
#ifdef FAST_CHECK
        }
#endif

        //
        // double precision interval bisection
        //
#ifdef FAST_CHECK
        if (q_s < DBL_TRUE_MIN ) {
            U_hi[n] = DBL_TRUE_MIN;
            U_lo[n] = 0.0;
        } else if (q_s > 1.0 - DBL_EPSILON/2.0) {
            U_hi[n] = 1.0;
            U_lo[n] = 1.0-DBL_EPSILON/2.0;
        } else {
#endif
            U_hi[n] = 1.0;
            U_lo[n] = 0.0;
            double U_mid = 0.5*(U_hi[n] + U_lo[n]);

            while ( (U_mid>U_lo[n]) & (U_mid<U_hi[n]) ) {
                if ((__float128) U_mid > q_s)
                    U_hi[n] = U_mid;
                else
                    U_lo[n] = U_mid;

                U_mid = 0.5*(U_hi[n] + U_lo[n]);
            }
#ifdef FAST_CHECK
        }
#endif
    }
}
