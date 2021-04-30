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

#ifndef ERFCXF_H
#define ERFCXF_H

#include <math.h>

//////////////////////////////////////////////////////////////////////
//
// Implementation of the scaled complementary error function erfcx(x)
// based on a Chebyshev approximation of the function
// (1+2x)*erfcx(x) after scaling the real axis to the interval [-1,1]
// using q = (x - 2)/(x + 2)
//
// Published details in:
// Mathematics of Computation, VOL. 36, NO. 153, January 1981
// DOI: 10.1090/S0025-5718-1981-0595058-X
//
// Max error in the half-plane is claimed to be ~2.0 ULPs.
// More details on implementation, error and source code by N. Juffa on:
// https://stackoverflow.com/questions/39777360/accurate-computation-of-scaled-complementary-error-function-erfcx
//
// If FMA operation is not supported and separate multiplication + addition
// is used the max error seems to increase by ~0.5 ULPs worst case.
//
//////////////////////////////////////////////////////////////////////

float erfcxf (float x) {

    float a, d, e, m, p, q, r, t;

#ifndef __FAST_CDF
    a = fmaxf(x, 0.0 - x); // NaN preserving absolute value computation
#else
    a = x;
#endif /* __FAST_CDF */

    /* Compute q = (a-2)/(a+2) accurately. [0,INF) -> [-1,1] */
    m = a - 2.0f;
    p = a + 2.0f;
    r = 1.0f / p;
    q = m * r;
    t = fmaf(q + 1.0f, -2.0f, a);
    e = fmaf(q, -a, t);
    q = fmaf(r, e, q);

    /* Approximate (1+2*a)*exp(a*a)*erfc(a) as p(q)+1 for q in [-1,1] */
    p =             0x1.f10000p-15f;  //  5.92470169e-5
    p = fmaf(p, q,  0x1.521cc6p-13f); //  1.61224554e-4
    p = fmaf(p, q, -0x1.6b4ffep-12f); // -3.46481771e-4
    p = fmaf(p, q, -0x1.6e2a7cp-10f); // -1.39681227e-3
    p = fmaf(p, q,  0x1.3c1d7ep-10f); //  1.20588380e-3
    p = fmaf(p, q,  0x1.1cc236p-07f); //  8.69014394e-3
    p = fmaf(p, q, -0x1.069940p-07f); // -8.01387429e-3
    p = fmaf(p, q, -0x1.bc1b6cp-05f); // -5.42122945e-2
    p = fmaf(p, q,  0x1.4ff8acp-03f); //  1.64048523e-1
    p = fmaf(p, q, -0x1.54081ap-03f); // -1.66031078e-1
    p = fmaf(p, q, -0x1.7bf5cep-04f); // -9.27637145e-2
    p = fmaf(p, q,  0x1.1ba03ap-02f); //  2.76978403e-1

    /* Divide (1+p) by (1+2*a) ==> exp(a*a)*erfc(a) */
    d = a + 0.5f;
    r = 1.0f / d;
    r = r * 0.5f;
    q = fmaf(p, r, r); // q = (p+1)/(1+2*a)
    t = q + q;
    e = (p - q) + fmaf(t, -a, 1.0f); // residual: (p+1)-q*(1+2*a)
    r = fmaf(e, r, q);

    // NOTE: case checks are not needed for binoinv project as erfcx code is
    // only ever called for 0 <= x < 30.
#ifndef __FAST_CDF
    /* Handle argument of infinity */
    if (a > 0x1.fffffep127f) r = 0.0f; // 3.40282347e+38 // handle INF argument

    /* Handle negative arguments: erfcx(x) = 2*exp(x*x) - erfcx(|x|) */
    if (x < 0.0f) {
        float s = x * x;
        d = fmaf(x, x, -s);
        e = expf(s);
        r = e - r;
        r = fmaf(e, d + d, r);
        r = r + e;
        if (e > 0x1.fffffep127f) r = e; // 3.40282347e+38 // avoid creating NaN
    }
#endif /* __FAST_CDF */
    return r;
}

#endif /* ERFCXF_H */
