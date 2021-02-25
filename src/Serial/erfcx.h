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

#ifndef ERFCX_H
#define ERFCX_H

#include <math.h>

//////////////////////////////////////////////////////////////////////
//
// Implementation of the scaled complementary error function erfcx(x)
// based on a Chebyshev approximation of the function
// (1+2x)*erfcx(x) after scaling the real axis to the interval [-1,1]
// using q = (x - 4)/(x + 4)
//
// Published details in:
// Mathematics of Computation, VOL. 36, NO. 153, January 1981
// DOI: 10.1090/S0025-5718-1981-0595058-X
//
// Max error in the half-plane is claimed to be ~3.0 ULPs.
// More details on implementation, error and source code on:
// https://stackoverflow.com/questions/39777360/accurate-computation-of-scaled-complementary-error-function-erfcx
//
// If FMA operation is not supported and separate multiplication + addition
// is used the max error seems to increase by ~0.5 ULPs worst case.
//
//////////////////////////////////////////////////////////////////////

double erfcx (double x) {

    double a, d, e, m, p, q, r, t;

#ifndef __FAST_CDF
    a = fmax(x, 0.0 - x); // NaN preserving absolute value computation
#else
    a = x;
#endif /* FAST_CDF */

    /* Compute q = (a-4)/(a+4) accurately. [0,INF) -> [-1,1] */
    m = a - 4.0;
    p = a + 4.0;
    r = 1.0 / p;
    q = m * r;
    t = fma(q + 1.0, -4.0, a);
    e = fma(q, -a, t);
    q = fma(r, e, q);

    /* Approximate (1+2*a)*exp(a*a)*erfc(a) as p(q)+1 for q in [-1,1] */
    p =            0x1.edcad78fc8044p-31;  //  8.9820305531190140e-10
    p = fma(p, q,  0x1.b1548f14735d1p-30); //  1.5764464777959401e-09
    p = fma(p, q, -0x1.a1ad2e6c4a7a8p-27); // -1.2155985739342269e-08
    p = fma(p, q, -0x1.1985b48f08574p-26); // -1.6386753783877791e-08
    p = fma(p, q,  0x1.c6a8093ac4f83p-24); //  1.0585794011876720e-07
    p = fma(p, q,  0x1.31c2b2b44b731p-24); //  7.1190423171700940e-08
    p = fma(p, q, -0x1.b87373facb29fp-21); // -8.2040389712752056e-07
    p = fma(p, q,  0x1.3fef1358803b7p-22); //  2.9796165315625938e-07
    p = fma(p, q,  0x1.7eec072bb0be3p-18); //  5.7059822144459833e-06
    p = fma(p, q, -0x1.78a680a741c4ap-17); // -1.1225056665965572e-05
    p = fma(p, q, -0x1.9951f39295cf4p-16); // -2.4397380523258482e-05
    p = fma(p, q,  0x1.3be1255ce180bp-13); //  1.5062307184282616e-04
    p = fma(p, q, -0x1.a1df71176b791p-13); // -1.9925728768782324e-04
    p = fma(p, q, -0x1.8d4aaa0099bc8p-11); // -7.5777369791018515e-04
    p = fma(p, q,  0x1.49c673066c831p-8);  //  5.0319701025945277e-03
    p = fma(p, q, -0x1.0962386ea02b7p-6);  // -1.6197733983519948e-02
    p = fma(p, q,  0x1.3079edf465cc3p-5);  //  3.7167515521269866e-02
    p = fma(p, q, -0x1.0fb06dfedc4ccp-4);  // -6.6330365820039094e-02
    p = fma(p, q,  0x1.7fee004e266dfp-4);  //  9.3732834999538536e-02
    p = fma(p, q, -0x1.9ddb23c3e14d2p-4);  // -1.0103906603588378e-01
    p = fma(p, q,  0x1.16ecefcfa4865p-4);  //  6.8097054254651804e-02
    p = fma(p, q,  0x1.f7f5df66fc349p-7);  //  1.5379652102610957e-02
    p = fma(p, q, -0x1.1df1ad154a27fp-3);  // -1.3962111684056208e-01
    p = fma(p, q,  0x1.dd2c8b74febf6p-3);  //  2.3299511862555250e-01

    /* Divide (1+p) by (1+2*a) ==> exp(a*a)*erfc(a) */
    d = a + 0.5;
    r = 1.0 / d;
    r = r * 0.5;
    q = fma(p, r, r); // q = (p+1)/(1+2*a)
    t = q + q;
    e = (p - q) + fma(t, -a, 1.0); // residual: (p+1)-q*(1+2*a)
    r = fma(e, r, q);

    // NOTE: case checks are not needed for binoinv project as erfcx code is
    // only ever called for 0 <= x < 30.
#ifndef __FAST_CDF
    /* Handle argument of infinity */
    if (a > 0x1.fffffffffffffp1023) r = 0.0; // 1.797693134862316e+308 //

    /* Handle negative arguments: erfcx(x) = 2*exp(x*x) - erfcx(|x|) */
    if (x < 0.0) {
        double s = x * x;
        d = fma(x, x, -s);
        e = exp(s);
        r = e - r;
        r = fma(e, d + d, r);
        r = r + e;
        if (e > 0x1.fffffffffffffp1023) r = e; // 1.797693134862316e+308 // avoid creating NaN
    }
#endif /* FAST_CDF */
    return r;
}

#endif /* ERFCX_H */
