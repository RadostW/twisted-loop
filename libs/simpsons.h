// Sourced from wikipedia -- Adaptive Simpsons Rule

#ifndef SIMPSONS_H
#define SIMPSONS_H

#include <math.h>  // include file for fabs and sin
#include <stdio.h> // include file for printf and perror
#include <errno.h>

namespace si
{

/** Adaptive Simpson's Rule, Recursive Core */
double adaptiveSimpsonsAux(double (*f)(double), double a, double b, double eps,
                          double whole, double fa, double fb, double fm, int rec) {
    double m   = (a + b)/2,  h   = (b - a)/2;
    double lm  = (a + m)/2,  rm  = (m + b)/2;
    // serious numerical trouble: it won't converge
    if ((eps/2 == eps) || (a == lm)) { errno = EDOM; return whole; }
    double flm = f(lm),      frm = f(rm);
    double left  = (h/6) * (fa + 4*flm + fm);
    double right = (h/6) * (fm + 4*frm + fb);
    double delta = left + right - whole;

    if (rec <= 0 && errno != EDOM) errno = ERANGE;  // depth limit too shallow
    // Lyness 1969 + Richardson extrapolation; see article
    if (rec <= 0 || fabs(delta) <= 15*eps)
        return left + right + (delta)/15;
    return adaptiveSimpsonsAux(f, a, m, eps/2, left,  fa, fm, flm, rec-1) +
           adaptiveSimpsonsAux(f, m, b, eps/2, right, fm, fb, frm, rec-1);
}

/** Adaptive Simpson's Rule Wrapper
 *  (fills in cached function evaluations) */
double adaptiveSimpsons(double (*f)(double),     // function ptr to integrate
                       double a, double b,      // interval [a,b]
                       double epsilon,         // error tolerance
                       int maxRecDepth) {     // recursion cap
    errno = 0;
    double h = b - a;
    if (h == 0) return 0;
    double fa = f(a), fb = f(b), fm = f((a + b)/2);
    double S = (h/6)*(fa + 4*fm + fb);
    return adaptiveSimpsonsAux(f, a, b, epsilon, S, fa, fb, fm, maxRecDepth);
}


}
#endif /* SIMPSONS_H */
