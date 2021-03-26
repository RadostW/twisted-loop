// Copyright (c) Radost Waszkiewicz - 2021
// This code is licensed under MIT license
// Computes writhe of parametric curve given as
// std::vector<Point>, with Point representing 3D vectors.

#ifndef WRITHE_H
#define WRITHE_H

#include <cassert>
#include <vector>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include "point.h"

namespace wr
{

typedef std::vector<Point> Curve;

// Based on algorithm:
// https://en.wikipedia.org/wiki/Writhe#Numerically_approximating_the_Gauss_integral_for_writhe_of_a_curve_in_space
double writhe(Curve curve)
{
    int n = curve.size();
    double accum = 0;

    for(int i=0; i<n; i++)for(int j=0; j<n; j++)
    {
        if(i == j) accum += 0;
        else if((n+i-j)%n == 1 || (n+j-i)%n == 1) accum +=0;
        else if(i>j)
        {
            Point r1 = curve[i%n];
            Point r2 = curve[(i+1)%n];
            Point r3 = curve[j%n];
            Point r4 = curve[(j+1)%n];

            Point r13 = r3 - r1;
            Point r14 = r4 - r1;
            Point r24 = r4 - r2;
            Point r23 = r3 - r2;

            Point n1 = (r13^r14).Normalized();
            Point n2 = (r14^r24).Normalized();
            Point n3 = (r24^r23).Normalized();
            Point n4 = (r23^r13).Normalized();


            // min/max to account for very small/large values of cross products
            double OmegaStar = asin( std::max(std::min( n1 * n2 ,1.),-1.) )+
                               asin( std::max(std::min( n2 * n3 ,1.),-1.) )+
                               asin( std::max(std::min( n3 * n4 ,1.),-1.) )+
                               asin( std::max(std::min( n4 * n1 ,1.),-1.) );
            
            Point r34 = r4 - r3;
            Point r12 = r2 - r1;

            double Omega = OmegaStar * (((r34^r12)*r13)>0.0? 1.0 : -1.0);
            if(!isfinite(Omega))
            {
                printf("Debug info:\n");
                printf("norms:\n");
                printf("r13^r14: %lf\n",(r13^r14).Norm());
                printf("r13^r24: %lf\n",(r13^r24).Norm());
                printf("r24^r23: %lf\n",(r24^r23).Norm());
                printf("r23^r13: %lf\n",(r23^r13).Norm());
                printf("arcsines:\n");
                printf("asin(n1*n2): %lf\n",asin(n1*n2));
                printf("asin(n2*n3): %lf\n",asin(n2*n3));
                printf("asin(n3*n4): %lf\n",asin(n3*n4));
                printf("asin(n4*n1): %lf\n",asin(n4*n1));
                assert(isfinite(Omega));
            }

            accum += Omega;
        }
    }    

    return accum / (2*M_PI);    
   
}

}
#endif /* WRITHE_H */
