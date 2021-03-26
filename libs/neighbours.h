#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include "point.h"
#include<cmath>
#include<cstdio>

namespace wr
{

double nearNeighbours(Curve curve,double near)
{

    int n = curve.size();
    double step = 1.0/n;
    std::vector<double> arccoords(n,0.0);

    for(int i=1;i<n;i++)
    {
        Point r1 = curve[(n+i-1)%n];
        Point r2 = curve[(n+i)%n];

        arccoords[i] = arccoords[i-1] + (r2-r1).Norm();
    }

    double accum = 0;

    for(int i=0; i<n; i++)
    {

        Point r2 = curve[(n+i)%n];

        for(int j=0; j<n; j++)
        {
            Point r5 = curve[(n+j)%n];

            if((r5-r2).Norm() < near && 
                    fmod(arccoords[i] - arccoords[j]+arccoords[n-1],arccoords[n-1]) > 0.501*M_PI * near &&
                    fmod(arccoords[j] - arccoords[i]+arccoords[n-1],arccoords[n-1]) > 0.501*M_PI * near
            )
            {
                accum += 0.5*fmod(2*arccoords[n-1] + arccoords[(n+i+1)%n] - arccoords[(n+i-1)%n],arccoords[n-1]);
                printf("near");
                break;
            }
        }
    }
    return accum;
}

}

#endif /* NEIGHBOURS_H */
