// Copyright (c) Radost Waszkiewicz - 2021
// This code is licensed under MIT license

#ifndef POINT_UTILS_H
#define POINT_UTILS_H

#include"points.h" // Point, Curve
#include<random> // Normal distribution

namespace pt
{

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,1.0);

    double GetNormal()
    {
      return distribution(generator);
    }

    Point GetSpherical()
    {
      auto p = Point();
      while(p.Norm() < 0.0001)
      {
        p.x = GetNormal();
        p.y = GetNormal();
        p.z = GetNormal();
      }
      return p.Normalized();
    }

    Point GetSphericalNormal()
    {
        auto p = Point();
        p.x = GetNormal();
        p.y = GetNormal();
        p.z = GetNormal();
        return p;
    }

    Curve Symmetrized(Curve c)
    {
        Curve phase1(c);
        int n = phase1.size();
        for(int i=0; i<n; i++)
        {
            pt::Point l = c[i];
            pt::Point r = c[(n-i)%n];
            phase1[i] = pt::Point( 0.5*(l.x-r.x), 0.5*(l.y-r.y), 0.5*(l.z+r.z) );
        }
        if(n%2!=0)
        {
            return phase1;
        }
        else
        {
            Curve phase2(phase1);
            for(int i=0; i<n; i++)
            {
                pt::Point l = phase1[i];
                pt::Point r = phase1[(n-i+(n/2))%n];
                phase2[i] = pt::Point( 0.5*(l.x+r.x), 0.5*(l.y-r.y), 0.5*(l.z-r.z) );
            }
            return phase2;
        }
    }

    double DebeyeEnergyDiscretized(Curve discretized, double DebeyeLength)
    {
        double accum = 0;
        int n = discretized.size();
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<i;j++)
            {
                double d = (discretized[i] - discretized[j]).Norm();
                if(d != 0)
                {
                    // (1.0 / n) adjustment to decrease constant term
                    accum += (1.0 / (n*n-1))*((1.0  / d)*exp(-d/DebeyeLength) - (1.0 / n));
                }
                else
                {
                    accum += 1000.0;
                }
            }
        }
        return accum;
    }

    // Based on algorithm:
    // https://en.wikipedia.org/wiki/Writhe#Numerically_approximating_the_Gauss_integral_for_writhe_of_a_curve_in_space
    double WritheDiscretized(Curve discretized)
    {
        int n = discretized.size();
        double accum = 0;

        for(int i=0; i<n; i++)for(int j=0; j<n; j++)
        {
            if(i == j) accum += 0;
            else if((n+i-j)%n == 1 || (n+j-i)%n == 1) accum +=0;
            else if(i>j)
            {
                Point r1 = discretized[i%n];
                Point r2 = discretized[(i+1)%n];
                Point r3 = discretized[j%n];
                Point r4 = discretized[(j+1)%n];

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

    double NearNeighboursDiscretized(Curve discretized, double Near)
    {
        //TODO: more efficient implementation with octal hashes

        int len = discretized.size();

        std::vector<double> arccoords(len);
        for(int i=1; i < len; i++)
        {
            arccoords[i] = arccoords[i-1] + (discretized[i] - discretized[i-1]).Norm();
        }

        double accum = 0;
        for(int i=0; i < len; i++)
        {
            Point l = discretized[i];
            for(int j=0; j < i; j++)
            {
                Point r = discretized[j];

                double distance = (l-r).Norm();

                if( distance < Near
                  && fmod(arccoords[i]-arccoords[j]+arccoords[len-1], arccoords[len-1]) > 0.501*M_PI * Near 
                  && fmod(arccoords[j]-arccoords[i]+arccoords[len-1], arccoords[len-1]) > 0.501*M_PI * Near)
                {
                    accum += 0.5*fmod(
                                       2*arccoords[len-1] + arccoords[(len+i+1)%len] - arccoords[(len+i-1)%len],
                                       arccoords[len-1]
                                     ); // length of segment near i
                    break;
                }
            }
        }
        return accum;
    }
}

#endif /* POINT_UTILS_H */
