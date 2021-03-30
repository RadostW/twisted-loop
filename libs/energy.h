#ifndef ENERGY_H
#define ENERGY_H

#include "randompts.h"
#include "simpsons.h"

namespace wr
{

    sp::spline gammaX, gammaY, gammaZ;
    wr::Curve wrSplineCurve;
    const int wrpoints = 200; // Number of subsampling points for writhe calc.

    wr::Curve nearSplineCurve;
    int nearSamplingN = 1000; // Number of subsampling points for steric interactions

    int splinePoints; // Stores number of points of initCurve

    double dl(double s);

    void curveInit(const wr::Curve& initCurve)
    {
        int n = initCurve.size();
        splinePoints = n;
        std::vector<double> xi(n),yi(n),zi(n);
        for(int i=0;i<n;i++)
        {
            xi[i] = initCurve.at(i).x;
            yi[i] = initCurve.at(i).y;
            zi[i] = initCurve.at(i).z;
        }
        gammaX.setPoints(xi);
        gammaY.setPoints(yi);
        gammaZ.setPoints(zi);

        wrSplineCurve = wr::Curve(wrpoints);
        double wrdx = 1.0/wrpoints;
        for(int i=0; i < wrpoints; i++)
        {
            wrSplineCurve.at(i) = wr::Point(gammaX(wrdx*i),gammaY(wrdx*i),gammaZ(wrdx*i));
        }

        nearSplineCurve = wr::Curve(nearSamplingN);
        double nedx = 1.0/nearSamplingN;
        for(int i=0; i < nearSamplingN; i++)
        {
            nearSplineCurve.at(i) = wr::Point(gammaX(nedx*i),gammaY(nedx*i),gammaZ(nedx*i));
        }

        return;
    }

    // Remember to run curveInit for each curve.
    double curvature(double s)
    {
        wr::Point gammaPrime(gammaX.prime(s),gammaY.prime(s),gammaZ.prime(s));
        wr::Point gammaBis(gammaX.bis(s),gammaY.bis(s),gammaZ.bis(s));
        return (gammaPrime^gammaBis).Norm()/pow(gammaPrime.Norm(),3);
    }

    // Remember to run curveInit for each curve.
    double writhe()
    {
        wr::writhe(wrSplineCurve);
    }

    double nearNeighbours(double near)
    {
        std::vector<double> arccoords(nearSamplingN);
        for(int i=1; i < nearSamplingN; i++)
        {
            arccoords[i] = arccoords[i-1] + (nearSplineCurve[i] - nearSplineCurve[i-1]).Norm();
        }

        double accum = 0;
        for(int i=0; i<nearSamplingN; i++)
        {
            wr::Point r2 = nearSplineCurve[(nearSamplingN+i)%nearSamplingN];
            for(int j=0; j<nearSamplingN; j++)
            {
                wr::Point r5 = nearSplineCurve[(nearSamplingN+j)%nearSamplingN];

                if((r5-r2).Norm() < near && 
                        fmod(arccoords[i] - arccoords[j]+arccoords[nearSamplingN-1],arccoords[nearSamplingN-1]) > 0.501*M_PI * near &&
                        fmod(arccoords[j] - arccoords[i]+arccoords[nearSamplingN-1],arccoords[nearSamplingN-1]) > 0.501*M_PI * near
                )
                {
                    accum += 0.5*fmod(2*arccoords[nearSamplingN-1] + arccoords[(nearSamplingN+i+1)%nearSamplingN] - arccoords[(nearSamplingN+i-1)%nearSamplingN],arccoords[nearSamplingN-1]);
                    //printf("near");
                    break;
                }
            }
        }
        return accum;

    }

    double squaredCurvature(double s)
    {
        double tmp = curvature(s);
        return tmp*tmp;
    }

    double dl(double s)
    {
        wr::Point gammaPrime(gammaX.prime(s),gammaY.prime(s),gammaZ.prime(s));
        return gammaPrime.Norm();
    }

    double squaredCurvatureIntegrand(double s)
    {
        return squaredCurvature(s)*dl(s);
    }

    // Remember to run curveInit for each curve
    double totalSquaredCurvature()
    {
        double accum = 0;
        int n = splinePoints;
        double dx = 1.0/n;
        for(int i=0;i<n;i++)
        {
            accum += si::adaptiveSimpsons(squaredCurvatureIntegrand,i*dx,(i+1)*dx, 0.0001, 10);
        }
        return accum;
    }

    // Remember to run curveInit for each curve
    double totalLength()
    {
        return si::adaptiveSimpsons(dl,0,1, 0.0001, 10);
    }

    double inverseDistances()
    {
        double accum = 0;
        int n = wrpoints;
        for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)
        {
            accum += 1.0 / (0.1+(wrSplineCurve[i]-wrSplineCurve[j]).Norm());
        }
        return accum;
    }






















    double energy(Curve c,double Link,double Radius,double Length,double Omega,double Electrostatic,double SoftContact)
    {
        curveInit(c);
        double curveLength = totalLength();
        double squaredCurvature = totalSquaredCurvature();
        double Tw = Link - writhe();
        double TwistDensity = 2*M_PI*Tw / curveLength;
        double contact = nearNeighbours(2*Radius);
        double softContact = nearNeighbours(5*Radius);
        double inverseDist = inverseDistances();

        double lengthPenalty = 8000;
        double contactPenalty = 10000;

        return lengthPenalty*(curveLength-Length)*(curveLength-Length) + squaredCurvature + 
                curveLength*TwistDensity*TwistDensity*Omega + contactPenalty*contact +
                Electrostatic * inverseDist + softContact*SoftContact;
    }

    Curve symmetrized(Curve c);
    Curve nudge(Curve c,double ammount)
    {
        Curve tmp(c);
        int n = c.size();
        for(int i=0;i<n;i++)
        {
                tmp[i] = tmp[i]+GetSphericalNormal()*ammount;
        }
        return symmetrized(tmp);
    }


    Curve symmetrized(Curve c)
    {
        Curve phase1(c);
        int n = phase1.size();
        for(int i=0; i<n; i++)
        {
            Point l = c[i];
            Point r = c[(n-i)%n];
            phase1[i] = wr::Point( 0.5*(l.x-r.x), 0.5*(l.y-r.y), 0.5*(l.z+r.z) );
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
                Point l = phase1[i];
                Point r = phase1[(n-i+(n/2))%n];
                phase2[i] = wr::Point( 0.5*(l.x+r.x), 0.5*(l.y-r.y), 0.5*(l.z-r.z) );
            }
            return phase2;
        }
    }

}

#endif /* ENERGY_H */

