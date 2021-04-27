// Copyright (c) Radost Waszkiewicz - 2021
// This code is licensed under MIT license
// Computes wirithe, lenght, and total curvature of a 3d curve

#ifndef SPLINE_UTILS_H
#define SPLINE_UTILS_H

#include"splines.h" // Spline
#include"points.h" // Point
#include"simpsons.h" // Integration
#include"point_utils.h" // WritheDiscretized, DebeyeEnergyDiscretized, NearNeighboursDiscretized

namespace sp
{
    const double INTEGRATION_PRECISION = 0.0001; // Smaller -> more accurate
    const int MAX_INTEGRATION_RECURSION = 10; // Larger -> more accurate
    const int WRITHE_DISCRETIZATION_POINTS = 200; // Larger -> more accurate
    const double NEAR_NEIGHBOURS_RESOLUTION = 0.1; // Smaller -> more accurate

    // Hacked to get pure function pointer to Simpsons integrator
    const Spline* _globalx;
    const Spline* _globaly;
    const Spline* _globalz;


    // Internal. Passed to integrator
    double _length_density(double s)
    {
        pt::Point gammaPrime(_globalx->prime(s),_globaly->prime(s),_globalz->prime(s));
        return gammaPrime.Norm();
    }

    double Length(const Spline &x, const Spline &y, const Spline &z)
    {
        _globalx = &x;
        _globaly = &y;
        _globalz = &z;
        return si::adaptiveSimpsons(_length_density, 0, 1, INTEGRATION_PRECISION, MAX_INTEGRATION_RECURSION);
    }

    // Internal. Passed to integrator
    double _squared_curvature_density(double s)
    {
        pt::Point gammaPrime(_globalx->prime(s),_globaly->prime(s),_globalz->prime(s));
        pt::Point gammaBis(_globalx->bis(s),_globaly->bis(s),_globalz->bis(s));
        return (gammaPrime^gammaBis).Norm()/pow(gammaPrime.Norm(),3);        
    }

    double TotalCurvature(const Spline &x, const Spline &y, const Spline &z)
    {
        _globalx = &x;
        _globaly = &y;
        _globalz = &z;
        return si::adaptiveSimpsons(_squared_curvature_density, 0, 1, INTEGRATION_PRECISION, MAX_INTEGRATION_RECURSION);        
    }

    double DebeyeEnergy(const Spline &x, const Spline &y, const Spline &z,double DebeyeLength,int BasePairs)
    {
        pt::Curve discretized(BasePairs);
        double ds = 1.0 / BasePairs;
        for(int i=0;i<BasePairs;i++)
        {
            discretized[i] = pt::Point(x(i*ds),y(i*ds),z(i*ds));
        }
        return pt::DebeyeEnergyDiscretized(discretized,DebeyeLength);
    }

    double Writhe(const Spline &x, const Spline &y, const Spline &z)
    {
        pt::Curve discretized(WRITHE_DISCRETIZATION_POINTS);
        double ds = 1.0 / WRITHE_DISCRETIZATION_POINTS;
        for(int i=0;i<WRITHE_DISCRETIZATION_POINTS;i++)
        {
            discretized[i] = pt::Point(x(i*ds),y(i*ds),z(i*ds));
        }
        return pt::WritheDiscretized(discretized);
    }

    double NearNeighbours(const Spline &x, const Spline &y, const Spline &z, double Near)
    {
        int points = NEAR_NEIGHBOURS_RESOLUTION * 1.0 / Near;
        pt::Curve discretized(points);
        double ds = 1.0 / points;
        for(int i=0;i<points;i++)
        {
            discretized[i] = pt::Point(x(i*ds),y(i*ds),z(i*ds));
        }
        return pt::NearNeighboursDiscretized(discretized,Near);
    }

}

#endif
