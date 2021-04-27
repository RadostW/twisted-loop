// Copyright (c) Radost Waszkiewicz - 2021
// This code is lincensed under MIT license
// Computes energy of twisted, bent loop of dna
// With energy density consisting of 5 parts
// 1. Length penalty (E = large*(L - L_0)^2)
// 2. Steric interactions (E = large*(intersection_length))
// 3. Bending energy (dE = EI * \kappa^2 ds)
// 4. Twisting energy (dE = \omega EI * \Omega^2 ds)
// 5. Electrostatic interaction (dE = const * exp(-\lambda distance) / distance)

#ifndef ENERGY_H
#define ENERGY_H

#include<cmath> // M_PI
#include"spline_utils.h" // TotalCurvature, NearNeighbours, DebeyeEnergy
#include"splines.h" // Spline
#include"points.h" // Curve
#include<vector> // std::vector

namespace lp
{
    class Loop
    {
        public:
            sp::Spline sx;
            sp::Spline sy;
            sp::Spline sz;
            Loop(pt::Curve c)
            {
                cx = std::vector<double>(c.size());
                cy = std::vector<double>(c.size());
                cz = std::vector<double>(c.size());
    
                for(int i=0;i<c.size();i++)
                {
                    cx[i] = c[i].x;
                    cy[i] = c[i].y;
                    cz[i] = c[i].z;
                }

                sx.SetPoints(cx);
                sy.SetPoints(cy);
                sz.SetPoints(cz);
            }
        private:
            std::vector<double> cx;
            std::vector<double> cy;
            std::vector<double> cz;

    };

    double LoopEnergy(Loop loop, 
                         double Lk, // Linking number
                         double omega, // Twisting to bending coef. ratio
                         double ElectroElasticNumber, // (N^2 q^2 / EI 4 \pi \epsilon)
                         int BasePairs,
                         double DebeyeLength, // Screening length for electrostatics
                         double StericPenalty,
                         double StericDiameter, // [multiples of beam length]
                         double LengthPenalty)
    {
        double Tw = Lk - sp::Writhe(loop.sx,loop.sy,loop.sz);
        double Omega = 2*M_PI*Tw;
        double Length =  sp::Length(loop.sx,loop.sy,loop.sz);
        double SquaredCurvatureIntegral = sp::TotalCurvature(loop.sx,loop.sy,loop.sz);
        double IntersectionLength = sp::NearNeighbours(loop.sx,loop.sy,loop.sz,StericDiameter);
        double DebeyeEnergy = sp::DebeyeEnergy(loop.sx,loop.sy,loop.sz,DebeyeLength,BasePairs);

        return LengthPenalty*(Length - 1)*(Length - 1) 
               + StericPenalty * IntersectionLength
               + SquaredCurvatureIntegral
               + omega*Tw*Tw 
               + ElectroElasticNumber * DebeyeEnergy;
    }

}

#endif
