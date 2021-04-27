// Copyright (c) Radost Waszkiewicz - 2021
// This code is licensed under MIT license
// Produces periodic cubic splines interpolation
// Algebra follows:
// https://kluge.in-chemnitz.de/opensource/spline/

#ifndef SPLINES_H
#define SPLINES_H

#include <cassert>
#include <vector>
#include <cmath>
#include <math.h>
#include <iostream>
#include "../external/Eigen/Dense"

namespace sp
{

    class Spline
    {
        private:
            std::vector<double> a;
            std::vector<double> b;
            std::vector<double> c;
            std::vector<double> y;

            Eigen::LLT<Eigen::MatrixXd> designLLT;
            Eigen::MatrixXd design; // Eigen matricies reference each other, need to keep the original.
            bool designLLTinitialized = false;

            void designFit(int n)
            {
                double dx = 1.0/n;
                design = Eigen::MatrixXd::Zero(n,n); // Eigen matricies need cleaning.
                for(int i=0; i<n; i++)
                {
                    design(i,i) = dx*4.0/3.0;
                    design(i,(i+1)%n) = dx*1.0/3.0;
                    design(i,(i-1+n)%n) = dx*1.0/3.0;
                }

                designLLT = Eigen::LLT<Eigen::MatrixXd>(design);
                designLLT.compute(design);
                designLLTinitialized = true;

                return;
            }
        public:
            void SetPoints(std::vector<double> yn)
            {
                y = yn;

                int n = yn.size();
                double dx = 1.0/n;

                if(!designLLTinitialized || designLLT.matrixL().rows() != yn.size())
                {
                    //printf("Designing\n");
                    designFit(n);
                }

                Eigen::VectorXd rhs(n);
                for(int i=0;i<n;i++)
                {
                    rhs(i) = (y[(n+i-1)%n] - 2.0*y[i] + y[(n+i+1)%n])/dx;
                }

                //std::cout<<rhs;
                auto sol = designLLT.solve(rhs);
                //std::cout<<sol;

                a = std::vector<double>(n);
                b = std::vector<double>(n);
                c = std::vector<double>(n);
                for(int i=0;i<n;i++)
                {
                    b[i] = sol(i);
                }
                for(int i=0;i<n;i++)
                {
                    a[i] = (b[(i+1)%n]-b[i])/(3.0*dx);
                    c[i] = (y[(i+1)%n]-y[i])/(dx) - dx*(2.0*b[i]+b[(i+1)%n])/(3.0);
                }
            }
            double operator() (double x) const
            {
                double q = fmod(fmod(x,1.0)+1.0,1.0);
                int n = y.size();
                int idx = q*n;

                double h = q - idx*(1.0/n);

                return ((a[idx]*h + b[idx])*h + c[idx])*h + y[idx];
            }
            double prime(double x) const
            {
                double q = fmod(fmod(x,1.0)+1.0,1.0);
                int n = y.size();
                int idx = q*n;

                double h = q - idx*(1.0/n);

                return 3.0*a[idx]*h*h + 2.0*b[idx]*h + c[idx];
            }
            double bis(double x) const
            {
                double q = fmod(fmod(x,1.0)+1.0,1.0);
                int n = y.size();
                int idx = q*n;

                double h = q - idx*(1.0/n);

                return 3.0*2.0*a[idx]*h + 2.0*b[idx];
            }
            double ter(double x) const
            {
                double q = fmod(fmod(x,1.0)+1.0,1.0);
                int n = y.size();
                int idx = q*n;

                double h = q - idx*(1.0/n);

                return 3.0*2.0*a[idx];
            }

    };

}
#endif /* SPLINES_H */
