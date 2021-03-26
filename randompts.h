#ifndef RANDOMPTS_H
#define RANDOMPTS_H

#include<cmath>
#include<vector>
#include<array>
#include<random>
#include<stdio.h>
#include<math.h>
#include"point.h"

namespace wr
{

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0,1.0);
std::uniform_real_distribution<double> uniformDistribution(0.0,1.0);

double GetUniform()
{
    return uniformDistribution(generator);
}

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

}

#endif /* RANDOMPTS_H */
