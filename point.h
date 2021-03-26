// Copyright (c) Radost Waszkiewicz - 2021
// This code is licensed under MIT license

#ifndef TRHEEDPOINT_H
#define TRHEEDPOINT_H

namespace wr
{

class Point
{
 public:
  double x;
  double y;
  double z;  
  double Norm()
  {
     return sqrt(x*x+y*y+z*z);
  }
  Point Normalized()
  {
     double n = Norm();
     Point ret = Point();
     ret.x = x/n; 
     ret.y = y/n; 
     ret.z = z/n;
     return ret;
  } 
  Point()
  {
    x=0;y=0;z=0;
  }
  Point operator+(Point rhs)
  {
    auto ret = Point();
    ret.x = (*this).x + rhs.x;
    ret.y = (*this).y + rhs.y;
    ret.z = (*this).z + rhs.z;
    return ret;
  }
  Point operator-(Point rhs)
  {
    auto ret = Point();
    ret.x = (*this).x - rhs.x;
    ret.y = (*this).y - rhs.y;
    ret.z = (*this).z - rhs.z;
    return ret;
  }
  Point operator^(Point rhs) // cross product
  {
    auto ret = Point();
    ret.x = rhs.z * (*this).y - rhs.y * (*this).z;
    ret.y = rhs.x * (*this).z - rhs.z * (*this).x;
    ret.z = rhs.y * (*this).x - rhs.x * (*this).y;
    return ret;
  }

  double operator*(Point rhs) // dot product
  {
    return rhs.x * (*this).x + rhs.y * (*this).y + rhs.z * (*this).z;
  }

  Point operator*(double rhs)
  {
    auto ret = Point();
    ret.x = (*this).x * rhs;
    ret.y = (*this).y * rhs;
    ret.z = (*this).z * rhs;
    return ret;
  }

  Point(double x_,double y_,double z_)
  {
    x = x_;
    y = y_;
    z = z_;
  }

};

}

#endif /* 3DPOINT_H */
