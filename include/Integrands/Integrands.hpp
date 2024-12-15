#pragma once
#include "eigen-3.4.0/Eigen/Dense"
using namespace Eigen;

class Integrand {

public:
  double x0, alpha;
  int l;

  virtual double eval(double x) = 0;

};

class Integrand3D {

public:
  Array3d xyz;
  double alpha;

  virtual double eval(Array3d& xyz) = 0;
  virtual double eval(double x, double y, double z) = 0;

};
