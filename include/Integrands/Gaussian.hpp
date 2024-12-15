#pragma once

#include "Integrands/Integrands.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include "Simulation/Simulation.hpp"
#include <vector>
#include "eigen-3.4.0/Eigen/Dense"
#include "Math/Math.hpp"


class Gauss : public Integrand {

public:
  double x0, alpha;
  int l;

  Gauss(double x0, double alpha, int l) : x0(x0), alpha(alpha), l(l) {}
  double eval(double x);

};

class GaussProd : public Integrand  {

public:

  std::vector<double> coefs;
  std::vector<Gauss> gaussians;

  GaussProd(Gauss g1, Gauss g2);
  GaussProd(GaussProd g1s, Gauss g2);
  GaussProd(GaussProd g1s, GaussProd g2s);
  GaussProd(std::vector<Gauss> gaussians);

  GaussProd(Gauss g1, double c1, Gauss g2, double  c2);
  GaussProd(GaussProd g1s, double c1, Gauss g2, double c2);
  GaussProd(GaussProd g1s, double c1, GaussProd g2s, double c2);
  GaussProd(std::vector<Gauss> gaussians, std::vector<double> coeffs);

  double eval(double x);
};


std::vector<Gauss> ReadGauss(std::string filename);


class Shell3D {

public:
  Array3d xyz;
  double alpha;
  int L;

  std::vector<Array3i> lmns;

  Shell3D(Array3d xyz, double alpha, int L) : xyz(xyz), alpha(alpha), L(L) {
    GetFunctions();
  };
  Shell3D(double x, double y, double z, double alpha, int L) : alpha(alpha) , L(L) {

    this->xyz = {x, y, z};
    GetFunctions();
    }

  void GetFunctions();

  double eval(Array3d& xyz);
  double eval(double x, double y, double z);
};

std::vector<Shell3D> ReadShells(std::string filename);

ArrayXXd ShellOverlap(Shell3D& s1, Shell3D& s2);
