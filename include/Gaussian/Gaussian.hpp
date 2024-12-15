#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <iostream>

#include "Math/Math.hpp"
#include "PeriodicTable/PeriodicTable.hpp"

#include "../../external/eigen-3.4.0/Eigen/Dense"
using namespace Eigen;



class _Primative {

    public:
        double coef, alpha, N;
        _Primative(double coef, double alpha) : coef(coef), alpha(alpha) {};
        _Primative(const _Primative& p);
};

class AO {

    public:
        Array3d center;
        int owning_atom_idx;
        int l,m,n;
        std::vector<_Primative> primatives;

        AO(int l, int m, int n): l(l), m(m), n(n) {};
        AO(const AO& ao);

        void add_primative(double coef, double alpha);
        int L() {return this->l + this->m + this->n;};

};



typedef std::unordered_map<int, std::vector<AO>> BasisSet;
BasisSet read_basis_set(std::string file);
 
double PrimativeOverlap(_Primative& A, _Primative& B, Array3d A_center, Array3d B_center, int A_l, int B_l, int A_m, int B_m, int A_n, int B_n);
double AO_overlap(AO& A, AO& B);

double __inner_overlap(double alpha, double beta, double Ra, double Rb, int la, int lb);
double __inner_overlap_derivative(double alpha, double beta, double Ra, double Rp, int la, int lb);
