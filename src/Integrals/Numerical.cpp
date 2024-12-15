#include "Integrals/Numerical.hpp"

class Trapzd {

public:
  int n;
  double a, b, s;

  Trapzd(const double a, const double b) : a(a), b(b), n(0) {};

  double next(Integrand& I) {

    double x, tnm, sum, del;
    int it, j;
    this->n++;

    if (n == 1) {
      return (this->s = 0.5 * (this->b - this->a) * (I.eval(a) + I.eval(b)));

    } else {
      for (it=1, j=1; j < n-1; j++) it <<= 1;

      tnm = it;
      del = (this->b - this->a) / tnm;
      x = this->a + 0.5 * del;

      for (sum=0.0, j=0; j < it ; j++, x+= del) sum += I.eval(x);
      s = 0.5 * (this->s + (this->b - this->a) * sum / tnm);
      return s;
    }
  }
};



double QSIMP(Integrand& I, const double r, double tol) {

  const int JMAX=20;
  double s, st, ost=0.0, os=0.0;
  Trapzd t(I.x0 - r, I.x0 + r);

  for (int j=0; j < JMAX; j++) {

    st = t.next(I);
    s = (4.0 * st - ost) / 3.0;

    if (j > 5) {
      if (
           abs(s - os) < tol * abs(os) ||
           (s == 0.0 && os == 0.0)
         ) {

        return s;
      }

      os = s;
      ost = st;
    }
  }

  return s;

}
