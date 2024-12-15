#include "Math/Math.hpp"

int factorial(int x) {

  x = abs(x);

  int result = 1;

  for (int i=2; i <= x; i++) {
    result *= i;
  }

  return result;
}

int double_factorial(int x) {

  x = abs(x);
  int result = 1;

  if (x % 2 == 0) {
    for (int i=2; i <= x; i+=2) {
      result *= i;
    }
    
  } else {
    for (int i=3; i <= x; i+=2) {
      result *= i;
    }
  }

  return result;
}

double double_factorial(double x) {
  return (double)double_factorial((int)x);
}



int binom(int m, int n) {
  return factorial(m) / (factorial(n) * factorial(m - n));
}

double binom(double m, double n) {
  return (double)binom((int)m, (int)n);
}


