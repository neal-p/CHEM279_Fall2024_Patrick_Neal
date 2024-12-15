#include <iostream>
#include "Integrands/Gaussian.hpp"
#include "Integrals/Numerical.hpp"
#include <vector>


int main(int argc, char** argv) {

  if (argc != 2) {
    std::cout << "You must provide an input file!" << std::endl;
    std::cout << "Usage: $ problem1 <input gauss file>" << std::endl;
    exit(1);
  }

  std::vector<Gauss> gaussians = ReadGauss(argv[1]);

  for (Gauss g : gaussians) {
    std::cout << "read gaussian centered at " << g.x0 << " with alpha=" << g.alpha << " and l=" << g.l << std::endl;
  }

  GaussProd gs(gaussians);
  std::cout << "1d numerical overlap integral between Gaussian functions is " << QSIMP(gs, 10.0, 1e-10) << std::endl;

  return 0;
}
