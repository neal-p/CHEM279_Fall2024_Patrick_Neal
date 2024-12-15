#include <iostream>
#include "Integrands/Gaussian.hpp"


int main(int argc, char** argv) {

  if (argc != 2) {
    std::cout << "You must provide an input file!" << std::endl;
    std::cout << "Usage: $ problem1 <input gauss file>" << std::endl;
    exit(1);
  }

  std::vector<Shell3D> shells = ReadShells(argv[1]);

  for (Shell3D s : shells) {

    std::cout << "Shell:" << std::endl;
    std::cout << "  alpha=" << s.alpha << std::endl;
    std::cout << "  xyz=" << s.xyz.transpose() << std::endl;
    std::cout << "  L=" << s.L << std::endl;
    std::cout << "  Composed of: " << std::endl;

    for (Array3i lmn : s.lmns) {

      std::cout << "        lmn: " << lmn.transpose() << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::endl;
  }

  ArrayXXd overlap = ShellOverlap(shells[0], shells[1]);

  std::cout << "Overlap matrix:" << std::endl;
  std::cout << overlap << std::endl;

  return 0;
}

