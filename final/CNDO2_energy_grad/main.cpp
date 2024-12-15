#include <stdlib.h>
#include <stdexcept>
#include <stdio.h>
#include <armadillo>
#include <vector>
#include "AO.h"
#include "CNDO.h"

using namespace std;

int main(int argc, char *argv[])
{

  if (argc != 3)
  {
    printf("usage hw5 filename basis_dir, for example ./hw5 example.txt basis/\n");
    return EXIT_FAILURE;
  }
  string fname(argv[1]);
  string basis_dir(argv[2]);
  try
  {
    Molecule_basis mol(fname, basis_dir);

    CNDO ourSCF(mol, 50, 1e-5);
    int ok = ourSCF.init();
    if(ok != 0) return EXIT_FAILURE;
    ok = ourSCF.run();
    if(ok != 0) return EXIT_FAILURE;


    double Energy = ourSCF.getEnergy();
    arma::mat gradient = ourSCF.getGradient();

    std::cout << "Energy" << std::endl;
    std::cout << Energy << std::endl;
    gradient.print("Gradient");

  }
  catch (invalid_argument &e)
  {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

