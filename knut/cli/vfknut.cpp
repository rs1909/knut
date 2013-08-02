#define KNUTSYS_H 

#include "expr.h"
#include "exprsystem.h"
#include "matrix.h"
#include "knerror.h"

#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{
  if (argc != 2)
  {
    std::cerr << "Bad number of arguments.\n";
    exit(-1);
  }
  try
  {
    KNSystem sys (argv[1], false);

    std::string cxx;
    sys.toCxx (cxx);
    std::cout << cxx << "\n";
  //   std::ofstream file("dd.cpp");
  //   file << cxx;
  }
  catch (KNException& ex)
  {
    std::cout << ex.str();
  }
  return 0;
}
