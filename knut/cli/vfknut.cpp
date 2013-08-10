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
  std::string contents;
  std::ifstream in (argv[1], std::ios::in | std::ios::binary);
  if (in)
  {
    in.seekg (0, std::ios::end);
    contents.resize (in.tellg());
    in.seekg (0, std::ios::beg);
    in.read (&contents[0], contents.size());
    in.close ();
  }
  try
  {
    KNSystem sys (argv[1], false);

    std::string cxx;
    sys.toCxx (cxx);
    std::cout << cxx << "\n";
    cxx.clear ();
    if (ExpTree::Expression::fromXML (cxx, contents))
    {
      std::cout << "/*\n";
      std::cout << cxx;
      std::cout << "\n*/\n";
    }
  }
  catch (KNException& ex)
  {
    std::cout << ex.str();
  }
  return 0;
}
