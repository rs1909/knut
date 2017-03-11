#define KNUTSYS_H

#include "expr.h"
#include "exprsystem.h"
#include "matrix.h"
#include "knerror.h"

#include <iostream>
#include <fstream>
#include <fenv.h>

int main(int argc, char **argv)
{
#ifndef __APPLE__
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

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
    KNExprSystem sys (contents, false);

    std::string cxx;
    sys.toCxx (cxx);
// #ifndef DEBUG
    std::cout << cxx << "\n";
// #endif
//     cxx.clear ();
//     if (ExpTree::Expression::fromXML (cxx, contents))
//     {
//       std::cout << "/*\n";
//       std::cout << cxx;
//       std::cout << "\n*/\n";
//     }
  }
  catch (KNException& ex)
  {
    std::cout << ex.exprStr (contents);
  }
#ifdef DEBUG
  for (auto it : ExpTree::Node::instances)
  {
    if (it != nullptr) { it->print(std::cout); std::cout << "\t"; }
    else std::cout << "NULL\t";
  }
  std::cout << "\n";
  std::cout << "Node::instances = " << ExpTree::Node::instances.size() << "\n";
#endif
  return 0;
}
