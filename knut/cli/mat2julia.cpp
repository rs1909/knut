#include "mat4data.h"
#include <fstream>

int main(int argc, char **argv)
{
  std::ofstream ofs;
  if (argc == 3)
  {
    ofs.open (argv[2], std::ofstream::out);
    auto data = KNDataFile (argv[1]);
    data.exportToJulia (ofs);
  } 
  else if  (argc == 2)
  {
    auto data = KNDataFile (argv[1]);
    data.exportToJulia (std::cout);
  } else
  {
    std::cerr << "Bad number of arguments.\n";
    exit(-1);
  }
  return 0;
}

