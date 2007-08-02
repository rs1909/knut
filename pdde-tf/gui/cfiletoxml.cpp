#include "constqtgui.h"


int main(int argc, const char** argv)
{
  // parameters
  NConstantsQtGui params;
  const char*     constFile = 0;
  const char*     xmlFile = 0;
  const char*     outFile = 0;
  const char*     inFile = 0;
  const char      usage[] = "\tUsage: cfiletoxml -c cfile -i input.mat -o output.mat -x cfile.xml -b branchfile\n";

  // argument parsing
  for (int acnt = 1; acnt < argc;  acnt++)
  {
    if (argv[acnt][0] == '-')
    {
      switch (argv[acnt][1])
      {
        case 'c':
          constFile = argv[++acnt];
          break;
        case 'i':
          inFile = argv[++acnt];
          break;
        case 'o':
          outFile = argv[++acnt];
          break;
        case 'x':
          xmlFile = argv[++acnt];
        case 'b':
          ++acnt;
          break;
        default:
          std::cout << usage;
          exit(0);
          break;
      }
    }
    else
    {
      std::cout << usage;
      exit(0);
    }
  }
  if (constFile)
  {
    params.loadFile(std::string(constFile));
//     params.printFile(std::cout);
    if (inFile) params.setInputFileText(std::string(inFile));
    if (outFile) params.setOutputFileText(std::string(outFile));
    if (xmlFile) params.saveXmlFile(std::string(xmlFile));
    else params.saveXmlFile(std::string(constFile)+".xml");
  }else std::cout << usage;
  return 0;
}
