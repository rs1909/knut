#include "knerror.h"
#include "config.h"

void knutException::removePath()
{
  std::string::size_type loc = file.find(KNUT_SOURCE_DIR);
  if (loc != std::string::npos)
  {
    file.erase(loc, loc + strlen(KNUT_SOURCE_DIR));
  }
  if (file[0] = '/') file.erase(0,1);
}
