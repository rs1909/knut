
/* Possible C replacement for the UNIX FORTRAN routine ETIME */
/* for use by SHELXS and SHELXL; may need to append underscore */
/* if FORTRAN compiler requires it:  float etime_(et)  */

#include <stdio.h>
#include <time.h>

float etime_(et)
float et[2];
{
  clock_t clicks;
  clicks = clock();
  return ( (float) clicks / CLOCKS_PER_SEC );
}
