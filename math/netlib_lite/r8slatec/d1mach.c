 /* C source for I1MACH */
 /* Note that some values may need changing -- see the comments below. */
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <math.h>
 
#include "f77name.h"   /* fortran linkage for i1mach, r1mach, d1mach */
 
 /* Standard C source for D1MACH -- remove the * in column 1 */
 
double F77NAME(d1mach)(int *i)
{
       switch(*i){
         case 1: return DBL_MIN;
         case 2: return DBL_MAX;
         case 3: return DBL_EPSILON/FLT_RADIX;
         case 4: return DBL_EPSILON;
         case 5: return log10((double)FLT_RADIX);
         }
       fprintf(stderr, "invalid argument: d1mach(%d)\n", *i);
       exit(EXIT_FAILURE); return 0; /* some compilers demand return values */
}
