/* isize_get_lsode */
/* get size of int */
 
#include <stdio.h>
 
/* function names link differently depending on OS */
 
#include "f77name.h"
 
int F77NAME(isize_get_lsode_r8 )(idum)
     int *idum;   /* integer form -- ascii character to write */
{
  return sizeof(int);
}
