#ifndef jga_legendre      // ************************************************
#define jga_legendre 1    // ***  jga/legendre.h              13-viii-98  ***
                          // ***                                          ***
#include "jga.h"          // ***                                          ***
                          // ***                                          ***
// ***  ----------------  // ***  --------------------------------------  ***
// ***                                                                    ***
// ***  Author: F. J. Garcia de Abajo, CSIC, Spain                        ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Legendre polynomials. Notation of Messiah: a factor (-1)^m        ***
// ***  different with respect to Abramowitz                              ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  y=legendre(l,m,x), for real and complex arguments x               ***
// ***                                                                    ***
// **************************************************************************


// **************************************************************************
// *** real arguments                                                     ***
// **************************************************************************

numero legendre(int l, int m, numero x)
{
  if(m<0 || m>l || x<-1 || 1<x) {
    printf(" Error: Bad argument in legendre: l=%1d  m=%1d  x=%f\n",l,m,x);
    return 0;
  }

  int  i,ll;
  numero  pmm=1, somx2, fct=1, pmmp1, pll;

  if(m>0) {
    somx2=sqrt((1-x)*(1+x));
    for(i=1; i<=m; i++, fct+=2)  pmm=pmm*fct*somx2;
  }
  if(l==m)    return pmm;
  pmmp1=x*(2*m+1)*pmm;
  if(l>m+1)  for(ll=m+2; ll<=l; ll++) {
               pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
               pmm=pmmp1;  pmmp1=pll;
             }
  return  pmmp1;
}


// **************************************************************************
// *** complex arguments                                                  ***
// **************************************************************************

#ifdef jga_complex

complex legendre(int l, int m, complex x)
{
  if(m<0 || m>l) {
    printf(" Error: Bad argument in legendre: l=%1d  m=%1d\n",l,m);
    return complex(0,0);
  }

  int  i,ll;
  complex  pmm=1, somx2, fct=1, pmmp1, pll;

  if(m>0) {
    somx2=sqrt((1-x)*(1+x));
    for(i=1; i<=m; i++, fct+=2)  pmm=pmm*fct*somx2;
  }
  if(l==m)    return pmm;
  pmmp1=x*(2*m+1)*pmm;
  if(l>m+1)  for(ll=m+2; ll<=l; ll++) {
               pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
               pmm=pmmp1;  pmmp1=pll;
             }
  return  pmmp1;
}

#endif  // complex

#endif  // ******************************************************************
