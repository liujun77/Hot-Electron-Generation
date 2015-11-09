#include "lib/jga.h"
#include "lib/complex.h"

numero eps_inf=9.5;
numero eps_g=0.071;
numero wp=9.06;

complex eps(numero w)
{
  return eps_inf-wp*wp/w/(w+i_c*eps_g);
}

numero ab(numero w){
  return w*imag((eps(w)-1)/(eps(w)+2));
}

int main()
{
  numero w;
  FILE *fout;
  fout=fopen("absorb.txt","w");
  for(w=2.6;w<=2.9;w+=0.002){
    fprintf(fout,"%g %g\n",w,ab(w));
  }
  fclose(fout);
  return 0;
}

