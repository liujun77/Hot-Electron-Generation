// This program reads vik's files and remove the zeros

#include "lib/jga.h"
#include "lib/files.h"
#include "lib/algebra.h"

int main(int argc, char **argv)
{
 int i,j,N,nr,Nt;
 int *li=NULL;
 int *ni=NULL;
 numero *Energy=NULL;
 
 numero rin,rout,Rs,rmin,rmax,temp;
 
 char name_in[100],name_out[100];
 strcpy(name_in,argv[1]);    //Liu
 strcpy(name_out,argv[2]);   //Liu
 
 FILE *infile,*outfile;
 infile =fopen(name_in,"r"); 
 outfile=fopen(name_out,"w");
 printf("input the radius\n");
 scanf("%lg%lg",&rin,&rout); 
 //A=alread_numero(infile);   //Liu
 Rs=3;//alread_numero(infile);  //Liu
 N=alread_int(infile);
 rmin=alread_numero(infile);
 rmax=alread_numero(infile);  
 nr=alread_int(infile); 

 Energy=new numero[N]; li=new int[N]; ni=new int[N];

 Nt=0;  
 for(i=0; i<N; i++){
   li[i]     =alread_int(infile);
   ni[i]     =alread_int(infile);
   Energy[i] =alread_numero(infile);
   if(Energy[i]!=0) Nt++; 
 }
  
 fprintf(outfile,"%g %g %d %d %g %g %d\n",rin,rout,1,Nt,rmin,rmax,nr);
 for(i=0; i<N; i++){
  if(Energy[i]!=0) fprintf(outfile,"%d %d %.10g\n",li[i],ni[i],Energy[i]);
 }
 
 for(i=0; i<N; i++)
 for(j=0; j<nr; j++){
   temp=alread_numero(infile);
   if(Energy[i]!=0) fprintf(outfile,"%.10g\n",temp);
 }
 
 fclose(infile);   fclose(outfile);

 delete [] Energy; delete [] li; delete [] ni;
 return 0;

}
