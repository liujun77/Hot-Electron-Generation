#define HOT_ELE_HELP                                                         "\
**************************************************************************  \n\
hot_ele.cpp        version 2.0beta                      14-x-13 ->          \n\
**************************************************************************  \n\
                                                                            \n\
   Hot-electron calculations for gold nanoparticles                         \n\
   Alejandro Manjavacas Jun Liu                                             \n\
                                                                            \n\
**************************************************************************  \n\
                                                                            \n\
Gaussian atomic units are used inside the code.                             \n\
                                                                            \n\
Input file:                                                                 \n\
                                                                            \n\
  - Material parameters must be defined first, then convergence parameters  \n\
    (if needed), calculation options and finally calculation loop.          \n\
  - simple algebraic expressions are allowed on input (without spaces),     \n\
  - use 'screen' as output to have the results displayed on the screen      \n\
                                                                            \n\
**************************************************************************  \n\
List of commands:                                                           \n\
                                                                            \n\
------ Material parameters -----------------------------------------------  \n\
                                                                            \n\
-->  A value                                                                \n\
Nanoparticle radius expressed in nm. Default: 1.5 nm.                       \n\
                                                                            \n\
-->  V0 value                                                               \n\
                                                                            \n\
Value of the potential at the bottom of the well in eV. Default: -10 eV.    \n\
                                                                            \n\
-->  Material value                                                         \n\
                                                                            \n\
1 for Gold                                                                  \n\
2 for Silver (not implemented yet)                                          \n\
                                                                            \n\
Default: 1 (Gold)                                                           \n\
                                                                            \n\
-->  T value                                                                \n\
                                                                            \n\
Assign the Temperature in K. Default: T=300 K.                              \n\
                                                                            \n\
-->  tau value                                                              \n\
                                                                            \n\
Hot carrier decay time in ps. Default: tau=0.5 ps.                          \n\
(hbar/tau = 0.00131642 eV)                                                  \n\
                                                                            \n\
------ Input/output of wavefunctions and energies ------------------------- \n\
                                                                            \n\
-->  Energies write/read output/input                                       \n\
                                                                            \n\
This command is needed to run calculations unless the command Wavefunctions \n\
with the option read is used (only nanoparticles)                           \n\
                                                                            \n\
Write/read the energies in output/from input file.                          \n\
                                                                            \n\
The format of the file is:                                                  \n\
                                                                            \n\
    l0 n0 energy_0                                                          \n\
    l0 n1 energy_1                                                          \n\
    .  .  .                                                                 \n\
    .  .  .                                                                 \n\
    .  .  .                                                                 \n\
    ln nn energy_n                                                          \n\
                                                                            \n\
The energies are ordered following increasing l and then n.                 \n\
Energies are in eV                                                          \n\
                                                                            \n\
-->  Wavefunctions write/read output/input                                  \n\
                                                                            \n\
This command is optional (if you choose to read you do not need to define   \n\
A, material, V0).                                                           \n\
                                                                            \n\
Write/read the wavefunctions and energies in output/from input file.        \n\
                                                                            \n\
The format of the file is:                                                  \n\
                                                                            \n\
    A Material N rmax nr                                                    \n\
    l0 n0 energy_0                                                          \n\
    l0 n1 energy_1                                                          \n\
    .  .  .                                                                 \n\
    .  .  .                                                                 \n\
    .  .  .                                                                 \n\
    ln nn energy_n                                                          \n\
    wavefunction_0(dr)                                                      \n\
    wavefunction_0(2*dr)                                                    \n\
    wavefunction_0(3*dr)                                                    \n\
    .                                                                       \n\
    .                                                                       \n\
    .                                                                       \n\
    wavefunction_n(n*dr)                                                    \n\
                                                                            \n\
for nanoparticles, while for nanoshells we have:                            \n\
                                                                            \n\
    rin rout Material N rmin rmax nr                                        \n\
    l0 n0 energy_0                                                          \n\
    l0 n1 energy_1                                                          \n\
    .  .  .                                                                 \n\
    .  .  .                                                                 \n\
    .  .  .                                                                 \n\
    ln nn energy_n                                                          \n\
    wavefunction_0(rmin+dr)                                                 \n\
    wavefunction_0(rmin+2*dr)                                               \n\
    wavefunction_0(rmin+3*dr)                                               \n\
    .                                                                       \n\
    .                                                                       \n\
    .                                                                       \n\
    wavefunction_n(rmin+n*dr)                                               \n\
                                                                            \n\
The energies and wavefunctions are ordered following increasing l           \n\
and then n.                                                                 \n\
                                                                            \n\
Important: All magnitudes are in a.u.                                       \n\
                                                                            \n\
                                                                            \n\
-->  Nanoshell on/off                                                       \n\
                                                                            \n\
When on the calculation is performed for a nanoshell. The wavefunctions     \n\
and energies must be provided using Wavefunctions with the option read.     \n\
                                                                            \n\
Default: off.                                                               \n\
                                                                            \n\
------ Scanning parameters -----------------------------------------------  \n\
                                                                            \n\
-->  w nw val1_1 val2_1 n_1 ... val1_nw val2_nw n_nw                        \n\
                                                                            \n\
Assign nw ranges of frequencies in w (n_i uniformly distributed values from \n\
val1_i to val2_i, all in eV). For a range containing a single frequency,    \n\
use the assignment 'val1_i val1_i 1'.                                       \n\
                                                                            \n\
-->  grid_xz  x0 x1 nx z0 z1 nz                                             \n\
                                                                            \n\
Real space grid used to calculate the carrier spatial distribution:         \n\
from x0 to x1 and from z0 to z1 (in units of the particle radius)           \n\
with a total of nx and nz equally-spaced points.                            \n\
                                                                            \n\
------ Actual calculations -----------------------------------------------  \n\
                                                                            \n\
-->  begin-calculation                                                      \n\
-->  end-calculation                                                        \n\
                                                                            \n\
Different calculations can be made inside a begin ... end loop, so that     \n\
the system needs to be initialized only once for all of them. The           \n\
following commands define the calculations that can be actually made.       \n\
                                                                            \n\
-->  calc-e-spectrum output                                                 \n\
                                                                            \n\
Calculate the spectrum of hot electrons.                                    \n\
The output contains the following columns:                                  \n\
                                                                            \n\
       column     1:  w/eV                                                  \n\
       column     2:  Ef/eV                                                 \n\
       column     3:  Number of hot electrons with energy Ef per unit       \n\
                      time/a.u.                                             \n\
                                                                            \n\
-->  calc-h-spectrum output                                                 \n\
                                                                            \n\
Calculate the spectrum of hot holes.                                        \n\
The output contains the following columns:                                  \n\
                                                                            \n\
       column     1:  w/eV                                                  \n\
       column     2:  Ef/eV                                                 \n\
       column     3:  Number of hot holes with energy Ef per unit           \n\
                      time/a.u.                                             \n\
                                                                            \n\
-->  calc-num-spectrum output                                               \n\
                                                                            \n\
Calculate the total number of hot electrons.                                \n\
                                                                            \n\
The output contains the following columns:                                  \n\
                                                                            \n\
       column     1:  w/eV                                                  \n\
       column     2:  Total number of hot electrons per unit time/a.u.      \n\
                                                                            \n\
-->  calc-carrier-distribution output                                       \n\
                                                                            \n\
Calculate the spatial hot electron distribution at the positions defined    \n\
by the grid command which must be invoked before the calculation starts.    \n\
The output contains the following columns:                                  \n\
                                                                            \n\
       column     1:  w/eV,                                                 \n\
       column     2:  x/nm,                                                 \n\
       column     3:  z/nm,                                                 \n\
       column     4:  electron density/nm^3.                                \n\
       column     5:  hole density/nm^3.                                    \n\
                                                                            \n\
-->  calc-carrier-distribution-E E0 E1 output                               \n\
                                                                            \n\
Similar to calc-carrier-distribution but it only calculates the spatial     \n\
distribution of the carriers whose energies lie in the range E0 E1.         \n\
E0 and E1 must be defined in eV.                                            \n\
The output contains the following columns:                                  \n\
                                                                            \n\
       column     1:  w/eV,                                                 \n\
       column     2:  x/nm,                                                 \n\
       column     3:  z/nm,                                                 \n\
       column     4:  electron density/nm^3.                                \n\
       column     5:  hole density/nm^3.                                    \n\
                                                                            \n\
------ Convergence parameters --------------------------------------------- \n\
                                                                            \n\
-->  EPS value                                                              \n\
                                                                            \n\
Precision for the eignvalue solver. Default: 1e-10.                         \n\
                                                                            \n\
-->  de value                                                               \n\
                                                                            \n\
Energy step for the eignvalue solver in eV. Default: 1e-3 eV                \n\
                                                                            \n\
-->  nr value                                                               \n\
                                                                            \n\
Number of points for r integrals. Default: 1000.                            \n\
                                                                            \n\
-->  rmax value                                                             \n\
                                                                            \n\
Upper limit for the r integrals in units of A. Default: 5 (i.e. rmax=5A).   \n\
                                                                            \n\
-->  lmax value                                                             \n\
                                                                            \n\
Maximum value of l used in the initialization of the factors involving      \n\
factorials. Default: 400                                                    \n\
                                                                            \n\
------ Optional reporting, miscellaneous, and end of script --------------  \n\
                                                                            \n\
-->  time on|off                                                            \n\
                                                                            \n\
When turned on, information on total/partial execution times is given on a  \n\
file named as the input script-file followed by \".time\". Default: off.    \n\
                                                                            \n\
-->  end                                                                    \n\
                                                                            \n\
End the execution of the program. This command is necessary because         \n\
several calculations can be done within a single input file.                \n\
                                                                            \n\
**************************************************************************  \n"

#include "lib/jga.h"
#include "lib/complex.h"
#include "lib/qfint.h"
#include "lib/spline.h"
#include "lib/bessel_calculation.h"
#include "lib/legendre.h"
#include "lib/time.h"
#include "lib/algebra.h"
#include "lib/files.h"
//#include <omp.h>

// --- Material parameter ------------------------------------------------

numero A=30*nm;                 //Nanoparticle radius expressed in nm
numero V0=-10/au_eV;            //Value of the potential at the bottom of the
                                //well in eV

int Material=1;                 //1 for Gold

numero Rs;                      //Free electron radius expressed in Bohr radius
numero eps_inf;                 //eps inifinity Drude model
numero eps_g;                   //dephasing rate Drude model
                                //keep the decay rate of the plasmon
numero wp;                      //Gold plasma frequency
numero kb_au=3.1668151e-06;     //kB in 1/(au K)
numero KT=kb_au*300;            //Material temperature 
numero au_s=2.418884326502e-17; //1 au of time in sec
numero tau=0.5e-12/au_s;        //hot carrier dacay time
// -----------------------------------------------------------------------

numero EPS=1e-8;         //Precision for the eignvalue solver
numero de=0.001/au_eV;   //Energy step for the eignvalue solver
int    nr=1000;          //Number of points for the r integrals
numero rmax=5*A;         //Upper limit for the r integrals in units of A
int lmax=400;            //Upper limit for fact1 and fact2 initiated

// -----------------------------------------------------------------------

// Nanoshell comands
int nanoshell_flag=0;    //Flag for nanoshell calculations
numero rin=0;            //Inner nanoshell radius in nm
numero rout=0;           //Outer nanoshell radius in nm
numero rmin=0;           //Lower limit for the r integrals in units of rin

// -----------------------------------------------------------------------

numero  *aa;             //Vector to keep the normalized coeff wavefunc (r<A)
complex *bb;             //Vector to keep the normalized coeff wavefunc (r>A)
numero  *occ=NULL;       //Vector to keep the occupation of state i 
int     *li=NULL;        //Vector to keep the value of l for state i
int     *ni=NULL;        //Vector to keep the value of n for state i
numero  *Energy=NULL ;   //Vector to keep the energies
spline  *State=NULL;     //Vector to keep the wavefunction
                         //(when reading from a file)
int     *pos=NULL;       //Vector to keep the position once the energies are
                         //ordered from the lowest to the highest
int     N;               //Number of states (without spin and m degeneracy)
int     num_elec;        //Number of electrons
numero  FERMI;           //Fermi level of the system
numero  FERMI_fc;        //Factor to accoun for the possible degeneration of the
                         //Fermi level
int grid_flag=0;         //Flag for the xy grid
int energies_flag=0;     //Flag for Energies reading/writting
                         //1 calculates & writes the energies (free elec model)
                         //2 reads the energies
int wavefunctions_flag=0;//Flag for Wavefunctions reading/writting
                         //0 calculates wavefunctions with free electron model
                         //1 same as 0 but also writes the wavefunctions
                         //2 reads the wavefunctions

numero  *spec_e=NULL ;   //Vector for hot electron spectrum
                         //(used in init_carrier_distribution)
numero  *spec_h=NULL ;   //Vector for hot hole spectrum
                         //(used in init_carrier_distribution)

int n_core;              //Number of cores (parallel calc)

// --- time tracking -----------------------------------------------------

char *tname=new char[105];         // file to write real CPU
FILE *tfile;
int  tflag=0,tflag_1st=1;

void hot_ele_time(char *name)
{
  if(tflag) {
    if(tflag_1st) 
      tfile=fopen(tname,"w"); 
    else 
      tfile=fopen(tname,"a");
    fprintf(tfile,"%s",name);  time_print(tfile);  fclose(tfile);  tflag_1st=0;
} }

void hot_ele_time(char *name, numero xxx)
{
  if(tflag) {
    if(tflag_1st) 
      tfile=fopen(tname,"w"); 
    else 
      tfile=fopen(tname,"a");
    fprintf(tfile,"%s%g\n",name,xxx);  fclose(tfile);  tflag_1st=0;
} }

// --- Function to perform integrations ----------------------------------

long numero integ(numero a, numero b, long numero f(numero x), int n)
{
 if(n<7)  n=7;
 int i;  numero  sum=0, h=(b-a)/n;
 if(a==b)  return 0;
 
 sum=(  17 * (f(a)     + f(b    ))  + 59 * (f(a+h)   + f(b-h  ))
      + 43 * (f(a+2*h) + f(b-2*h))  + 49 * (f(a+3*h) + f(b-3*h))) / 48;
 for(i=4; i<=n-4; i++)  sum+=f(a+i*h);
 
 return sum*h;
}

numero integ(numero a, numero b, numero f(numero x), int n)
{
 if(n<7)  n=7;
 int i;  numero  sum=0, h=(b-a)/n;
 if(a==b)  return 0;
 
 sum=(  17 * (f(a)     + f(b    ))  + 59 * (f(a+h)   + f(b-h  ))
      + 43 * (f(a+2*h) + f(b-2*h))  + 49 * (f(a+3*h) + f(b-3*h))) / 48;
 for(i=4; i<=n-4; i++)  sum+=f(a+i*h);
 
 return sum*h;
}

// --- Set Material parameters -------------------------------------------

void set_Material(void)  //Define material parameters
{
  if(Material==1){   //Gold                        
    Rs=3.0;               
    eps_inf=9.5;          
    eps_g=0.071/au_eV;    
    wp=9.06/au_eV;       
  }
  if(Material==2){   //Silver (not implemented)
    Rs=3.0;               
    eps_inf=9.5;          
    eps_g=0.071/au_eV;    
    wp=9.06/au_eV;       
  }
}

// --- Drude model dielectric function -----------------------------------

complex eps(numero w) //eps = 1-9.06^2/omega^2
{
  return eps_inf-wp*wp/w/(w+i_c*eps_g);
}

// --- Energy calculation ------------------------------------------------
FILE *Energies_file;
char *Energies_name=new char[105];  // file to write/read Energy

numero func_energy(int l, numero energy) //Function used to calculate the energy
{
  numero alfa, k;
  complex val;
  alfa=sqrt(2*(energy-V0));
  k=sqrt(-2*energy);
  val=alfa*besseljp_j(l,alfa*A)-i_c*k*besselh1ip_h(l,k*A);
  return real(val);
}

numero search_root(numero x1, numero x2, int l) //search the root between x1 &
{                                               //x2 of function func_energy
  numero mid,fx1,fx2,fmid;
  do{
    mid=(x1+x2)/2;
    fx1=func_energy(l,x1);  
    fx2=func_energy(l,x2);
    fmid=func_energy(l,mid);
    if(fx1*fx2>0||(fabs((fx1-fx2)/(x1-x2))>1e6&&fmid>10))
    return 1;                //avoid the divergence case means there is no root
    else{ 
      if(fx1*fmid>0) x1=mid;
      else x2=mid;
    }
  }while(fabs(fmid)>EPS);
  return mid;
}

void get_free_energies(void)
{
  int l=-1,i,nnl=1,ii=0; N=0;
  numero y1,y2,energy,root;

  if(energies_flag==0)
  on_error("Energies",
           "You must specify if you want to write or read the energies.");

  if(energies_flag==1){ // Calculates and writes the energies
    Energies_file=fopen(Energies_name,"w");
    while(nnl!=0){
      nnl=0; l++;
      y1=func_energy(l,V0+de);
      for(energy=V0+de;energy<-de;energy+=de){
        y2=func_energy(l,energy+de);
        if(y1*y2<0){
          root=search_root(energy,energy+de,l);
          if(root<0){
            fprintf(Energies_file,"%d %d %.10g\n",l,nnl,root*au_eV);
            nnl++; N++;
        } }
       y1=y2;  
    } }
    fclose(Energies_file);                                          
  }
  
  Energies_file=fopen(Energies_name,"r");
  if(Energies_file==NULL) on_error("Energies",
                                   "cannot open energies file:", Energies_name);
  N=number_of_rows(0,Energies_name); //Calculates the number of rows
  Energy=new numero[N]; li=new int[N]; ni=new int[N];

  for(i=0; i<N; i++){
   li[i]    =alread_int(Energies_file);
   ni[i]    =alread_int(Energies_file);
   Energy[i]=alread_numero(Energies_file)/au_eV; 
  } 
  fclose(Energies_file);
}

// --- Normalization of radial wavefunction ------------------------------
int inn;

long numero f1norm(numero r)
{
 numero alfa=sqrt(2*(Energy[inn]-V0));
 if(r==0) return 0;
 else return r*r*besselj_n(li[inn],alfa*r)*besselj_n(li[inn],alfa*r);
}

long numero f2norm(numero r)
{
 numero k=sqrt(-2*Energy[inn]);
 long numero auxr=real(besselh1i(li[inn],k*r));
 long numero auxi=imag(besselh1i(li[inn],k*r));
 return r*r*(auxr*auxr+auxi*auxi);
}

void normalize_coeff(void)
{
  aa=new numero[N];   bb=new complex[N];
  long numero inte1, inte2, inte;
  numero k, alfa;
  for(inn=0;inn<N;inn++){
    alfa=sqrt(2*(Energy[inn]-V0));
    k=sqrt(-2*Energy[inn]);
    inte1=integ(0,   A,f1norm,nr);
    inte2=integ(A,rmax,f2norm,nr);
    inte=inte1+mod2(besselj_n(li[inn],alfa*A)/besselh1i(li[inn],k*A))*inte2;
    aa[inn]=sqrt(1/inte);
    bb[inn]=besselj_n(li[inn],alfa*A)/besselh1i(li[inn],k*A)*aa[inn];
  }
}

// --- Calculate the Fermi level -----------------------------------------

void set_Fermi(void) // Sets the Fermi level and the Fermi level degeneration
{
  int i,j,imax,ii;
  numero temp_energy; 
  int temp_li,temp_ni,temp_pos,num_temp;

  pos=new int[N]; for(i=0; i<N; i++) pos[i]=i;
  
  num_elec=int(A*A*A/Rs/Rs/Rs);       // Total number of electrons
 
  for(i=0; i<N; i++)
  for(j=0; j<N; j++) 
  if(Energy[j]>Energy[i]){
    temp_energy=Energy[i]; temp_li=li[i]; temp_ni=ni[i];  temp_pos=pos[i];
    Energy[i]=Energy[j];   li[i]=li[j];   ni[i]=ni[j];    pos[i]=pos[j];
    Energy[j]=temp_energy; li[j]=temp_li; ni[j]=temp_ni;  pos[j]=temp_pos;
  }

  num_temp=num_elec;
  i=-1;
  while(num_temp>0){
   i++;
   num_temp-=2*(2*li[i]+1); // 2 for spin and 2l+1 for m-degeneration
   FERMI=Energy[i];
  }
  imax=i; num_temp=0;
  for(i=0; i<imax; i++) num_temp+=2*(2*li[i]+1);
  FERMI_fc=(num_elec-num_temp)/2./(2.*li[imax]+1.);
}

numero f_FERMI(numero E, numero KT) //Fermi-Dirac distribution function
{
  if(E==FERMI) return FERMI_fc;
  if(ABS(KT)<1e-20) if(E>FERMI) return 0; else return 1;
  if((E-FERMI)/KT>50) return 0;  if((E-FERMI)/KT<-50) return 1;
  return 1/(exp((E-FERMI)/KT)+1);
}

// --- Radial wavefunction -----------------------------------------------

numero psiR(int i, numero r)
{
  if(wavefunctions_flag==2)  return State[pos[i]].val(r);
  else{
    numero alfa,k;
    alfa=sqrt(2*(Energy[i]-V0));
    k=sqrt(-2*Energy[i]);
    if(r<A) return aa[i]*besselj_n(li[i],alfa*r);
    else    return real(bb[i]*besselh1i(li[i],k*r));
  }
}

// --- Get Wavefunctions and Energies ---------------------------------

FILE *Wavefunctions_file;
char *Wavefunctions_name=new char[105];  // file to write/read Wavefunctions

void get_wavefunctions(void)
{
  numero dr,r,val;
  int i,j,nrr;
  if(wavefunctions_flag==0||wavefunctions_flag==1){  //Calculate wavefunctions
                                                     //with free electron model
    dr=rmax/nr;    

    hot_ele_time("Calculation of energies:                  ");
    get_free_energies(); 				               
    hot_ele_time("                                          "); 

    set_Material();

    hot_ele_time("Calculation of Fermi energy:              ");
    set_Fermi();
    hot_ele_time("                                          ");

    hot_ele_time("Calculation of normalization coeff:       ");
    normalize_coeff();
    hot_ele_time("                                          ");
   
    if(wavefunctions_flag==1) {                     // write the Wavefunctions
      Wavefunctions_file=fopen(Wavefunctions_name,"w");
      fprintf(Wavefunctions_file,"%g %d %d %g %d\n",A,Material,N,rmax,nr);
      for(i=0; i<N; i++)
      fprintf(Wavefunctions_file,"%d %d %.10g\n",li[i],ni[i],Energy[i]);
      for(i=0; i<N; i++) 
      for(j=1; j<=nr; j++){
        r=j*dr;
        fprintf(Wavefunctions_file,"%.10g\n",psiR(i,r));
      }
      fclose(Wavefunctions_file);
    }
  }

  if(wavefunctions_flag==2) {                       // read the wavefunctions
    Wavefunctions_file=fopen(Wavefunctions_name,"r");
    if(Wavefunctions_file==NULL)
     on_error("Wavefunctions", "cannot open wavefunctions file:"
              , Wavefunctions_name);

    if(nanoshell_flag==0){ //nanosphere
      A       =alread_numero(Wavefunctions_file); 
      Material=alread_int   (Wavefunctions_file); 
      N       =alread_int   (Wavefunctions_file); 
      rmax    =alread_numero(Wavefunctions_file); 
      nrr     =alread_int   (Wavefunctions_file);
    }
    if(nanoshell_flag==1){ //nanoshell
      rin     =alread_numero(Wavefunctions_file); 
      rout    =alread_numero(Wavefunctions_file); 
      Material=alread_int   (Wavefunctions_file); 
      N       =alread_int   (Wavefunctions_file); 
      rmin    =alread_numero(Wavefunctions_file); 
      rmax    =alread_numero(Wavefunctions_file); 
      nrr     =alread_int   (Wavefunctions_file);
    }
    
    dr=(rmax-rmin)/nrr;
    set_Material();
    Energy=new numero[N]; li=new int[N]; ni=new int[N];
    State=new spline[N];
   
    for(i=0; i<N; i++){
       li[i]=alread_int(Wavefunctions_file);
       ni[i]=alread_int(Wavefunctions_file); 
       Energy[i]=alread_numero(Wavefunctions_file);
    }
  
    for(i=0; i<N; i++){
      State[i].alloc(nrr); 
      for(j=1; j<=nrr; j++){
       r=rmin+j*dr;
       val=alread_numero(Wavefunctions_file);
       State[i].put(j-1,r,val/r);  //AM we divide by r because Vik's
      }                            //wavefunctions are u(r)=rR(r)
      State[i].init();
    }
    
    fclose(Wavefunctions_file);
   
    hot_ele_time("Calculation of Fermi energy:              ");
    set_Fermi();
    hot_ele_time("                                          ");
  }
}

// --- Angular transition matrix element ---------------------------------

numero tran_ang(int l, int lp, int m) //|<lm|cos(theta)|l'm'>|^2
{
  numero val=0; int lmin;
  lmin=(l<lp)?l:lp;
  if(abs(m)<=lmin){
    if(lp==l+1) val=1.0*(l-m+1)*(l+m+1)/(2*l+3)/(2*l+1);
    if(lp==l-1) val=1.0*(l-m  )*(l+m  )/(2*l-1)/(2*l+1);
  }
  return val;
}

// --- Radial transition matrix element ----------------------------------
numero *tran_V1=NULL; //keep the integral <i|r|j> r<A           for nanospheres
                      //or                <i|r|j> r<rin         for nanoshells
numero *tran_V2=NULL; //keep the integral <i|r|j> r>A           for nanospheres
                      //or                <i|r|j> rin<r<rout    for nanoshells
numero *tran_V3=NULL; //keep the integral <i|r^-2|j> r>A        for nanospheres
                      //or                <i|r^-2|j> rin<r<rout for nanoshells
numero *tran_V4=NULL; //keep the integral <i|r|j>    r>rout     for nanoshells
numero *tran_V5=NULL; //keep the integral <i|r^-2|j> r>rout     for nanoshells

int it,jt;

numero f1tran(numero r)
{
 if(r==0) return 0;
 else return r*r*r*psiR(it,r)*psiR(jt,r);
}

numero f2tran(numero r)
{
 return psiR(it,r)*psiR(jt,r);
}

void set_tran_r(void)
{
  tran_V1=new numero[N*N]; tran_V2=new numero[N*N]; tran_V3=new numero[N*N];
  if(nanoshell_flag){      tran_V4=new numero[N*N]; tran_V5=new numero[N*N];}
  
  if(nanoshell_flag==0) //Nanospheres
  for(it=0;it<N;it++)
  for(jt=0;jt<N;jt++){
    if(li[it]==li[jt]+1||li[it]==li[jt]-1){
      tran_V1[it*N+jt]=integ(0,   A,f1tran,nr);
      tran_V2[it*N+jt]=integ(A,rmax,f1tran,nr);
      tran_V3[it*N+jt]=integ(A,rmax,f2tran,nr);
    }
    else tran_V1[it*N+jt]=tran_V1[it*N+jt]=tran_V3[it*N+jt]=0;
  }
 
  if(nanoshell_flag==1) //Nanoshells
  for(it=0;it<N;it++)
  for(jt=0;jt<N;jt++)
  if(li[it]==li[jt]+1||li[it]==li[jt]-1){
    tran_V1[it*N+jt]=integ(rmin, rin,f1tran,nr);
    tran_V2[it*N+jt]=integ(rin ,rout,f1tran,nr);
    tran_V3[it*N+jt]=integ(rin ,rout,f2tran,nr);
    tran_V4[it*N+jt]=integ(rout,rmax,f1tran,nr);
    tran_V5[it*N+jt]=integ(rout,rmax,f2tran,nr);
  }
}

// --- Transition rate from i to j ---------------------------------------

complex Mfi(int i, int j, int m, numero w)
{
  complex val=0,A1,C1,D1,G1,hh;
  numero ang=tran_ang(li[i],li[j],m);

  if(nanoshell_flag==0) val=-3./(eps(w)+2.)*tran_V1[i*N+j]-tran_V2[i*N+j]
                            +(eps(w)-1.)/(eps(w)+2.)*A*A*A*tran_V3[i*N+j];
 
  if(nanoshell_flag==1){
    hh=1/(1.-2.*(eps(w)-1)/(eps(w)+2)/(2*eps(w)+1)*rin*rin*rin/rout/rout/rout);
    A1=-9*eps(w)/(eps(w)+2)/(2*eps(w)+1)*hh;
    C1=-3/(eps(w)+2)*hh;
    D1=-3*(eps(w)-1)/(eps(w)+2)/(2*eps(w)+1)*rin*rin*rin*hh;
    G1=(eps(w)-1)/(eps(w)+2)*(rout*rout*rout-rin*rin*rin)*hh;
    val=A1*tran_V1[i*N+j]+C1*tran_V2[i*N+j]
       +D1*tran_V2[i*N+j]-   tran_V4[i*N+j]
       +G1*tran_V5[i*N+j];
  }  
 return ang*val; 
}

// --- Spectrum of hot electrons for certain frequency -------------------

numero e_spectrum(int f, numero w)  //Sum up the initial states
{
  numero Ge=0,wfi; int m,i,lmin;
 
  for(i=0; i<N; i++){
    lmin=(li[f]<li[i])?li[f]:li[i];
    wfi=Energy[f]-Energy[i];
    for(m=-lmin; m<=lmin; m++)
    Ge+=f_FERMI(Energy[i],KT)*(1.0-f_FERMI(Energy[f],KT))
        *mod2(Mfi(i,f,m,w))*1.0/tau/((wfi-w)*(wfi-w)+1./tau/tau);
  }
  return 2*pi*Ge; 
}

// --- Spectrum of hot holes for certain frequency -------------------

numero h_spectrum(int f, numero w)  //Sum up the initial states
{
  numero Gh=0,wfi; int m,i,lmin;
 
  for(i=0; i<N; i++){
    lmin=(li[f]<li[i])?li[f]:li[i];
    wfi=Energy[f]-Energy[i];
    for(m=-lmin; m<=lmin; m++) 
    Gh+=(1.-f_FERMI(Energy[i],KT))*f_FERMI(Energy[f],KT)
        *mod2(Mfi(i,f,m,w))*1.0/tau/((wfi+w)*(wfi+w)+1./tau/tau);
  }
  return 2*pi*Gh; 
}

// --- Number of hot electrons -------------------------------------------

numero num_e(numero w)
{
  numero Gt=0; int i;
  for(i=0;i<N;i++)
  Gt+=e_spectrum(i,w);
  return Gt;
}

// --- Carrier spatial distribution --------------------------------------


void init_carrier_distribution(numero w)
{
  numero wfi; int m,i,f,lmin,ii=0;
  int tot_states=0;                 //Total number of states 
  for(i=0; i<N; i++) tot_states+=(li[i]+1);
  spec_e= new numero[tot_states];
  spec_h= new numero[tot_states];
 
  for(f=0; f<N; f++)
  for(m=0; m<=li[f]; m++){
    spec_e[ii]=spec_h[ii]=0;
    for(i=0; i<N; i++)
    if(li[i]>=abs(m)){
      wfi=Energy[f]-Energy[i];
      spec_e[ii]+=2*pi*f_FERMI(Energy[i],KT)*(1.0-f_FERMI(Energy[f],KT))
                  *mod2(Mfi(i,f,m,w))*1.0/tau/((wfi-w)*(wfi-w)+1./tau/tau);
      spec_h[ii]+=2*pi*f_FERMI(Energy[f],KT)*(1.0-f_FERMI(Energy[i],KT))
                  *mod2(Mfi(i,f,m,w))*1.0/tau/((wfi+w)*(wfi+w)+1./tau/tau);
    }
    ii++;
  }
}

numero get_carrier_distribution(int opt, numero E0, numero E1,
                                         numero R, numero phi)
{
  numero val=0,cphi=cos(phi); R=R*A;
  long numero factor; 
  int i,m,ii=0,eflag;
  for(i=0;i<N;i++){
    eflag=1; if(Energy[i]<E1&&Energy[i]>E0) eflag=0;
    if(opt==0&&eflag==0)
     val+=spec_e[ii]*psiR(i,R)*psiR(i,R)*legendre(li[i],0,cphi)
          *legendre(li[i],0,cphi)*(2*li[i]+1);
    if(opt==1&&eflag==0)
     val+=spec_h[ii]*psiR(i,R)*psiR(i,R)*legendre(li[i],0,cphi)
          *legendre(li[i],0,cphi)*(2*li[i]+1);
    ii++;
    for(m=1; m<=li[i]; m++){
      factor=(2*li[i]+1)*fact2_n[li[i]*(li[i]+1)/2+m];
      if(opt==0&&eflag==0)
       val+=2*spec_e[ii]*psiR(i,R)*psiR(i,R)*legendre(li[i],abs(m),cphi)
            *legendre(li[i],abs(m),cphi)*factor;
      if(opt==1&&eflag==0)
       val+=2*spec_h[ii]*psiR(i,R)*psiR(i,R)
            *legendre(li[i],abs(m),cphi)*legendre(li[i],abs(m),cphi)*factor;
      ii++;
    }
  }
  return val/4/pi;
}

numero get_carrier_distribution(int opt, numero R, numero phi)
{
  return get_carrier_distribution(opt,-1000.0,1000.0,R,phi);
}


// --- input-file script -------------------------------------------------

numero *w0=NULL,*w1=NULL,*ww=NULL;  int nnw=0,*nwi=NULL;  // energy range(s)

void program(char *fin_name)
{
  numero w,R,phi,x,x0,x1,nx,sx,z,z0,z1,nz,sz,vale,valh,E0,E1;
  int     nwtot=0,w_flag,iw;
  int     i,j;
  
  char command[100], qntt[100], name[100];
  FILE *fin,*fout;
  fpos_t pos_calc;

// --- Set Openmp  //Command for doing parallel calculation

//n_core=omp_get_num_procs();
//omp_set_num_threads(n_core);

// --- Try to open input file.

  if((fin=fopen(fin_name,"r"))==NULL) {
    printf("*** error: cannot open file %s\n", fin_name);
    exit(1);
  }
  strcpy(tname,fin_name);    strcat(tname,".time"); // file name for time track

// --- Input-file interpreter.

  do {
    read_name(fin,command);

// --- definitions based upon = sign, with no blank spaces in between

    if((i=algebra_locate(command,'='))>0) {   // = sign at position i.
      command[i]=0;
      algebra_define(command,command+i+1);
    } else

// --- A

    if(!strcmpC(command,"A")) A=alread_numero(fin)*nm; else
    
// --- V0

    if(!strcmpC(command,"V0")) V0=alread_numero(fin)/au_eV; else
    
// --- Material

    if(!strcmpC(command,"Material")) Material=alread_int(fin); else

// --- EPS

    if(!strcmpC(command,"EPS")) EPS=alread_numero(fin); else

// --- de

    if(!strcmpC(command,"de")) de=alread_numero(fin)/au_eV; else

// --- nr

    if(!strcmpC(command,"nr")) nr=alread_int(fin); else

// --- rmax

    if(!strcmpC(command,"rmax")) rmax=alread_numero(fin)*A; else

// --- lmax

    if(!strcmpC(command,"lmax")) lmax=alread_int(fin); else

// --- T

    if(!strcmpC(command,"T")) KT=kb_au*alread_numero(fin); else

// --- tau

    if(!strcmpC(command,"tau")) tau=alread_numero(fin)*1e-12/au_s; else
    
// --- Nanoshell

    if(!strcmpC(command,"Nanoshell"))  nanoshell_flag=read_int(fin);  else

// --- Energies

    if(!strcmpC(command,"Energies")) {
      read_name(fin,name);
      if(!strcmpC(name,"write")) energies_flag=1;
      else                       energies_flag=2;
      read_name(fin,Energies_name);
    } else

// --- Wavefunctions

    if(!strcmpC(command,"Wavefunctions")) {
      read_name(fin,name);
      if(!strcmpC(name,"write")) wavefunctions_flag=1;
      else                       wavefunctions_flag=2;
      read_name(fin,Wavefunctions_name);
    } else

// --- grid_xz: grid for induced-charge rendering

    if(!strcmpC(command,"grid_xz")) {
      x0=alread_numero(fin);  x1=alread_numero(fin);
      nx=alread_int(fin); if(nx<1) nx=1;
      z0=alread_numero(fin);  z1=alread_numero(fin);
      nz=alread_int(fin); if(nz<1) nz=1;
      grid_flag=1; sx=(x1-x0)/(nx-1); sz=(z1-z0)/(nz-1);
    } else
    
// --- w

    if(!strcmpC(command,"w")) {
      if(nnw>0) {delete [] w0;  delete [] w1;  delete [] ww;  delete [] nwi;}
      nnw=alread_int(fin);
      w0=new numero[nnw];  w1=new numero[nnw];  nwi=new int[nnw];  nwtot=0;
      for(i=0; i<nnw; i++) {
        w0[i]=alread_numero(fin)/au_eV;
        w1[i]=alread_numero(fin)/au_eV;
        nwi[i]=alread_int(fin);  nwtot+=nwi[i];
      }
      ww=new numero[nwtot];
      for(i=0,j=0; i<nnw; i++)  for(iw=0; iw<nwi[i]; iw++,j++)
        if(iw==0) ww[j]=w0[i]; else ww[j]=w0[i]+iw*(w1[i]-w0[i])/(nwi[i]-1);
    } else

// --- begin-calculation ... end-calculation loop and calc-... commands

    if(!strcmpC(command,"begin-calculation")) {
      fgetpos(fin,&pos_calc);      // file pointer pointing right after "begin"
      if(nanoshell_flag==1&&wavefunctions_flag!=2)
       on_error("calculation",
                "Nanoshell calculations require reading wavefunctions.");

      #define openfile(W) read_name(fin,name);          \
               if(!strcmpC(name,"screen")) fout=stdout; \
               else  if(W) fout=fopen(name,"w"); else fout=fopen(name,"a");
      #define closefile if(strcmpC(name,"screen")) fclose(fout);

      init_fact_factor(lmax); //Initializes the factors involving factorials
                              //needed in the calculation of bessel functions
                              //and in the Ylm
      get_wavefunctions();
     
      hot_ele_time("Initialization of transition elements:    ");
      set_tran_r();
      hot_ele_time("                                          ");
     
      for(iw=0; iw<nwtot; iw++) {                      // w loop
        w=ww[iw];
                  
        fsetpos(fin,&pos_calc); // set file pointer right after "begin"
        do{
          read_name(fin,qntt);
          
          // --- e-spectrum
          
          if(!strcmpC(qntt,"calc-e-spectrum")) {
            hot_ele_time("Calculation of e-spectrum:                ");
            openfile(iw==0)
            for(i=0; i<N; i++) fprintf(fout,"%g %g %g\n",w*au_eV
                                                        ,Energy[i]*au_eV
                                                        ,e_spectrum(i,w));
            closefile
            hot_ele_time("                                          ");
          }

          // --- h-spectrum

          if(!strcmpC(qntt,"calc-h-spectrum")) {
            hot_ele_time("Calculation of h-spectrum =               ");
            openfile(iw==0)
            for(i=0; i<N; i++) fprintf(fout,"%g %g %g\n",w*au_eV
                                                        ,Energy[i]*au_eV
                                                        ,h_spectrum(i,w));
            closefile
            hot_ele_time("                                          ");
          }

          // --- num-spectrum
          
          if(!strcmpC(qntt,"calc-num-spectrum")) {
            hot_ele_time("Calculation of num-spectrum:              ");
            openfile(iw==0)
  
            fprintf(fout,"%g %g\n",w*au_eV,num_e(w));
            closefile
            hot_ele_time("                                          ");
          }
          
          // --- carrier-distribution
          
          if(!strcmpC(qntt,"calc-carrier-distribution")) {
            hot_ele_time("Calculation of carrier distribution:      ");
            if(grid_flag<1) on_error("carrier-distribution",
                             "a grid must be defined for carrier-distribution");
            openfile(iw==0)
            init_carrier_distribution(w);
            for(int ix=0; ix<nx; ix++){
              if(ix==0) x=x0; else x=x0+ix*sx;
              for(int iy=0; iy<nz; iy++){
                if(iy==0) z=z0; else z=z0+iy*sz;
                R=sqrt(x*x+z*z); phi=atan2(x,z);
                vale=get_carrier_distribution(0,R,phi);
                valh=get_carrier_distribution(1,R,phi);
                fprintf(fout,"%g %g %g %g %g\n",w*au_eV,x,z,vale,valh);
              }
          //  fprintf(fout,"\n");  //Liu
            }

            closefile   
            hot_ele_time("                                          ");
          }
          
          if(!strcmpC(qntt,"calc-carrier-distribution-E")) {
            hot_ele_time("Calculation of carrier distribution:      ");
            E0=alread_numero(fin)/au_eV; E1=alread_numero(fin)/au_eV;
            if(grid_flag<1) on_error("carrier-distribution",
                             "a grid must be defined for carrier-distribution");
            openfile(iw==0)
            init_carrier_distribution(w);
            for(int ix=0; ix<nx; ix++){
              if(ix==0) x=x0; else x=x0+ix*sx;
              for(int iy=0; iy<nz; iy++){
                if(iy==0) z=z0; else z=z0+iy*sz;
                R=sqrt(x*x+z*z); phi=atan2(x,z);
                vale=get_carrier_distribution(0,E0,E1,R,phi);
                valh=get_carrier_distribution(1,E0,E1,R,phi);
                fprintf(fout,"%g %g %g %g %g\n",w*au_eV,x,z,vale,valh);
              }
          //  fprintf(fout,"\n");  //Liu
            }

            closefile   
            hot_ele_time("                                          ");
          }


        }while(strcmpC(qntt,"end-calculation"));
        
    } } else

// --- skipline

    if(!strcmpC(command,"skipline"))  printf("\n"); else

// --- time

    if(!strcmpC(command,"time"))  tflag=read_int(fin);  else

// --- end, exit: terminate the execution of the program

    if(strcmpC(command,"end") && strcmpC(command,"exit"))
      on_error("HOT_ELECTRON", "command not found:", command);

  } while(strcmpC(command,"end") && strcmpC(command,"exit"));
  fclose(fin);
  
  delete [] aa;      delete [] bb;
  delete [] occ;     delete [] Energy; 
  delete [] li;      delete [] ni;      
  delete [] tran_V1; delete [] tran_V2; 
  delete [] tran_V3; delete [] tran_V4;
  delete [] tran_V5;  
  delete [] pos;    
  delete [] State;   
  delete [] spec_e;  delete [] spec_h;   
  delete [] fact2_n; delete [] fact1_n;
}

int main(int argc, char **argv)
{
// --- Help on the use of the program if no parameters are given.

  if(argc<2) {
    printf("This program requires an input file.\n");
    printf("The correct call syntax is as follows:\n\n");
    printf("        hot_elec input\n\n");
    printf("Use the commands 'HOT_ELE help' or 'HOT_ELE h'\n");
    printf("for a detailed description on the input file format.\n");
    return 0;
  }

// --- Read first argument and print out help if that is what is wanted.

  if(!strcmp(argv[1],"h") || !strcmp(argv[1],"help"))
    {printf(HOT_ELE_HELP);  return 0;}
  int i;
  for(i=1; i<argc; i++){
    program(argv[i]);
  }
  return 0;
}
