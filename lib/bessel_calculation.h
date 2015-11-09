#include "jga.h"    
#include "complex.h"

//
//Library to calculate the spherical bessel function j_l(x) and the 
//spherical hankel function h_l^(1)(ix) as well as their first derivatives
//
//It is essentially a modification of bessel.h from jga library that allows 
//to compute these functions for values of l as large as 400 and for arguments
//from 0 up to ..
//the original routines do not work for such high values of l
//
//The exact routines taken from bessel.h are: chebev,beschb,besselJY,besseljy
//These routines have been modified changing numero by long numero and including
//a factor 1e300 to increase the value 
//
//We have added also a routine to initialize and calculate diferent factors involving
//factorials that are needed for the calculation of h_l^(1)(ix) and in the main code
//The routine is inspired by the similar one in jga.h but with a larger value of l_max 
//and also changing numero by long numero
//
//IMPORTANT 
//init_fact_factor must be called before using besselh1i, besselh1ip or besselh1ip_h
//
//Final routines:
//
//  besselj_n(int l,numero x)     returns j_l(x)
//
//  besseljp_n(int l, numero x)   returns j'_l(x)
//
//  besseljp_j(int l, numero x)   returns j'_l(x)/j_l(x)
//
//  besselh1i(int l, numero x)    returns h^(1)_l(ix) x real
//
//  besselh1ip(int l, numero x)   returns h'^(1)_l(ix) x real
//
//  besselh1ip_h(int l, numero x) returns h'^(1)_l(ix)/h^(1)_l(ix) x real
//
//


//************************************************************************
//--- Routines to initialize the expressions involving factorials -------*
//************************************************************************

//int l_max=801;               // maxium of l
long numero *fact1_n;        // fact1_n[l*(l+1)/2+m]=(l+m)!*/m!/(l-m)!
long numero *fact2_n;        // fact2_n[l*(l+1)/2+m]=(l-m)!/(l+m)!
int fact1_initialized_n=0;    // flag to show if fact has been init.

int init_fact_factor(int l_max)       // init fact1_n
{
  if(fact1_initialized_n)  return 0;
  int l,m,n_l;
  fact1_n = new long numero [(l_max+2)*(l_max+1)/2];  fact1_initialized_n=1;
  fact2_n = new long numero [(l_max+2)*(l_max+1)/2];
  for(l=0;l<=l_max;l++)
    for(m=0;m<=l;m++){
      n_l=l*(l+1)/2;
      if(m==0) {
        fact1_n[n_l+m]=1.;
	fact2_n[n_l+m]=1.;
      }
      else{
        fact1_n[n_l+m]=fact1_n[n_l+m-1]*(l+m)*(l-m+1)/m;
	fact2_n[n_l+m]=fact2_n[n_l+m-1]/(l+m)/(l-m+1);
      }
    }
  return 0;
}


//************************************************************************
//--- Routines taken from bessel.h  -------------------------------------*
//************************************************************************

long numero chebev(numero a, numero b, long numero c[], int m, numero x)
{
  long numero d=0,dd=0,sv,y,y2;  int j;

  if((x-a)*(x-b)>0)  on_warning("chebev", "x not in range");
  y=(2*x-a-b)/(b-a);  y2=2*y;
  for(j=m; j>=1; j--) {sv=d;  d=y2*d-dd+c[j];  dd=sv;}

  return y*d-dd+c[0]/2;
}

// --------------------------------------------------------------------------

void beschb (numero x, long numero &g1, long numero &g2, long numero &gamp, long numero &gamm)
{
  long numero xx,c1[7],c2[8];

  c1[0] = -1.142022680371172;  c1[1] = 0.006516511267076;
  c1[2] = 0.000308709017308;   c1[3] = -3.470626964e-6;
  c1[4] = 6.943764e-9;         c1[5] = 3.678e-11;
  c1[6] = -1.36e-13;
  c2[0] = 1.843740587300906;   c2[1] = -0.076852840844786;
  c2[2] = 0.001271927136655;   c2[3] = -4.971736704e-6;
  c2[4] = -3.3126120e-8;       c2[5] = 2.42310e-10;
  c2[6] = -1.7e-13;            c2[7] = -1.0e-15;

  xx=8*x*x-1;
  g1=chebev(-1,1,c1,6,xx);     g2=chebev(-1,1,c2,7,xx);
  gamp=g2-x*g1;                gamm=g2+x*g1;
}

// --------------------------------------------------------------------------

int besselJY(numero xnu, numero x, long numero &rj,  long numero &ry,
             long numero &rjp, long numero &ryp)
{
  int i,isign=1,l,nl, maxit=100000, salir;
  long numero xmin=2.0, eps=1.0e-10, fpmin=1.0e-30,
	 a,b,br,bi,c,cr,ci,d=0,del,del1,den,di,dlr,dli,
       	 dr,e,f,fct,fct2,fct3,ff,gam,gam1,gam2,gammi,gampl,h,
      	 p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,
       	 rymu,rymup,rytemp,sum,sum1,temp,w,x2,xi,xi2,xmu,xmu2;

  if(x<=0)  on_error("besselJY", "wrong argument x =", x, "");
  if(xnu<0)  on_error("besselJY", "wrong order n =", xnu, "");
  if(x<xmin)    nl=int(xnu+0.5);
  else {nl=int(xnu-x+1.5); nl=0>nl?0:nl;}
  xmu=xnu-nl;  xmu2=xmu*xmu;  xi=1/x;  xi2=2*xi;  w=xi2/pi;  h=xnu*xi;

  if(h<fpmin) h=fpmin;     c=h;  b=xi2*xnu;

  for(i=1, salir=1; i<=maxit && salir; i++) {
    b+=xi2;  d=b-d;             if(fabs(d)<fpmin) d=fpmin;
    c=b-1/c;                    if(fabs(c)<fpmin) c=fpmin;
    d=1/d;  del=c*d;  h=del*h;  if(d<0) isign=-isign;
    if(fabs(del-1)<eps) salir=0; else salir=1;
  }
  if(salir) on_warning("besselJY", "x too large", x, "");

  rjl1=rjl=isign*fpmin;   rjp1=rjpl=h*rjl;  fct=xnu*xi;

  for(l=nl; l; l--) {rjtemp=fct*rjl+rjpl;
                     fct-=xi;  rjpl=fct*rjtemp-rjl;  rjl=rjtemp;}

  if(rjl==0.0) rjl=eps;
  f=rjpl/rjl;

  if(x<xmin) {
    x2=x/2;  pimu=pi*xmu;
    if(fabs(pimu)<eps)  fct=1;  else fct=pimu/sin(pimu);
    d=-log(x2);  e=xmu*d;
    if(fabs(e)<eps)     fct2=1; else fct2=sinh(e)/e;

    beschb(xmu,gam1,gam2,gampl,gammi);
    ff=2/pi*fct*(gam1*cosh(e)+gam2*fct2*d);
    e=exp(e);   p=e/(gampl*pi);
    q=1/(e*pi*gammi);
    pimu2=pimu/2;
    if(fabs(pimu2)<eps) fct3=1; else fct3=sin(pimu2)/pimu2;
    r=pi*pimu2*fct3*fct3;
    c=1;  d=-x2*x2;  sum=ff+r*q;  sum1=p;

    for(i=1, salir=1; i<=maxit && salir; i++) {
      ff=(i*ff+p+q)/(i*i-xmu2);   c*=d/i;
      p/=(i-xmu);   q/=(i+xmu);   del=c*(ff+r*q);
      sum+=del;   del1=c*p-i*del;
      sum1+=del1;
      if(fabs(del)<(1+fabs(sum))*eps) salir=0; else salir=1;
    }

    if(salir)  on_warning("besselJY", "series failed to converge");

    rymu=-sum;  ry1=-sum1*xi2;  rymup=xmu*xi*rymu-ry1;
    rjmu=w/(rymup-f*rymu);

  } else {

    a=0.25-xmu2;  p=-xi/2;  q=1;  br=2*x;  bi=2;
    fct=a*xi/(p*p+q*q);  cr=br+q*fct;  ci=bi+p*fct;
    den=br*br+bi*bi;  dr=br/den;  di=-bi/den;  dlr=cr*dr-ci*di;
    dli=cr*di+ci*dr;  temp=p*dlr-q*dli;  q=p*dli+q*dlr;
    p=temp;
    for(i=2, salir=1; i<=maxit && salir; i++) {
      a+=2*(i-1);  bi+=2;  dr=a*dr+br;  di=a*di+bi;
      if(fabs(dr)+fabs(di)<fpmin) dr=fpmin;
      fct=a/(cr*cr+ci*ci);  cr=br+cr*fct;  ci=bi-ci*fct;
      if(fabs(cr)+fabs(ci)<fpmin) cr=fpmin;
      den=dr*dr+di*di;  dr=dr/den;  di=-di/den;
      dlr=cr*dr-ci*di;  dli=cr*di+ci*dr;  temp=p*dlr-q*dli;
      q=p*dli+q*dlr;  p=temp;
      if(fabs(dlr-1)+fabs(dli)<eps) salir=0; else salir=1;
    }
    if(salir)  on_warning("besselJY", "cf2 failed");
    gam=(p-f)/q;
    rjmu=sqrt(w/((p-f)*gam+q));
    if(rjl<0) rjmu=-rjmu;
    rymu=rjmu*gam;
    rymup=rymu*(p+q/gam);
    ry1=xmu*xi*rymu-rymup;
  }
  
  fct=rjmu/rjl;  //Liu
  rj=rjl1*fct;
  rjp=rjp1*fct;
  for(i=1; i<=nl; i++) {
    rytemp=(xmu+i)*xi2*ry1-rymu;
    rymu=ry1;  ry1=rytemp;
  }

  ry=rymu;  ryp=xnu*xi*rymu-ry1;

  return 0;
}

void besseljy(int l, numero x, long numero &j, long numero &y)
{
  long numero jp,yp;
  long numero factor=sqrt(pi/(2*x));

  besselJY(l+0.5,x,j,y,jp,yp);

  j*=factor;  y*=factor;
}

//************************************************************************
// --- Spherical Bessel functions ---------------------------------------*
//************************************************************************

long numero init_besselj_n(int l, numero x)   //Besselj for jp/j
{
  if(ABS(x)<1e-10)  if(l==0) return 1.0;     else
                    if(l==1) return x/3;     else
                    if(l==2) return x*x/15;  else  return 0;
  long numero jl,yl;
  besseljy(l,x,jl,yl);
  return jl;
}

numero besselj_n(int l, numero x)  
{
  long numero j;
  j=init_besselj_n(l, x);
  //if(j<1e-300) return 0;     //Liu
  return j;
}

numero besseljp_n(int l, numero x)
{
  if(l==0) return -besselj_n(1,x);
  else return besselj_n(l-1,x)-((l+1)/x)*besselj_n(l,x);
}

numero besseljp_j(int l, numero x)  //jp/j
{
  long numero jp,j;
  j=init_besselj_n(l,x);
  if(l==0) jp=-init_besselj_n(1,x);
  else jp=init_besselj_n(l-1,x)-((l+1)/x)*j;
  return jp/j;
}

//************************************************************************
// --- Spherical Hankel functions ---------------------------------------*
//************************************************************************

complex besselh1i(int l, numero x) //besselh1(l,ix)	int l, real x
{
  complex factor=i_l(l);
  int m,s;
  long numero sum=1,p;
  for(m=1; m<=l; m++){
    p=1.; for(s=1; s<=m; s++) p=p*2*x;
    sum+=fact1_n[l*(l+1)/2+m]/p;
  }
  l=l%2; if(l==0) factor=-1.*factor;
  return sum*factor*exp(-x)/x;
}

complex besselh1ip(int l, numero x) //besselh1'(l,ix)	int l, real x
{
  if(l==0) return -(1/x/x+1/x)*exp(-x)*i_c;
  else  return -0.5*besselh1i(l+1,x)+0.5*besselh1i(l-1,x)-0.5*besselh1i(l,x)/x/i_c;
}

complex besselh1ip_h(int l, numero x) //h1'(ix)/h1(ix) 
{
  int m,s;
  long numero sum1,sum2,p;
  sum1=0;
  sum2=0;
  for(m=0;m<=l;m++){
    p=fact1_n[l*(l+1)/2+m];
    if(x>0.5)
      for(s=1;s<=m;s++) p=p/x/2;
    else
      for(s=1;s<=l-m;s++) p=p*2*x;
    sum1+=(m+1)/x*p;
    sum2+=p;
  }
  return i_c*(1+sum1/sum2);
}
