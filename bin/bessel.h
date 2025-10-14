#ifndef jga_bessel        // ************************************************
#define jga_bessel 1      // ***  jga/bessel.h                            ***
                          // ***                                          ***
#include "jga.h"          // ***  BESSEL FUNCTIONS Kn(z), Jn(x), Yn(x)    ***
#include "complex.h"      // ***                                          ***
                          // ***                                          ***
// --------------------------------------------------------------------------
// The function besselK(n,z) returns the Bessel function K_n(z), for
// integer n and complex z. (Author: F. J. Garcia de Abajo, CSIC, Spain.)
//
// The routine besselJY(n,x,J,Y,Jp,Yp) accepts real parameters n and x
// on input and returns the Bessel functions J_n(x)=J, J'_n(x)=Jp,
// Y_n(x)=Y, and Y'_n(x)=Yp. (Adapted from Numerical Recipies.)
// --------------------------------------------------------------------------


// --------------------------------------------------------------------------
// J_n(x), Y_n(x), J'_n(x), Y'_n(x),  real n,  real x

double chebev(double a, double b, double c[], int m, double x)
{
  double d=0,dd=0,sv,y,y2;  int j;

  if((x-a)*(x-b)>0)  on_error("chebev", "x not in range");
  y=(2*x-a-b)/(b-a);  y2=2*y;
  for(j=m; j>=1; j--) {sv=d;  d=y2*d-dd+c[j];  dd=sv;}

  return y*d-dd+c[0]/2;
}

// --------------------------------------------------------------------------

void beschb (double x, double &g1, double &g2, double &gamp, double &gamm)
{
  double xx,c1[7],c2[8];

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

int besselJY(double xnu, double x, numero &rj,  numero &ry,
             numero &rjp, numero &ryp)
{
  int i,isign=1,l,nl, maxit=100000, salir;
  double xmin=2.0, eps=1.0e-10, fpmin=1.0e-30,
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
  if(salir)  on_error("besselJY", "x too large");

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

    if(salir)  on_error("besselJY", "series failed to converge");

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
    if(salir)  on_error("besselJY", "cf2 failed");
    gam=(p-f)/q;
    rjmu=sqrt(w/((p-f)*gam+q));
    if(rjl<0) rjmu=-rjmu;
    rymu=rjmu*gam;
    rymup=rymu*(p+q/gam);
    ry1=xmu*xi*rymu-rymup;
  }

  fct=rjmu/rjl;
  rj=rjl1*fct;
  rjp=rjp1*fct;
  for(i=1; i<=nl; i++) {
    rytemp=(xmu+i)*xi2*ry1-rymu;
    rymu=ry1;  ry1=rytemp;
  }

  ry=rymu;  ryp=xnu*xi*rymu-ry1;

  return 0;
}

// --------------------------------------------------------------------------
// I_n(x), K_n(x), I'_n(x), K'_n(x),  int n,  complex x

complex besselI0(complex x)
{
  complex  ax,y;

  if(mod(x)<3.75) {
    y=x/3.75;  y=y*y;
    return  1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*
                  (0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  } else {
    ax=(real(x)<0)?-x:x;  y=3.75/ax;                            // jga:
    return  (exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1      // this part
              +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2    // is wrong
              +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
              +y*0.392377e-2))))))));
  }
}

// --------------------------------------------------------------------------

complex besselI1(complex x)
{
  complex  ax,y,ans;

  if(mod(x)<3.75) {
    y=x/3.75;  y=y*y;
    return  x*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
               +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  } else {
    ax=(real(x)<0)?-x:x;  y=3.75/ax;         //jga: this part is wrong
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
                  +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    return  (exp(ax)/sqrt(ax))*ans;
  }
}

// --------------------------------------------------------------------------

complex besselK0(complex x)
{
  complex  y;

  if(mod(x)<=2) {
    y=x*x/4;
    return  (-log(x/2.0)*besselI0(x))+(-0.57721566+y*(0.42278420
               +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
               +y*(0.10750e-3+y*0.74e-5))))));
  } else {
    y=2/x;
    return  (exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
               +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
               +y*(-0.251540e-2+y*0.53208e-3))))));
} }

// --------------------------------------------------------------------------

complex besselK1(complex x)
{
  complex  y;

  if(mod(x)<=2) {
      y=x*x/4;
      return  (log(x/2.0)*besselI1(x))+(1.0/x)*(1.0+y*(0.15443144
                  +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
                  +y*(-0.110404e-2+y*(-0.4686e-4)))))));
  } else {
      y=2/x;
      return  (exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
                +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
                 +y*(0.325614e-2+y*(-0.68245e-3)))))));
} }

// --------------------------------------------------------------------------

complex besselK(int n, complex x)
{
  complex  tox,bkp,bkm,bk;
  int  j;

  if(n<0)   n=-n;
  if(n==0)  return  besselK0(x);
  if(n==1)  return  besselK1(x);

  tox=2/x;
  bkm=besselK0(x);
  bk=besselK1(x);

  for(j=1; j<n; j++) {
    bkp=bkm+j*tox*bk;    bkm=bk;    bk=bkp;
  }

  return  bk;
}

#endif  // ******************************************************************
