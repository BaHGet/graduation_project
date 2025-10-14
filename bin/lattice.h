#ifndef jga_lattice       // ************************************************
#define jga_lattice 1     // ***  jga/lattice.h                  11-v-01  ***
                          // ***                                          ***
#include "jga.h"          // ***                                          ***
                          // ***                                          ***
// **************************************************************************
// ***                                                                    ***
// ***  Author: F. J. Garcia de Abajo, CSIC, Spain                        ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
//    Creation of a list of n_g reciprocal lattice vectors selected among
//    the closest ones with respect to the origin. (i=0 -> Gx=Gy=0 and Qg=Q.)
//
//    (Qgx[i],Qgy[i])=Q+g  for i=0, ..., n_g  after calling lattice(gmax).
//    n_g is approx. gmax^2 (not quite, since all symm.-equiv. g's are incl.)
//
//    lattice(gmax)  -> gmax is the radius of the circle in which reciprocal
//    lattice vectors are included, given in units of the equivalent radius
//    g of the circle with the same area as the reciprocal lattice unit cell
//    -----------------------------------------------------------------------

numero  a0=0,b0=0,alpha,beta;     // lattice parameters
numero  Qx=0,Qy=0;                // parallel momentum

numero  area;                     // surface unit cell area
numero  Gax,Gay,Gbx,Gby;          // reciprocal-lattice unit vectors
int     n_g=0;                    // number of beams
numero  *Qgx,*Qgy;                // Q+g beam vectors
numero  *RQg,*phiQg;              // Q+g in polar coordinates
int     flag_polar_Qg=0;          // 1 if modQg and phiQg have been used

void lattice_free(void)
{
  if(n_g>0) {
    delete [] Qgx;  delete [] Qgy;  n_g=0;
    if(flag_polar_Qg) {delete [] RQg;  delete [] phiQg;}  flag_polar_Qg=0;
} }

void lattice(numero gmax, int polar_Qg)
{
  lattice_free();
  if(ABS(a0+b0)<1e-20)  on_error("lattice", "undefined lattice");

  // --- initialize lattice
                                                    // Cartesian coordinates
  numero  ax=a0*cos(alpha), bx=b0*cos(alpha+beta);  // of unit vectors a, b
  numero  ay=a0*sin(alpha), by=b0*sin(alpha+beta);  // of 2D real-space latt.

  numero  axb=ax*by-ay*bx;                          // a x b

  Gax= 2*pi*by/axb;  Gay=-2*pi*bx/axb;              // 2D reciprocal vectors
  Gbx=-2*pi*ay/axb;  Gby= 2*pi*ax/axb;              // with Ga*b=Gb*a=0 and
                                                    // Ga*a=Gb*b=2*pi;
  area=ABS(axb);                                    // surface unit cell area

  // --- reciprocal lattice vectors

  numero GG=2*pi*gmax/sqrt(pi*area), GG2;  GG2=GG*GG;

  int na=int(GG/sqrt(Gax*Gax+Gay*Gay)+1);
  int nb=int(GG/sqrt(Gbx*Gbx+Gby*Gby)+1);
  int dim=(2*na+1)*(2*nb+1);

  numero *gx;  gx=new numero [dim];                 // lattice vector coord.
  numero *gy;  gy=new numero [dim];
  numero *dg;  dg=new numero [dim];
  int   *occ;  occ=new int   [dim];

  int i,j,ij=0;  n_g=0;                             // number of beams;

  for(i=-na; i<=na; i++)
  for(j=-nb; j<=nb; j++, ij++)  {
    gx[ij]=i*Gax+j*Gbx;
    gy[ij]=i*Gay+j*Gby;
    dg[ij]=sqr(gx[ij])+sqr(gy[ij]);
    if(dg[ij]<=GG2)  {occ[ij]=1;  n_g++;}  else  occ[ij]=0;
  }
 
  // --- beam vectors

  Qgx=new numero [n_g];  Qgy=new numero [n_g];
  if(polar_Qg) {
    RQg=new numero [n_g];  phiQg=new numero [n_g];  flag_polar_Qg=1;
  }

  numero val;
  int    n;

  for(n=0; n<n_g; n++) {
    val=infinite;
    for(i=0; i<dim; i++) if(occ[i]) if(dg[i]<val) {val=dg[i]; j=i;}
    Qgx[n]=Qx+gx[j];  Qgy[n]=Qy+gy[j];  occ[j]=0;
    if(polar_Qg) {
      RQg[n]=sqrt(sqr(Qgx[n])+sqr(Qgy[n]));  phiQg[n]=atan2(Qgy[n],Qgx[n]);
  } }

  delete [] gx;  delete [] gy;  delete [] dg;  delete [] occ;
}

void lattice(numero gmax) {lattice(gmax,0);}

void triangular_lattice(numero gmax, numero a)
  {a0=b0=a;  alpha=-pi/6;  beta=pi/3;  lattice(gmax,1);}

void square_lattice(numero gmax, numero a)
  {a0=b0=a;  alpha=0;  beta=pi/2;  lattice(gmax,1);}


int excursion_1BZ(int i, int nQ, numero &q)
{
  numero r, Mx,My,Kx,Ky,Xx,Xy, fiq=atan2(Gay,Gax), G=sqrt(sqr(Gax)+sqr(Gay));

  if(ABS(a0/b0-1)>1e-5) return 0;

  if(ABS(ABS(beta)-pi/3)<=1e-5) { // triangular lattice
    Mx=Gax/2;  My=Gay/2;
    Kx=G/sqrt(3)*cos(fiq+pi/6);
    Ky=G/sqrt(3)*sin(fiq+pi/6);

    if(nQ<=1) q=0;                     //    G    M        K            G
    else      q=i*(3+sqrt(3))/(nQ-1);  // q: 0 -> 3^0.5 -> 1+3^0.5+1 -> 3+3^0.5
    if(q>sqrt(3)+1) {r=1-(q-sqrt(3)-1)/2;  Qx=r*Kx;  Qy=r*Ky;}          else
    if(q>sqrt(3))   {r=q-sqrt(3);  Qx=Mx+(Kx-Mx)*r;  Qy=My+(Ky-My)*r;}  else
                    {r=q/sqrt(3);  Qx=r*Mx;          Qy=r*My;}
    return 1;
  }

  if(ABS(ABS(beta)-pi/2)<=1e-5) { // square lattice,  s2=2^0.5
    Xx=Gax/2;  Xy=Gay/2;
    Mx=G/sqrt(2)*cos(fiq+pi/4);
    My=G/sqrt(2)*sin(fiq+pi/4);

    if(nQ<=1) q=1;                     //    G    X    M      G
    else      q=i*(2+sqrt(2))/(nQ-1);  // q: 0 -> 1 -> 2 -> 2+2^0.5
    if(q>2) {r=1-(q-2)/sqrt(2);  Qx=r*Mx;          Qy=r*My;}          else
    if(q>1) {r=q-1;              Qx=Xx+(Mx-Xx)*r;  Qy=Xy+(My-Xy)*r;}  else
            {r=q;                Qx=r*Xx;          Qy=r*Xy;}
    return 1;
  }

  return 0;
}

#endif  // ******************************************************************
