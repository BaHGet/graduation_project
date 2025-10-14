#ifndef jga_matrix        // ************************************************
#define jga_matrix 1      // ***  jga/matrix.h                            ***
                          // ***                                          ***
#include "jga.h"          // ***                                          ***
#include "complex.h"      // ***                                          ***
                          // ***                                          ***
// ***  ----------------  // ***  --------------------------------------  ***
// ***                                                                    ***
// ***  Author: F. J. Garcia de Abajo, CSIC, Spain                        ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Matrix algebra (complex matrices)                                 ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Operations:                                                       ***
// ***                                                                    ***
// ***    For a nxm matrix,  0<=i<n  and  0<=j<m                          ***
// ***    A*x  x*A  A/x  -A              A(nxm)                           ***
// ***    A+B  A-B  A=B                  A(nxm) and B(nxm)                ***
// ***    x+A  A+x  A-x  x-A  A=x  x/A   A(nxn)                           ***
// ***    A*B                            A(nxm) and B(mxq)                ***
// ***    A/B                        <-  A(nxm) and B(mxm)                ***
// ***                                                                    ***
// ***    --A                        <-  A(nxn) inverse                   ***
// ***                                                                    ***
// ***    A.a(i,j,x)                 <-  A(i,j) = x                       ***
// ***    A.a(i,j)                   <-  gives A(i,j)                     ***
// ***    A.a(i,x)                   <-  like A.a(i,0,x)                  ***
// ***    A.a(i)                     <-  like A.a(i,0)                    ***
// ***                                                                    ***
// ***    A.add(i,j,x)               <-  add x to A.a(i,j)                ***
// ***    A.add(i,x)                 <-  add x to A.a(i,0)                ***
// ***    A.alloc(n,m)               <-  memory allocation (nxm)          ***
// ***    A.alloc(n)                 <-  like A.alloc(n,1)                ***
// ***                                                                    ***
// ***    A.free()                   <-  memory de-allocation             ***
// ***    r*A  A*r                   <-  A(nxm) and r (el. by el.)        ***
// ***                                                                    ***
// **************************************************************************


// **************************************************************************
// *** declaration                                                        ***
// **************************************************************************

class matrix {
public:

  int n,m,nm, sel;                        // sel=0  -> general
  complex *p;                             // sel=1  -> diagonal

  matrix(void) {n=m=nm=sel=0;}
  void free(void)  {if(nm) {delete [] p;}  n=m=nm=0;}
  matrix(int i, int j)  {nm=0;  alloc(i,j);}
  matrix(int i) {nm=0;  alloc(i,1);}
  matrix(const matrix &A);
  ~matrix(void) {free();}
  void alloc(int i, int j);
  void alloc(int i)  {alloc(i,1);}
  void alloc_diagonal(int i) {alloc(i,1);  sel=1;}
  void a(int i, int j, complex x) {p[i*m+j]=x;}
  void a(int i, int j, numero x) {p[i*m+j]=x;}
  void a(int i, complex x) {p[i]=x;}
  void a(int i, numero x) {p[i]=x;}
  void add(int i, int j, complex x) {p[i*m+j]+=x;}
  void add(int i, int j, numero x) {p[i*m+j]+=x;}
  void add(int i, complex x) {p[i]+=x;}
  void add(int i, numero x) {p[i]+=x;}
  complex a(int i, int j);
  complex a(int i) {if(i>-1 && i<nm)  return p[i];  else  return complex(0,0);}
  matrix &operator=(complex x);
  matrix &operator=(numero x);
  matrix &operator=(const matrix &A);
  matrix &operator--(void);                // matrix inversion
  complex gaussj(const matrix &Mb, int n, int m, int opt);
  complex gaussj(const matrix &Mb, int n, int m)   // return ln(|Mb|)
    {return gaussj(Mb,n,m,0);}
//------------------------------------------------------------
  friend matrix operator+(complex x, matrix A);
  friend matrix operator+(numero x, matrix A) {return complex(x,0)+A;};
  friend matrix operator+(matrix A, complex x) {return x+A;}
  friend matrix operator+(matrix A, numero x) {return x+A;}
  friend matrix operator+(matrix A, const matrix &B);
  friend matrix operator-(matrix A);
  friend matrix operator-(complex x, matrix A);
  friend matrix operator-(numero x, matrix A) {return complex(x,0)-A;};
  friend matrix operator-(matrix A, complex x) {return (-x)+A;}
  friend matrix operator-(matrix A, numero x) {return (-x)+A;}
  friend matrix operator-(matrix A, const matrix &B);
  friend matrix operator*(complex x, matrix A);
  friend matrix operator*(numero x, matrix A);
  friend matrix operator*(matrix A, complex x) {return x*A;}
  friend matrix operator*(matrix A, numero x) {return x*A;}
  friend matrix operator*(const matrix &A, const matrix &B);
  friend matrix operator*(numero *r, matrix A);
  friend matrix operator*(const matrix &A, numero *r) {return r*A;}
  friend matrix operator/(complex x, matrix A) {return x*(--A);}
  friend matrix operator/(numero x, matrix A) {return x*(--A);}
  friend matrix operator/(matrix A, complex x) {return (1/x)*A;}
  friend matrix operator/(matrix A, numero x) {return (1/x)*A;}
  friend matrix operator/(matrix A, matrix B) {return A*(--B);}
};


// **************************************************************************
// *** implementation                                                     ***
// **************************************************************************

matrix::matrix(const matrix &A)
{
  int ij;
  nm=0;  if(A.nm) alloc(A.n,A.m);  sel=A.sel;

  for(ij=0; ij<nm; ij++)  p[ij]=A.p[ij];
}

complex matrix::a(int i, int j)
{
  if(i>-1 && i<n) {
    if(sel==1)  if(i==j)          return p[i];
    if(sel==0)  if (j>-1 && j<m)  return p[i*m+j];
  }
  return 0;
}

void matrix::alloc(int i, int j)
{
  if(nm)  delete [] p;   n=i; m=j; nm=i*j;  sel=0;
  if(nm) {p=new complex [nm];  for(i=0; i<nm; i++) p[i]=0;}
}

matrix operator+(complex x, matrix A)
{
  int i,j,ij;

  if(A.sel==0) {
    if(A.n!=A.m)  on_error("matrix","A is not square in x+A");
    for(i=ij=0; i<A.n; i++)  for(j=0; j<A.m; j++,ij++)  if(i==j)  A.p[ij]+=x;
  } else

  if(A.sel==1)  for(i=0; i<A.n; i++)  A.p[i]+=x;

  return A;
}

matrix operator+(matrix A, const matrix &B)
{
  int i,j,ij;

  if(A.sel==0 && B.sel==0) {
    if(A.n!=B.n || A.m!=B.m) on_error("matrix","diff. dimensions in A+B (1)");
    for(ij=0; ij<A.nm; ij++)  A.p[ij]+=B.p[ij];
  } else

  if(A.sel==0 && B.sel==1) {
    if(A.n!=B.n || A.m!=B.n) on_error("matrix","diff. dimension in A+B (2)");
    for(i=ij=0; i<A.n; i++)  for(j=0; j<A.m; j++,ij++)
      if(i==j)  A.p[ij]+=B.p[i];
  } else

  if(A.sel==1 && B.sel==1) {
    if(A.n!=B.n) on_error("matrix","different dimension in A+B (3)");
    for(i=0; i<A.n; i++)  A.p[i]+=B.p[i];
  } else

  if(A.sel==1 && B.sel==0)  return B+A;

  return A;
}

matrix operator-(matrix A)
{
  int ij;
  for(ij=0; ij<A.nm; ij++)  A.p[ij]=-A.p[ij];
  return A;
}

matrix operator-(complex x, matrix A)
{
  int i,j,ij;

  if(A.sel==0) {
    if(A.n!=A.m)  on_error("matrix","A is not square in x-A");
    for(i=ij=0; i<A.n; i++)  for(j=0; j<A.m; j++,ij++)
      if(i==j)  A.p[ij]=x-A.p[ij];  else  A.p[ij]=-A.p[ij];
  } else

  if(A.sel==1)  for(i=0; i<A.n; i++)  A.p[i]=x-A.p[i];

  return A;
}

matrix operator-(matrix A, const matrix &B)
{
  int i,j,ij;

  if(A.sel==0 && B.sel==0) {
    if(A.n!=B.n || A.m!=B.m) on_error("matrix","diff. dimensions in A-B (1)");
    for(ij=0; ij<A.nm; ij++)  A.p[ij]-=B.p[ij];
  } else

  if(A.sel==0 && B.sel==1) {
    if(A.n!=B.n || A.m!=B.n) on_error("matrix","diff. dimension in A-B (2)");
    for(i=ij=0; i<A.n; i++)  for(j=0; j<A.m; j++,ij++)
      if(i==j)  A.p[ij]-=B.p[i];
  } else

  if(A.sel==1 && B.sel==1) {
    if(A.n!=B.n) on_error("matrix","different dimension in A-B (3)");
    for(i=0; i<A.n; i++)  A.p[i]-=B.p[i];
  } else

  if(A.sel==1 && B.sel==0)  return -(B-A);

  return A;
}

matrix operator*(complex x, matrix A)
{
  int ij;
  for(ij=0; ij<A.nm; ij++)  A.p[ij]=x*A.p[ij];
  return A;
}

matrix operator*(numero x, matrix A)
{
  int ij;
  for(ij=0; ij<A.nm; ij++)  A.p[ij]=x*A.p[ij];
  return A;
}

matrix operator*(const matrix &A, const matrix &B)
{
  int i,j,k,ij;  complex val;  matrix C;
  complex *Cp,*Ap0,*Ap,*Bp,*Bp0;

  if(A.sel==0 && B.sel==0) {
    if(A.m!=B.n)  on_error("matrix","incompatible dimensions in A*B");
    C.alloc(A.n,B.m);
    for(i=0, Cp=C.p, Ap0=A.p; i<A.n; i++, Ap0+=A.m)
    for(j=0, Bp0=B.p; j<B.m; j++, Cp++, Bp0++) {
      for(k=0, val=0,  Ap=Ap0, Bp=Bp0;  k<A.m; k++, Ap++, Bp+=B.m)
        val+=(*Ap)*(*Bp);
      *Cp=val;
  } } else

  if(A.sel==0 && B.sel==1) {
    if(A.m!=B.n)  on_error("matrix","incompatible dimensions in A*B");
    C.alloc(A.n,A.m);
    for(i=ij=0; i<A.n; i++) for(j=0; j<A.m; j++,ij++) C.p[ij]=A.p[ij]*B.p[j];
  } else

  if(A.sel==1 && B.sel==0) {
    if(A.n!=B.n)  on_error("matrix","incompatible dimensions in A*B");
    C.alloc(B.n,B.m);
    for(i=ij=0; i<B.n; i++) for(j=0; j<B.m; j++,ij++) C.p[ij]=A.p[i]*B.p[ij];
  } else

  if(A.sel==1 && B.sel==1) {
    if(A.n!=B.n)  on_error("matrix","incompatible dimensions in A*B");
    C.alloc_diagonal(A.n);
    for(i=0; i<A.n; i++) C.p[i]=A.p[i]*B.p[i];
  }

  return C;
}

matrix &matrix::operator--(void)        // matrix inversion and return ln|A|
{
  int ij;  if(sel==1) {for(ij=0; ij<nm; ij++) p[ij]=1/p[ij];  return *this;}
  matrix B;
  gaussj(B,n,0);

  return *this;
}

matrix &matrix::operator=(const matrix &A)
{
  int ij;
  if(A.nm)  alloc(A.n,A.m);  sel=A.sel;
  for(ij=0; ij<nm; ij++)  p[ij]=A.p[ij];
  return *this;
}

matrix &matrix::operator=(complex x)
{
  int i,j,ij;

  if(sel==0)
    for(i=ij=0; i<n; i++) for(j=0; j<m; j++,ij++)
      if(i==j) p[ij]=x; else  p[ij]=0;

  else if(sel==1) for(i=0; i<n; i++) p[i]=x;

  return *this;
}

matrix &matrix::operator=(numero x)
{
  int i,j,ij;

  if(sel==0)
    for(i=ij=0; i<n; i++) for(j=0; j<m; j++,ij++)
      if(i==j) p[ij]=x; else  p[ij]=0;

  else if(sel==1) for(i=0; i<n; i++) p[i]=x;

  return *this;
}

matrix operator*(numero *r, matrix A)
{
  int ij;  complex *p;  p=A.p;

  for(ij=0; ij<A.nm; ij++, r++, p++)  *p=(*r)*(*p);

  return A;
}

complex matrix::gaussj(const matrix &Mb, int nn, int mm, int opt)
{
          if(n!=m || n<nn)  on_error("matrix","A not square in gaussj(A)");
  if(mm)  if(m!=Mb.n || Mb.m<mm)
            on_error("matrix","incompatible dimensions in gaussj");

  int  i,icol,irow,j,k,l,ll, *indxc,*indxr,*ipiv,*iA,*iB;
  numero  big,val;
  complex  dum,pivinv, *pirow,*picol,*pp,*ppp,*pa,*pb,*ppb, det=0;

  ipiv=new int [nn];  indxc=new int [nn];  indxr=new int [nn];
  iA=new int [nn];  iB=new int [nn];
  for(i=0; i<nn; i++)  {ipiv[i]=0;  iA[i]=i*m;  iB[i]=i*Mb.m;}

  for(i=0; i<nn; i++) {                    // main loop over columns
    big=0;
    for(j=0,pp=p; j<nn; j++,pp+=m) {       // search pivot element
      if(ipiv[j]!=1)
      for(k=0,ppp=pp; k<nn; k++,ppp++) {
        if(!ipiv[k]) {
          val=mod(*ppp);
          if(val>=big) {
            big=val;
            irow=j;
            icol=k;
        } }  else  if(ipiv[k]>1)
                     on_error("matrix","1) singular matrix in gaussj");
    } }
    ipiv[icol]=ipiv[icol]+1;
    if(irow!=icol) {
      det=det+complex(0,pi);
      pirow=p+iA[irow];  picol=p+iA[icol];
      for(l=0; l<nn; l++,pirow++,picol++) {
        dum=*pirow;  *pirow=*picol;  *picol=dum;
      }
      if(mm)  {
        pirow=Mb.p+iB[irow];  picol=Mb.p+iB[icol];
        for(l=0; l<mm; l++,pirow++,picol++) {
          dum=*pirow;  *pirow=*picol;  *picol=dum;
    } } }
    indxr[i]=irow;
    indxc[i]=icol;   pp=p+iA[icol]+icol;
    if(real(*pp)==0 && imag(*pp)==0) {
      if(opt)  return -infinity;
      else     on_error("matrix","2) singular matrix in gaussj");
    }
    det=det+log(*pp);
    pivinv=1/(*pp);
    *pp=1;
    for(l=0, pp=p+iA[icol]; l<nn; l++,pp++)  *pp=(*pp)*pivinv;
    if(mm)  for(l=0, pp=Mb.p+iB[icol]; l<mm; l++,pp++)  *pp=(*pp)*pivinv;
    for(ll=0,pp=p,ppb=Mb.p; ll<nn; ll++,pp+=m,ppb+=mm)  if(ll!=icol)  {
        dum=*(pp+icol);
        *(pp+icol)=0;  pa=p+iA[icol];  pb=Mb.p+iB[icol];
        for(l=0,ppp=pp; l<nn; l++,ppp++,pa++)  *ppp=*ppp-(*pa)*dum;
        if(mm)  for(l=0,ppp=ppb; l<mm; l++,ppp++,pb++)  *ppp=*ppp-(*pb)*dum;
    }
  }

  for(l=nn-1; l>=0; l--)  if(indxr[l]!=indxc[l])
  for(k=0,pp=p+indxr[l],ppp=p+indxc[l]; k<nn; k++,pp+=m,ppp+=m) {
      dum=*pp;  *pp=*ppp;  *ppp=dum;
  }

  delete [] ipiv;  delete [] indxc;  delete [] indxr;
  delete [] iA;  delete [] iB;

  val=imag(det);  val-=int(val/(2*pi))*2*pi;            // 2*pi phase out
  while(val<-2*pi) val+=2*pi;  while(val>2*pi) val-=2*pi;

  return complex(real(det), val);
}

#endif  // ******************************************************************
