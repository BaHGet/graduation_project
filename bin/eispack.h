#ifndef jga_eispack       // ************************************************
#define jga_eispack 1     // ***  jga/eispack.h                           ***
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
// ***  Subset of EISPACK library to find eigenvalues and eigenvectors    ***
// ***  of arbitrary square complex matrices.                             ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  void eispack(matrix &a, matrix &e, matrix &v);                    ***
// ***                                                                    ***
// **************************************************************************


// **************************************************************************
//    cbabk2 forms the eigenvectors of a complex general matrix by back
//    transforming those of the corresponding balanced matrix determined
//    by cbal
//
//    On input:
//
//       N  IS THE ORDER OF THE MATRIX,
//
//       LOW AND IGH ARE INTEGERS DETERMINED BY cbal,
//
//       SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS
//             AND SCALING FACTORS USED BY  cbal,
//
//       M  IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,
//
//       A  CONTAINS THE EIGENVECTORS TO BE
//          BACK TRANSFORMED IN THEIR FIRST M COLUMNS
//
//    On output:
//
//       A  CONTAINS THE TRANSFORMED EIGENVECTORS IN THEIR FIRST M COLUMNS
//    -----------------------------------------------------------------------

void cbabk2(int N, int M, int LOW, int IGH, numero *SCALE, complex **A)
{
  complex s;
  numero  val;
  int I,J,K,II;

  if(M>0) {
    if(IGH!=LOW)
    for(I=LOW; I<=IGH; I++) {
      val=SCALE[I];
      for(J=1; J<=M; J++) A[I][J]=A[I][J]*val;
    }

    for(II=1; II<=N; II++) {
      I=II;
      if(I<LOW || I>IGH) {
        if(I<LOW) I=LOW-II;
        K=int(SCALE[I]);
        if(K!=I)
        for(J=1; J<=M; J++) {s=A[I][J]; A[I][J]=A[K][J]; A[K][J]=s;}
} } } }

// **************************************************************************
//    THIS SUBROUTINE BALANCES A COMPLEX MATRIX AND ISOLATES
//    EIGENVALUES WHENEVER POSSIBLE
//
//    ON INPUT:
//
//       N IS THE ORDER OF THE MATRIX
//
//       A  CONTAINS THE COMPLEX MATRIX TO BE BALANCED
// 
//    ON OUTPUT:
//
//       A  CONTAINS THE BALANCED MATRIX
//
//       LOW AND IGH ARE TWO INTEGERS SUCH THAT A[I][J]
//          ARE EQUAL TO ZERO IF
//          (1) I IS GREATER THAN J AND
//          (2) J=1,...,LOW-1 OR I=IGH+1,...,N
//
//       SCALE CONTAINS INFORMATION DETERMINING THE
//             PERMUTATIONS AND SCALING FACTORS USED
//
//    SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
//    HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
//    WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
//    OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J). THEN
//       SCALE(J) = P(J),    FOR J = 1,...,LOW-1
//                = D(J,J)       J = LOW,...,IGH
//                = P(J)         J = IGH+1,...,N.
//    THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
//    THEN 1 TO LOW-1.
//
//    NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY
//    -----------------------------------------------------------------------

void cbal(int N, complex **A, int &LOW, int &IGH, numero *SCALE)
{
  int I,J,K,L,M,JJ,IEXC, NOCONV;
  complex FF;
  numero C,G,F,R,S,B2,RADIX;

  RADIX=2;  B2=RADIX*RADIX;  K=1;  L=N;

  goto gt100;

gt20:
  SCALE[M]=J;
  if(J!=M) {
    for(I=1; I<=L; I++) {FF=A[I][J];  A[I][J]=A[I][M];  A[I][M]=FF;}
    for(I=K; I<=N; I++) {FF=A[J][I];  A[J][I]=A[M][I];  A[M][I]=FF;}
  }

  if(IEXC==2)  K++;  else  {
    if(L==1) goto gt280;
    L--;

gt100:
    for(JJ=1; JJ<=L; JJ++) {
      J=L+1-JJ;
      for(I=1; I<=L; I++)
        if(I!=J)  if(real(A[J][I])!=0 || imag(A[J][I])!=0) goto gt120;
      M=L;  IEXC=1;  goto gt20;
gt120: ;
  } }

  for(J=K; J<=L; J++) {
    for(I=K; I<=L; I++)
      if(I!=J)  if(real(A[I][J])!=0 || imag(A[I][J])!=0) goto gt170;
    M=K;  IEXC=2;
    goto gt20;
gt170: ;
  }

  for(I=K; I<=L; I++) SCALE[I]=1;

  do {NOCONV=0;

    for(I=K; I<=L; I++) {C=R=0;
      for(J=K; J<=L; J++) if(J!=I) {
        C+=ABS(real(A[J][I]))+ABS(imag(A[J][I]));    // jga
        R+=ABS(real(A[I][J]))+ABS(imag(A[I][J]));    // jga
	//        C+=mod(A[J][I]);  R+=mod(A[I][J]);
      }

      if(C!=0 && R!=0) {
        G=R/RADIX;
        F=1;  S=C+R;

        do {
          if(C>=G) goto gt220;
          F=F*RADIX;
          C=C*B2;
        } while(1);

gt220:  G=R*RADIX;
        do {
          if(C<G) goto gt240;
          F=F/RADIX;
          C=C/B2;
        } while(1);

gt240:  if((C+R)/F<0.95*S) {
          G=1/F;
          SCALE[I]=SCALE[I]*F;
          NOCONV=1;

          for(J=K; J<=N; J++) A[I][J]=A[I][J]*G;
          for(J=1; J<=L; J++) A[J][I]=A[J][I]*F;

  } } } } while(NOCONV);

gt280:  LOW=K;  IGH=L;
}

// **************************************************************************
//    GIVEN A  COMPLEX  GENERAL  MATRIX, THIS  SUBROUTINE  REDUCES  A 
//    SUBMATRIX SITUATED IN ROWS AND COLUMNS LOW THROUGH IGH TO UPPER 
//    HESSENBERG FORM BY STABILIZED ELEMENTARY SIMILARITY TRANSFORMS.  
// 
//    ON INPUT:
//
//       N        IS THE ORDER OF THE MATRIX
//       LOW,IGH  ARE INTEGERS DETERMINED BY THE BALANCING SUBROUTINE
//                cbal. IF  cbal  HAS NOT BEEN USED, SET LOW=1, IGH=N
//       A        CONTAINS THE COMPLEX INPUT MATRIX
//
//    ON OUTPUT:
//
//       A        CONTAINS THE HESSENBERG MATRIX.THE MULTIPLIERS WHICH WERE
//                USED IN THE  REDUCTION  ARE STORED IN THE REMAINING
//                TRIANGLES UNDER THE HESSENBERG MATRIX,
//       III      CONTAINS INFORMATION ON THE ROWS AND COLUMNS INTER-
//                CHANGED IN THE REDUCTION. ONLY ELEMENTS LOW THROUGH
//                IGH ARE USED
//
//    ARITHMETIC IS COMPLEX
//    -----------------------------------------------------------------------

int comhes(int N, int LOW, int IGH, complex **A, int *III)
{
  int I,J,M,LA,KP1,MM1,MP1;
  complex  X,Y,Z3;

  LA=IGH-1;  KP1=LOW+1;
  if(LA<KP1)  return 0;

  for(M=KP1; M<=LA; M++) {
    MM1=M-1;
    X=0;  I=M;

    for(J=M; J<=IGH; J++)
    //      if(ABS(real(A[J][MM1]))+ABS(imag(A[J][MM1]))>           // jga
    //           ABS(real(X))+ABS(imag(X))) {X=A[J][MM1];  I=J;}    // jga
      if(mod2(A[J][MM1])>mod2(X)) {X=A[J][MM1];  I=J;}
    III[M]=I;
    if(I!=M) {
      for(J=MM1; J<=N; J++)  {Y=A[I][J];  A[I][J]=A[M][J];  A[M][J]=Y;}
      for(J=1; J<=IGH; J++)  {Y=A[J][I];  A[J][I]=A[J][M];  A[J][M]=Y;}
    }
    if(real(X)!=0 || imag(X)!=0) {
      MP1=M+1;
      for(I=MP1; I<=IGH; I++) {     // heaviest part here
        Y=A[I][MM1];
        if(real(Y)!=0 || imag(Y)!=0) {
          Z3=Y/X;  Y=Z3;  A[I][MM1]=Y;
	  for(J=M; J<=N; J++)  A[I][J]-=Y*A[M][J];
          for(J=1; J<=IGH; J++)  A[J][M]+=Y*A[J][I];
  } } } }

  return 0;
}

// **************************************************************************
//    THIS SUBROUTINE FINDS  THE  EIGENVALUES AND  EIGENVECTORS  OF  A
//    COMPLEX UPPER HESSENBERG  MATRIX BY THE MODIFIED  LR METHOD. THE
//    EIGENVECTORS  OF A COMPLEX  GENERAL MATRIX  CAN ALSO BE FOUND IF
//    comhes HAS BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG
//    FORM
//
//    ON INPUT:
//
//       N       IS THE ORDER OF THE MATRIX
//       LOW,IGH ARE INTEGERS DETERMINED BY THE  BALANCING  SUBROUTINE
//               cbal.  IF  cbal  HAS NOT BEEN USED,  SET LOW=1, IGH=N
//       INT     CONTAINS INFORMATION ON THE ROWS AND  COLUMNS  INTER-
//               CHANGED IN THE REDUCTION BY comhes,IF PERFORMED. ONLY
//               ELEMENTS LOW THROUGH IGH ARE USED.IF THE EIGENVECTORS
//               OF THE HESSENBERG MATRIX ARE DESIRED,SET III[J]=J FOR
//               THESE ELEMENTS
//       H       CONTAINS THE COMPLEX UPPER HESSENBERG MATRIX. THEIR LOWER
//               TRIANGLES  BELOW THE SUBDIAGONAL CONTAIN THE MULTIPLIERS
//               WHICH   WERE  USED  IN  THE  REDUCTION BY  comhes, IF
//               PERFORMED.  IF  THE  EIGENVECTORS  OF  THE HESSENBERG
//               MATRIX ARE DESIRED,THESE ELEMENTS MUST BE SET TO ZERO
//
//     ON OUTPUT:
//
//       H       THE   UPPER HESSENBERG PORTIONS OF H HAVE BEEN
//               DESTROYED, BUT  THE LOCATION H(1,1) CONTAINS THE NORM
//               OF THE TRIANGULARIZED MATRIX,
//       W       CONTAINS THE EIGENVALUES.
//       Z       CONTAINS THE EIGENVECTORS.THE EIGENVECTORS ARE UNNORMALIZED.
//               IF AN ERROR EXIT IS  MADE, NONE OF THE  EIGENVECTORS  HAS
//               HAS BEEN FOUND
//
//       Returns 0 for normal execution or a number j corresponding to an
//               eigenvalue that has not been determined after 200 iterations
//               and the eigenvalues should be correct for indices j+1,...,N
//    -----------------------------------------------------------------------

numero MACHEP_int=25;

int comlr2(int N, int LOW, int IGH, int *III,
           complex **H, complex *W, complex **Z)
{
  int I,J,K,L,M,EN,II,JJ,LL,MM,NN,IM1,IP1,ITS,MP1,ENM1,IEND;
  numero NORM, MACHEP=pow(2.0,-MACHEP_int); // old value MACHEP=pow(2.0,-47.0);
  complex S,T,X,Y,ZZ,Z3;

  for(I=1; I<=N; I++)
  for(J=1; J<=N; J++)  if(I==J) Z[I][I]=1; else Z[I][J]=0;

  IEND=IGH-LOW-1;
  if(IEND>0)
    for(II=1; II<=IEND; II++) {
      I=IGH-II;
      IP1=I+1;
      for(K=IP1; K<=IGH; K++) Z[K][I]=H[K][I-1];
      J=III[I];
      if(I!=J) {
        for(K=I; K<=IGH; K++) {Z[I][K]=Z[J][K];  Z[J][K]=0;}
        Z[J][I]=complex(1, imag(Z[J][I]));
    } }

  for(I=1; I<=N; I++)  if(I<LOW || I>IGH) W[I]=H[I][I];

  EN=IGH;  T=0;

gt220:
  if(EN<LOW) goto gt680;
  ITS=0;
  ENM1=EN-1;

  do {

    for(LL=LOW; LL<=EN; LL++) {
      L=EN+LOW-LL;
      if(L==LOW) goto gt300;
      if(ABS(real(H[L][L-1]))+ABS(imag(H[L][L-1]))<=            // jga
         MACHEP*(  ABS(real(H[L-1][L-1]))+ABS(imag(H[L-1][L-1]))
                 + ABS(real(H[L][L]))+ABS(imag(H[L][L])))) goto gt300;
      // if(mod(H[L][L-1])<=MACHEP*(mod(H[L-1][L-1])+mod(H[L][L]))) goto gt300;
    }

gt300:
    if(L==EN) goto gt660;
    if(ITS==200) return EN;
    if(ITS==0 || ITS%10!=0) {
      S=H[EN][EN];
      X=H[ENM1][EN]*H[EN][ENM1];
      if(real(X)==0 && imag(X)==0) goto gt340;
      Y=(H[ENM1][ENM1]-S)/2;
      ZZ=Z3=sqrt(Y*Y+X);
      if(real(Y)*real(ZZ)+imag(Y)*imag(ZZ)<0)  ZZ=-ZZ;
      Z3=X/(Y+ZZ);
      S-=Z3;
      goto gt340;
    }

    S=complex(ABS(real(H[EN][ENM1]))+ABS(real(H[ENM1][EN-2])),
              ABS(imag(H[EN][ENM1]))+ABS(imag(H[ENM1][EN-2])));

gt340:
    for(I=LOW; I<=EN; I++) H[I][I]-=S;

    T+=S;  ITS++;

    X=complex(ABS(real(H[ENM1][ENM1]))+ABS(imag(H[ENM1][ENM1])),imag(X));
    Y=complex(ABS(real(H[EN][ENM1]))+ABS(imag(H[EN][ENM1])), imag(Y));
    ZZ=complex(ABS(real(H[EN][EN]))+ABS(imag(H[EN][EN])), imag(ZZ));
    for(MM=L; MM<=ENM1; MM++) {
      M=ENM1+L-MM;
      if(M==L) goto gt420;
      Y=complex(ABS(real(H[M][M-1]))+ABS(imag(H[M][M-1])), real(Y));
      X=complex(real(X), real(ZZ));
      ZZ=complex(real(X), imag(ZZ));
      X=complex(ABS(real(H[M-1][M-1]))+ABS(imag(H[M-1][M-1])), imag(X));
      if(real(Y)<=MACHEP*real(ZZ)/imag(Y)*(real(ZZ)+real(X)+imag(X)))
        goto gt420;
    }

gt420:
    MP1=M+1;

    for(I=MP1; I<=EN; I++) {
      IM1=I-1;
      X=H[IM1][IM1];  Y=H[I][IM1];
      if(ABS(real(X))+ABS(imag(X))<ABS(real(Y))+ABS(imag(Y))) { // jga
	//      if(mod(X)<mod(Y)) {
        for(J=IM1; J<=N; J++) {ZZ=H[IM1][J];  H[IM1][J]=H[I][J];  H[I][J]=ZZ;}
        Z3=X/Y;
        W[I]=complex(1, imag(W[I]));
        goto gt480;
      }
      Z3=Y/X;
      W[I]=complex(-1, imag(W[I]));

gt480:
      H[I][IM1]=ZZ=Z3;
      for(J=I; J<=N; J++)  H[I][J]-=ZZ*H[IM1][J];
    }

    for(J=MP1; J<=EN; J++) {        // most time is consumed in this loop
      X=H[J][J-1];  H[J][J-1]=0;
      if(real(W[J])>0) {
        for(I=1; I<=J; I++)     {ZZ=H[I][J-1]; H[I][J-1]=H[I][J]; H[I][J]=ZZ;}
        for(I=LOW; I<=IGH; I++) {ZZ=Z[I][J-1]; Z[I][J-1]=Z[I][J]; Z[I][J]=ZZ;}
      }
      for(I=1; I<=J; I++)  H[I][J-1]+=X*H[I][J];
      for(I=LOW; I<=IGH; I++)  Z[I][J-1]+=X*Z[I][J];
    }

  } while(1);

gt660:
  H[EN][EN]+=T;
  W[EN]=H[EN][EN];
  EN=ENM1;
  goto gt220;

gt680:
  NORM=0;
  for(I=1; I<=N; I++)  for(J=I; J<=N; J++)
    NORM+=ABS(real(H[I][J]))+ABS(imag(H[I][J]));  // jga
  //    NORM+=mod(H[I][J]);

  H[1][1]=NORM;
  if(N==1 || NORM==0)  return 0;

  for(NN=2; NN<=N; NN++) {        // N^3, only once?
    EN=N+2-NN;
    X=W[EN];
    ENM1=EN-1;
    for(II=1; II<=ENM1; II++) {
      I=EN-II;
      ZZ=H[I][EN];
      if(I!=ENM1) {
        IP1=I+1;
        for(J=IP1; J<=ENM1; J++)  ZZ+=H[I][J]*H[J][EN];
      }
      Y=X-W[I];
      if(real(Y)==0 && imag(Y)==0) Y=complex(MACHEP*NORM, imag(Y));
      H[I][EN]=Z3=ZZ/Y;
  } }
  
  ENM1=N-1;

  for(I=1; I<=ENM1; I++) {
    if(I>=LOW && I<=IGH) goto gt840;
    IP1=I+1;
    for(J=IP1; J<=N; J++)  Z[I][J]=H[I][J];
  }

gt840:
  for(JJ=LOW; JJ<=ENM1; JJ++) {        // N^3, only once, but heavy
    J=N+LOW-JJ;
    M=(J-1<IGH)?J-1:IGH;
    for(I=LOW; I<=IGH; I++) {
      ZZ=Z[I][J];
      for(K=LOW; K<=M; K++)  ZZ+=Z[I][K]*H[K][J];
      Z[I][J]=ZZ;
  } }

  return 0;
}

// **************************************************************************
//    THIS   SUBROUTINE  COMPUTES  ALL  EIGENVALUES  AND  CORRESPONDING
//    EIGENVECTORS  OF  AN  ARBITRARY   COMPLEX  MATRIX.  THE MATRIX IS
//    BALANCED BY EXACT NORM  REDUCING  SIMILARITY  TRANSFORMATIONS AND
//    THEN  IS  REDUCED  TO  COMPLEX  HESSENBERG   FORM  BY  STABILIZED
//    ELEMENTARY SIMILARITY TRANSFORMATIONS. A MODIFIED LR ALGORITHM IS
//    USED TO COMPUTE THE EIGENVALUES OF THE HESSENBERG MATRIX.
//
//      ON INPUT:
//
//         N        IS THE ORDER OF THE MATRIX. N MAY BE 1.
//         A        NxN ARBITRARY COMPLEX MATRIX WHOSE EIGENSYSTEM IS TO BE
//                  COMPUTED.
//
//       ON OUTPUT:
//
//         EV         CONTAIN THE REAL AND IMAGINARY PARTS RESPECTIVELY
//                    OF THE COMPUTED EIGENVALUES.  THE EIGENVALUES ARE
//                    NOT ORDERED IN ANY WAY.
//         VEC        CONTAINS IN THE LEADING N BY N  SUBARRAYS THE COMPUTED
//                    EIGENVECTORS.  THE J-TH COLUMNS  OF VEC CONTAIN THE
//                    EIGENVECTOR  ASSOCIATED  WITH EV(J). THE EIGENVECTORS
//                    ARE NOT NORMALIZED IN ANY WAY.
//         It stops if the j-TH EIGENVALUE HAS NOT BEEN FOUND IN
//                    200 ITERATIONS. THE FIRST J-1 ELEMENTS OF EV
//                    CONTAIN THOSE EIGENVALUES ALREADY FOUND. NO
//         A          IS DESTROYED.
//    -----------------------------------------------------------------------

void cnaa(int N, complex **A, complex *EV, complex **VEC)
{
  int IGH,LOW;
  int *III, j;     III=new int [N+1];
  numero * SCALE;   SCALE=new numero [N+1];

  cbal(N,A,LOW,IGH,SCALE);
  comhes(N,LOW,IGH,A,III);

  if(j=comlr2(N,LOW,IGH,III,A,EV,VEC))
    on_error("cnaa", "eigenvalue", j, "not found in 200 iterations");

  cbabk2(N,N,LOW,IGH,SCALE,VEC);

  delete [] III;  delete [] SCALE;
}

// **************************************************************************
//    eispack(a,e,v) computes all the eigenvalues and eigenvectors of the
//    square matrix a. On input, a is an arbitrary square complex matrix.
//    On output, a is unaltered, e is the array of eigenvalues corresponding
//    to the normalized eigenvectors given by the columns of v.
//    Dimensions:   a(n,n), v(n,n), e(n)=e(n,1)
//    -----------------------------------------------------------------------

void eispack(numero *a, numero *v, numero *e, int n)
{
  complex **aa, *ee, **vv;
  numero  norm;
  int     i,j;

  aa=new complex* [n+1];  ee=new complex [n+1];  vv=new complex* [n+1];
  for(i=0; i<n; i++) {
    aa[i+1]=new complex [n+1];
    vv[i+1]=new complex [n+1];
    for(j=0; j<n; j++)  aa[i+1][j+1]=a[i*n+j];
  }

  cnaa(n,aa,ee,vv);

  for(i=0; i<n; i++) {  // store eigensystem in e and v
    e[i]=real(ee[i+1]);
    for(j=0; j<n; j++)  v[i*n+j]=real(vv[i+1][j+1]);
    delete [] aa[i+1];  delete [] vv[i+1];
  }

 for(j=0; j<n; j++) {   // normalize eigenvectors
    norm=0;
    for(i=0; i<n; i++)  norm+=sqr(v[i*n+j]);  norm=sqrt(norm);
    for(i=0; i<n; i++)  v[i*n+j]=v[i*n+j]/norm;
  }

  delete [] aa;  delete [] ee;  delete [] vv;
}

#ifdef jga_matrix

void eispack(matrix &a, matrix &e, matrix &v)
{
  complex **aa, *ee, **vv;
  numero  norm;
  int     n=a.n, i,j;

  if(a.n<1 || a.n!=a.m)  on_error("eispack", "wrong intput dimensions");

  aa=new complex* [n+1];  ee=new complex [n+1];  vv=new complex* [n+1];
  for(i=0; i<n; i++) {
    aa[i+1]=new complex [n+1];
    vv[i+1]=new complex [n+1];
    for(j=0; j<n; j++)  aa[i+1][j+1]=a.a(i,j);
  }

  cnaa(n,aa,ee,vv);

  e.alloc(n);  v.alloc(n,n);
  for(i=0; i<n; i++) {  // store eigensystem in e and v
    e.a(i,ee[i+1]);
    for(j=0; j<n; j++)  v.a(i,j,vv[i+1][j+1]);
    delete [] aa[i+1];  delete [] vv[i+1];
  }

 for(j=0; j<n; j++) {   // normalize eigenvectors
    norm=0;
    for(i=0; i<n; i++)  norm+=mod2(v.a(i,j));  norm=sqrt(norm);
    for(i=0; i<n; i++)  v.a(i,j, v.a(i,j)/norm);
  }

  delete [] aa;  delete [] ee;  delete [] vv;
}

void eispack_order(matrix &e, matrix &v, int sel)
{
  matrix ee,vv;
  numero r;
  int    n=v.n, i,j,l;
  int    *o,*oo;                          // sel=8 -> small  Re{e}
                                          // sel=1 -> small |Re{e}|
  o=new int [n];  oo=new int [n];         // sel=2 -> small |Im{e}|
  for(j=0; j<n; j++) oo[j]=1;             // sel=3 -> small |e|
                                          // sel=4 -> large |Re{e}|
  for(i=0; i<n; i++) {                    // sel=5 -> large |Im{e}|
    r=infinity;                           // sel=6 -> large |e|
    for(j=0; j<n; j++) if(oo[j])          // sel=7 -> large Im{e}
      switch(sel) {
        case 1: if( ABS(real(e.a(j)))<r) {r= ABS(real(e.a(j))); l=j;}  break;
        case 2: if( ABS(imag(e.a(j)))<r) {r= ABS(imag(e.a(j))); l=j;}  break;
        case 3: if( ABS(mod2(e.a(j)))<r) {r= ABS(mod2(e.a(j))); l=j;}  break;
        case 4: if(-ABS(real(e.a(j)))<r) {r=-ABS(real(e.a(j))); l=j;}  break;
        case 5: if(-ABS(imag(e.a(j)))<r) {r=-ABS(imag(e.a(j))); l=j;}  break;
        case 6: if(-ABS(mod2(e.a(j)))<r) {r=-ABS(mod2(e.a(j))); l=j;}  break;
        case 7: if(     -imag(e.a(j))<r) {r=    -imag(e.a(j));  l=j;}  break;
        case 8: if(      real(e.a(j))<r) {r=     real(e.a(j) ); l=j;}  break;
      }
    o[i]=l;  oo[l]=0;
  }

  ee=e;  vv=v;
  for(j=0; j<n; j++) {
    e.a(j, ee.a(o[j]));
    for(i=0; i<n; i++)  v.a(i,j, vv.a(i,o[j]));
  }

  ee.free();  vv.free();  delete [] oo;  delete [] o;
}

#endif  // matrix

#endif  // ******************************************************************
