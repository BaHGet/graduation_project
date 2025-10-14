#ifndef jga_plot          // ************************************************
#define jga_plot 1        // ***  jga/plot.h           1999 -> 14-iii-00  ***
                          // ***                                          ***
#include "jga.h"          // ***                                          ***
#include "gd/a.h"         // ***  GIF generation routines                 ***
                          // ***                                          ***
// ***  ----------------  // ***  --------------------------------------  ***
// ***                                                                    ***
// ***  Author: F. J. Garcia de Abajo, CSIC, Spain                        ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  These routines make use of the gd library to plot functions       ***
// ***                                                                    ***
// ***  The corners of the image are:  upper-left  (0,0)                  ***
// ***                                 upper-right (dimx-1,0)             ***
// ***                                 lower-left  (0,dimy-1)             ***
// ***                                 lower-right (dimx-1,dimy-1)        ***
// ***  in physical pixels.                                               ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  plot_curves(*name, *x, *y, n, *ni, *ci, dimx, dimy,               ***
// ***              *labelx, *labely);                                    ***
// ***  plot_curves(*name, *x, *y, n,  ni, *ci, dimx, dimy,               ***
// ***              *labelx, *labely);                                    ***
// ***  plot_polar(*name, *r, *phi, n,  ni, *ci, dim, maxr, axis,frame);  ***
// ***  plot_polar(*name, *r, *phi, n,  ni, *ci, dim);                    ***
// ***  plot_cluster(*name, *x, *y, *z, *t, n, *ai, *ci, th, phi,         ***
// ***               dimx, dimy, &minx, &maxx, &miny, &maxy, axis, size); ***
// ***  plot_cluster(*name, *x, *y, *z, *t, n, *ai, *ci, th, phi,         ***
// ***               dimx, dimy, axis, size);                             ***
// ***  plot_cluster(*name, *x, *y, *z, *t, n, *ai, *ci, th, phi,         ***
// ***               dimx, dimy);   // print axis by default              ***
// ***  plot_hologram(*name, *intensity, tth0, tth1, th0, th1, nth,       ***
// ***                                   fi0, fi1, nfi, dim);             ***
// ***  plot_hologram(*name, *intensity, tth0, tth1, th0, th1, nth,       ***
// ***                                   fi0, fi1, nfi, dim, ilog);       ***
// ***                                                                    ***
// ***  int plot_color(*color);    // return color index                  ***
// ***  int plot_color(gdImagePtr &plot_graph, int color);                ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***      *name      output file name                                   ***
// ***      *x,*y,*z   array with x,y,z coordinates                       ***
// ***      n          number of curves                                   ***
// ***      *ni        number of points in each curve (one *x per curve)  ***
// ***      ni         number of points (all curves share same *x)        ***
// ***      *ci        colors (one per curve, atom, etc.)                 ***
// ***      dimx,dimy  dimension of image (pixels in hor/ver direction)   ***
// ***      *label?    axis labels                                        ***
// ***      dim        square image with dim pixels in each side          ***
// ***      min?,max?  axis limits                                        ***
// ***      th,phi     point of observation for 3D plots                  ***
// ***      th0,th1,nth   polar angle range and number of angles          ***
// ***      fi0,fi1,nfi   azimuthal angle range and number of angles      ***
// ***      *intensity    array with (th,phi) intensities                 ***
// ***      tth1,tth0     range of representation                         ***
// ***      axis       print axes?                                        ***
// ***                                                                    ***
// **************************************************************************

// **************************************************************************
// *** plot_color                                                         ***
// **************************************************************************

int plot_color(char *name)
{
  int a=256;

  if(!strcmp(name,"black"))     return  0;
  if(!strcmp(name,"white"))     return  a*a*a-1;
  if(!strcmp(name,"red"))       return  255*a*a;
  if(!strcmp(name,"green"))     return  255*a;
  if(!strcmp(name,"blue"))      return  255;
  if(!strcmp(name,"violet"))    return  153*a*a+255;
  if(!strcmp(name,"yellow"))    return  (255*a+255)*a;
  if(!strcmp(name,"magenta"))   return  255*a*a+255;
  if(!strcmp(name,"orange"))    return  (255*a+153)*a;   return 0;
}

// --------------------------------------------------------------------------

int plot_color(gdImagePtr &plot_graph, int color)
{
  int r,g,b;   r=color/(65536);  g=color%65536;  b=g%256;  g=g/256;

  return gdImageColorAllocate(plot_graph, r, g, b);
}

// **************************************************************************
// *** plot_curves -> plots functions in gif format                       ***
// *** (*x,*y) contains the points of each curve in consecutive order     ***
// *** (one curve after another, and for each, one point after another)   ***
// **************************************************************************

void set_maxmin(numero &ma, numero &mi)
{
  numero ama=ABS(ma), ami=ABS(mi);

  if(ma*mi>0) if(ami<ama/6) mi=0; else
              if(ama<ami/6) ma=0;
}

// --------------------------------------------------------------------------

void plot_curves(char *name,      // output file
          numero *x, numero *y,   // points in curves
          int n,                  // number of curves
          int *ni,                // number of points in each curve
          int *ci,                // colors of the curves (from 216; <0 thick)
          int dimx, int dimy,     // physical dimension of image in pixels
          char *labelx,           // x-axis label
          char *labely)           // y-axis label
{
  int i,j;
  numero *xx,*yy;  int xv,yv, xv0,yv0;
  numero minx=infinity,maxx=-infinity,miny=infinity,maxy=-infinity;
  int smar=12, lmar=24;   // margins

  // --- init ---

  gdImagePtr plot_graph;
  plot_graph = gdImageCreate(dimx, dimy);  // physical dimensions of graph

  // --- colors ---

  int *colors;  colors=new int [n];

  int white = gdImageColorAllocate(plot_graph, 255, 255, 255);
  int black = gdImageColorAllocate(plot_graph,   0,   0,   0);

  for(i=0; i<n; i++)  colors[i] = plot_color(plot_graph, ci[i]);

  gdImageColorTransparent(plot_graph, white);

  // --- limits of x and y ---

  xx=x; yy=y;
  for(i=0; i<n; i++)
  for(j=0; j<ni[i]; j++) {
    if(*xx<minx) minx=*xx;
    if(*xx>maxx) maxx=*xx;
    if(*yy<miny) miny=*yy;
    if(*yy>maxy) maxy=*yy;
    xx++;  yy++;
  }

  set_maxmin(maxx,minx);
  set_maxmin(maxy,miny);

  // --- plotting axes ---

  gdImageLine(plot_graph,lmar,     smar,     dimx-smar,smar,     black);
  gdImageLine(plot_graph,dimx-smar,smar,     dimx-smar,dimy-lmar,black);
  gdImageLine(plot_graph,dimx-smar,dimy-lmar,lmar,     dimy-lmar,black);
  gdImageLine(plot_graph,lmar,     dimy-lmar,lmar,     smar,     black);

  dimx-=lmar+smar;  dimy-=lmar+smar;  // dimx, dimy refer to the interior

  // --- grid lines ---

  if(maxy>0 && miny<0) {
    yv=int(floor(smar+dimy+dimy*miny/(maxy-miny)));
    gdImageDashedLine(plot_graph,lmar,yv,lmar+dimx,yv,black);
  }

  if(maxx>0 && minx<0) {
    xv=int(floor(lmar-dimx*minx/(maxx-minx)));
    gdImageDashedLine(plot_graph,xv,smar,xv,smar+dimy,black);
  }

  // --- axis labels ---

  gdImageString(plot_graph, gdFontMediumBold,
                     lmar+(dimx-7*strlen(labelx))/2,
                     dimy+smar+smar,
                     (unsigned  char *) labelx, black);
  gdImageStringUp(plot_graph, gdFontMediumBold,
                     0, smar+(dimy+7*strlen(labely))/2,
                     (unsigned  char *) labely, black);
  gdImageString(plot_graph, gdFontMediumBold, lmar, dimy+smar,
                     (unsigned  char *) ftostr(minx,2), black);
  gdImageString(plot_graph,gdFontMediumBold,
                     dimx+lmar-6*strlen(ftostr(maxx,2)), dimy+smar,
                     (unsigned  char *) ftostr(maxx,2), black);
  gdImageStringUp(plot_graph, gdFontMediumBold,
                     smar, dimy+smar,
                     (unsigned  char *) ftostr(miny,2), black);
  gdImageStringUp(plot_graph, gdFontMediumBold,
                     smar, smar+6*strlen(ftostr(maxy,2)),
                     (unsigned  char *) ftostr(maxy,2), black);

  // --- plotting functions ---

  xx=x; yy=y;

  for(i=0; i<n; i++)
  for(j=0; j<ni[i]; j++) {
    xv=int(floor(lmar+     dimx*((*xx)-minx)/(maxx-minx)));
    yv=int(floor(smar+dimy-dimy*((*yy)-miny)/(maxy-miny)));
    if(j>0) gdImageLine(plot_graph, xv0,yv0,xv,yv, colors[i]);
    xv0=xv;  yv0=yv;
    xx++;  yy++;
  }

  // --- output ---

  gdImageInterlace(plot_graph, 1);    // interlace when viewing

  FILE *fout;
  fout=fopen(name, "wb");
  gdImageGif(plot_graph, fout);
  fclose(fout);

  // --- ending ---

  gdImageDestroy(plot_graph);
  delete [] colors;
}

// **************************************************************************
// *** the same as above, but all curves with the same points xi          ***
// **************************************************************************

void plot_curves(char *name,      // output file
          numero *x, numero *y,   // points in curves
          int n,                  // number of curves
          int ni,                 // number of points in each curve
          int *ci,                // colors of the curves (from 216; <0 thick)
          int dimx, int dimy,     // physical dimension of image in pixels
          char *labelx,           // x-axis label
          char *labely)           // y-axis label
{
  numero *xx,*yy;
  int *nni,i,j;

  nni=new int [n];
  xx=new numero [n*ni];  yy=new numero [n*ni];

  for(i=0; i<n; i++) {
    nni[i]=ni;
    for(j=0; j<ni; j++) {
      xx[i*ni+j]=x[j];  yy[i*ni+j]=y[i*ni+j];
  } }

  plot_curves(name, xx,yy, n,nni,ci, dimx,dimy, labelx,labely);
}

// **************************************************************************
// *** plot_polar -> make polar plot                                      ***
// *** (*r,*phi) contains the polar points of each curve,                 ***
// *** one curve after another, and for each, one point after another     ***
// **************************************************************************

void plot_polar(char *name,       // output file
          numero *r, numero *phi, // points in curves
          int n,                  // number of curves
          int ni,                 // number of points in each curve
          int *ci,                // colors of the curves (from 216; <0 thick)
          int dim,                // physical dimension of image in pixels
          numero maxr,            // maximum r
          int axis, int frame)    // plot axis and frame? (1=yes, 0=no)
{
  int i,j;
  numero *rr,*fi;  int xv,yv,xv0,yv0;
  if(dim%2)  dim--;

  // --- init ---

  gdImagePtr plot_graph;
  plot_graph = gdImageCreate(dim+1, dim+1);  // physical dimensions of graph

  // --- colors ---

  int *colors;  colors=new int [n];

  int white = gdImageColorAllocate(plot_graph, 255, 255, 255);
  int black = gdImageColorAllocate(plot_graph,   0,   0,   0);

  for(i=0; i<n; i++)  colors[i] = plot_color(plot_graph, ci[i]);

  gdImageColorTransparent(plot_graph, white);

  // --- plotting axes ---

  if(axis) {
    gdImageLine(plot_graph,0,dim/2,dim,dim/2, black);   // x-axis
    gdImageLine(plot_graph,dim/2,0,dim/2,dim, black);   // y-axis
  }
  if(frame) {
    gdImageLine(plot_graph,0,dim,dim,dim, black);       // frame
    gdImageLine(plot_graph,0,0,dim,0, black);
    gdImageLine(plot_graph,0,0,0,dim, black);
    gdImageLine(plot_graph,dim,0,dim,dim, black);
  }

  // --- plotting functions ---

  rr=r;

  for(i=0; i<n; i++) {
    fi=phi;
    for(j=0; j<ni; j++) {
      xv=int(dim/2+floor(dim/2*cos(*fi)*(*rr)/(1.1*maxr)));
      yv=int(dim/2-floor(dim/2*sin(*fi)*(*rr)/(1.1*maxr)));
      if(j>0) gdImageLine(plot_graph, xv0,yv0,xv,yv, colors[i]);
      xv0=xv;  yv0=yv;
      rr++;  fi++;
  } }

  // --- output ---

  gdImageInterlace(plot_graph, 1);    // interlace when viewing

  FILE *fout;
  fout=fopen(name, "wb");
  gdImageGif(plot_graph, fout);
  fclose(fout);

  // --- ending ---

  gdImageDestroy(plot_graph);
  delete [] colors;
}

// **************************************************************************
// *** plot_polar -> make polar plot                                      ***
// *** (*r,*phi) contains the polar points of each curve,                 ***
// *** one curve after another, and for each, one point after another     ***
// **************************************************************************

void plot_polar(char *name,       // output file
          numero *r, numero *phi, // points in curves
          int n,                  // number of curves
          int ni,                 // number of points in each curve
          int *ci,                // colors of the curves (from 216; <0 thick)
          int dim)                // physical dimension of image in pixels
{
  int i,j;
  numero *rr;
  numero maxr=-infinity;

  rr=r;
  for(i=0; i<n; i++)
  for(j=0; j<ni; j++)  {if(*rr>maxr) maxr=*rr;  rr++;}

  plot_polar(name,r,phi,n,ni,ci,dim, maxr, 1,0);
}

// **************************************************************************
// *** plot_cluster -> plots clusters of spherical components             ***
// *** (*x,*y,*z) contains the centers of the cluster points              ***
// *** (one curve after another, and for curve, one point after another)  ***
// **************************************************************************

void plot_cluster(char *name,     // output file
          numero *x, numero *y,   // objects at points (x,y,z)
          numero *z, int *t,      // t is the type of atom
          int n,                  // number of cluster elements
          numero *ai,             // radius of elements type i
          int *ci,                // colors of i (from 216; <0 thick)
          numero th, numero phi,  // direction from which it is observed
          int dimx, int dimy,     // physical dimension of image in pixels
          numero &minx,           // minimum value of x in the window
          numero &maxx,           // maximum value of x in the window
          numero &miny,           // minimum value of y in the window
          numero &maxy,           // maximum value of y in the window
          int axis,               // 1 to print axes (0 otherwise)
          int size)               // 1 for atom size decrease (0 otherwise)
{
  int i,j,k,l, i0,j0,k0;
  numero *xx,*yy,*zz, x0,y0,z0,n1x,n1y,n1z,n2x,n2y,n2z, val, xv,yv;
  numero bx[2],by[2],bz[2], dd,d, maxz,minz;
  int ix0,iy0,ix,iy;
  int nn=20;    // number of points per circle
  numero *cth, *sth;

  cth=new numero [nn];  sth=new numero [nn];
  for(j=0; j<nn; j++) {
    cth[j]=cos((2*pi*j)/nn);
    sth[j]=sin((2*pi*j)/nn);
  }

  th*=pi/180;  phi=phi*pi/180;    // phi=phi*pi/180-pi/2;

  // --- init ---

  gdImagePtr plot_graph;
  plot_graph = gdImageCreate(dimx, dimy);  // physical dimensions of graph

  // --- colors ---

  int *colors;  colors=new int [n];

  int white = gdImageColorAllocate(plot_graph, 255, 255, 255);
  int black = gdImageColorAllocate(plot_graph,   0,   0,   0);

  for(i=0; i<n; i++)  colors[i] = plot_color(plot_graph, ci[i]);

  gdImageColorTransparent(plot_graph, white);

  // --- projection ---

  xx=new numero [n];  yy=new numero [n];  zz=new numero [n];

  x0=sin(th)*cos(phi);
  y0=sin(th)*sin(phi);
  z0=cos(th);

  n1x=-sin(phi);
  n1y= cos(phi);
  n1z= 0;

  n2x=-cos(th)*cos(phi);
  n2y=-cos(th)*sin(phi);
  n2z= sin(th);

  for(i=0; i<n; i++) {
    xx[i]=x[i]*n1x+y[i]*n1y+z[i]*n1z;
    yy[i]=x[i]*n2x+y[i]*n2y+z[i]*n2z;
    zz[i]=x[i]*x0 +y[i]*y0 +z[i]*z0;
  }

  // --- limits of x and y ---

  if(minx>maxx || miny>maxy) {

    minx=infinity;  maxx=-infinity;
    miny=infinity;  maxy=-infinity;

    for(i=0; i<n; i++) {
      if(xx[i]-ai[i]<minx) minx=xx[i]-ai[i];
      if(xx[i]+ai[i]>maxx) maxx=xx[i]+ai[i];
      if(yy[i]-ai[i]<miny) miny=yy[i]-ai[i];
      if(yy[i]+ai[i]>maxy) maxy=yy[i]+ai[i];
    }

    val=ABS((maxy-miny)-(maxx-minx))/2;

    if(maxy-miny>maxx-minx) {maxx+=val;  minx-=val;}
    else                    {maxy+=val;  miny-=val;}

    val=0.2*(maxx-minx);  maxx+=val;  minx-=val;
    val=0.2*(maxy-miny);  maxy+=val;  miny-=val;
  }

  // --- resize atoms ---

  if(size) {

    minz=infinity;  maxz=-infinity;
    for(i=0; i<n; i++) {if(zz[i]<minz) minz=zz[i]; if(zz[i]>maxz) maxz=zz[i];}

    for(i=0; i<n; i++) {
      ai[i]=ai[i]*(3+2*(zz[i]-minz)/(maxz-minz))/5;
  } }

  // --- frame ---

  //  gdImageLine(plot_graph,0,     0,     dimx-1,0,     black);
  //  gdImageLine(plot_graph,dimx-1,0,     dimx-1,dimy-1,black);
  //  gdImageLine(plot_graph,dimx-1,dimy-1,0,     dimy-1,black);
  //  gdImageLine(plot_graph,0,     dimy-1,0,     0,     black);

  // --- axes ---

  if(axis) {

    bx[0]=by[0]=bz[0]= infinity;    // --- limits of x, y, and z ---
    bx[1]=by[1]=bz[1]=-infinity;

    for(i=0; i<n; i++) {
      if(x[i]-ai[i]<bx[0]) bx[0]=x[i]-ai[i];
      if(y[i]-ai[i]<by[0]) by[0]=y[i]-ai[i];
      if(z[i]-ai[i]<bz[0]) bz[0]=z[i]-ai[i];
      if(x[i]+ai[i]>bx[1]) bx[1]=x[i]+ai[i];
      if(y[i]+ai[i]>by[1]) by[1]=y[i]+ai[i];
      if(z[i]+ai[i]>bz[1]) bz[1]=z[i]+ai[i];
    }

    d=infinity;                     // --- most distant corner ---
    for(i=0; i<2; i++)
    for(j=0; j<2; j++)
    for(k=0; k<2; k++) {
      dd=x0*bx[i]+y0*by[j]+z0*bz[k];
      if(dd<d) {d=dd; i0=i; j0=j; k0=k;}
    }
    xv=bx[i0]*n1x+by[j0]*n1y+bz[k0]*n1z;
    yv=bx[i0]*n2x+by[j0]*n2y+bz[k0]*n2z;
    ix0=int(floor((dimx-1)*(xv-minx)/(maxx-minx)));
    iy0=int(floor((dimy-1)*(1-(yv-miny)/(maxy-miny))));

    i=1-i0;  j=1-j0;  k=1-k0;

    xv=bx[i]*n1x+by[j0]*n1y+bz[k0]*n1z;       // --- x axis
    yv=bx[i]*n2x+by[j0]*n2y+bz[k0]*n2z;
    ix=int(floor((dimx-1)*(xv-minx)/(maxx-minx)));
    iy=int(floor((dimy-1)*(1-(yv-miny)/(maxy-miny))));
    gdImageDashedLine(plot_graph,ix0,iy0,ix,iy, black);

    xv=bx[i0]*n1x+by[j]*n1y+bz[k0]*n1z;       // --- y axis
    yv=bx[i0]*n2x+by[j]*n2y+bz[k0]*n2z;
    ix=int(floor((dimx-1)*(xv-minx)/(maxx-minx)));
    iy=int(floor((dimy-1)*(1-(yv-miny)/(maxy-miny))));
    gdImageDashedLine(plot_graph,ix0,iy0,ix,iy, black);

    xv=bx[i0]*n1x+by[j0]*n1y+bz[k]*n1z;       // --- z axis
    yv=bx[i0]*n2x+by[j0]*n2y+bz[k]*n2z;
    ix=int(floor((dimx-1)*(xv-minx)/(maxx-minx)));
    iy=int(floor((dimy-1)*(1-(yv-miny)/(maxy-miny))));
    gdImageDashedLine(plot_graph,ix0,iy0,ix,iy, black);
  }

  // --- plotting cluster components ---

  gdPoint *points;  points=new gdPoint [2*nn];

  for(i=0; i<n; i++) {
    val=infinity/10;
    for(j=0; j<n; j++) if(zz[j]<val) {val=zz[j]; k=j;}

    for(j=0; j<nn; j++) {
      xv=xx[k]+ai[k]*cth[j];
      yv=yy[k]+ai[k]*sth[j];
      points[j].x=int(floor((dimx-1)*(xv-minx)/(maxx-minx)));
      points[j].y=int(floor((dimy-1)*(1-(yv-miny)/(maxy-miny))));
    }

    gdImageFilledPolygon(plot_graph, points, nn, colors[k]);
    gdImagePolygon      (plot_graph, points, nn, black);

    zz[k]=infinity;
  }

  delete [] points;

  // --- output ---

  gdImageInterlace(plot_graph, 1);    // interlace when viewing

  FILE *fout;
  fout=fopen(name, "wb");
  gdImageGif(plot_graph, fout);
  fclose(fout);

  // --- ending ---

  gdImageDestroy(plot_graph);
  delete [] colors;  delete [] cth;  delete [] sth;
  delete [] xx;  delete [] yy;  delete [] zz;
}

// **************************************************************************
// *** the same as plot_cluster, but always automat. setting the limits   ***
// **************************************************************************

void plot_cluster(char *name, numero *x, numero *y, numero *z, int *t, int n,
          numero *ai, int *ci, numero th, numero phi, int dimx, int dimy,
          int axis, int size)
{
  numero minx=1,maxx=0,miny=1,maxy=0;
  plot_cluster(name,x,y,z,t,n,ai,ci,th,phi,
               dimx,dimy,minx,maxx,miny,maxy,axis,size);
}

// --------------------------------------------------------------------------

void plot_cluster(char *name, numero *x, numero *y, numero *z, int *t, int n,
          numero *ai, int *ci, numero th, numero phi, int dimx, int dimy)
{
  numero minx=1,maxx=0,miny=1,maxy=0;
  plot_cluster(name, x,y,z,t,n,ai,ci,th,phi, dimx,dimy, 1,0);
}

// **************************************************************************
// *** plot_hologram -> plots angular distributions in gif format         ***
// *** *intensity contains the values to plot                             ***
// **************************************************************************

void normalize_angle(numero &a)     // normalizes angles to be between
{                                   // 0 and 360 degrees
  while(a<0) a+=360;
  while(a>360) a-=360;
}

// --------------------------------------------------------------------------

void plot_hologram(char *name,    // output file
          numero *intensity,      // points in curves
          numero tth0,            // initial angle th actually plotted
          numero tth1,            // final angle th actually plotted
          numero th0,             // initial angle th
          numero th1,             // final angle th          * all angles in
          int    nth,             // number of th values     * degrees
          numero phi0,            // initial angle phi
          numero phi1,            // final angle phi
          int    nphi,            // number of phi values
          int    dim,             // the image will contain dim x dim pixels
          int    ilog)            // 0/1 for linear/logarithmic represent.
{
  int i,j,ij,ith,iphi, *colors,ccc,dimx,dimy;  colors=new int [255];
  numero mmax,mmin,average,absave,val, th,phi,x,y,ii,xth,xphi, l10=log(10);
  numero *intensity0;

  // --- init ---

  if(0<=phi0 && 0<=phi1 && phi0<=90 && phi1<=90)
    dimx=dimy=dim;  else

  if(0<=phi0 && 0<=phi1 && phi0<=180 && phi1<=180)
    {dimx=dim; dimy=dim/2;} else

  if(-90<=phi0 && -90<=phi1 && phi0<=90 && phi1<=90)
    {dimx=dim/2; dimy=dim;} else

  dimx=dimy=dim;

  gdImagePtr plot_graph;
  plot_graph = gdImageCreate(dimx, dimy);  // physical dimensions of graph

  // --- colors ---

  int ncolors=256;

  int white = gdImageColorAllocate(plot_graph, 255, 255, 255);
  for(i=0; i<ncolors; i++)
    colors[i] = gdImageColorAllocate(plot_graph,
                               i,i,i);
//                             51*((i/36)%6), 51*((i/6)%6), 51*(i%6));
  gdImageColorTransparent(plot_graph, white);

  // --- performing a logarithmic transformation ---

  if(ilog) {
    intensity0=intensity;
    intensity=new numero [nth*nphi];
    for(i=0; i<nth; i++)
    for(j=0; j<nphi; j++) {
      val=intensity0[i*nphi+j];
      if(val<1e-30)  val=1e-30;
      intensity[i*nphi+j]=log(val)/l10;
  } }

  // --- normalizing data ... ---

  mmax=-infinity;  mmin=infinity;  average=0;
  for(i=0; i<nth; i++)
  for(j=0; j<nphi; j++) {
    val=intensity[i*nphi+j];
    if(val>mmax)  mmax=val;
    if(val<mmin)  mmin=val;
    average+=val;
  }
  average=average/(nth*nphi);

  absave=0;
  for(i=0; i<nth; i++)
  for(j=0; j<nphi; j++)  absave+=ABS(intensity[i*nphi+j]-average);
  absave=absave/(nth*nphi);

  if(mmax>average+2*absave)  mmax=average+2*absave;
  if(mmin<average-2*absave)  mmin=average-2*absave;

  // --- ... and adjusting the limits of intensity appropriately ---

  for(i=0; i<nth; i++)
  for(j=0; j<nphi; j++)
    intensity[i*nphi+j]=(intensity[i*nphi+j]-mmin)/(mmax-mmin);

  // --- plotting hologram ---

  for(i=0; i<dimx; i++)
  for(j=0; j<dimy; j++) {

    if(0<=phi0 && 0<=phi1 && phi0<=90 && phi1<=90) {
      x= i/(dimx-1.0);
      y=(dimy-j)/(dimy-1.0);
    } else
    if(0<=phi0 && 0<=phi1 && phi0<=180 && phi1<=180) {
      x= (2*i-(dimx-1))/(dimx-1.0);
      y=(dimy-j)/(dimy-1.0);
    } else
    if(-90<=phi0 && -90<=phi1 && phi0<=90 && phi1<=90) {
      x= i/(dimx-1.0);
      y=-(2*j-(dimy-1))/(dimy-1.0);  // j runs from up to bottom in plot
    } else {
      x= (2*i-(dimx-1))/(dimx-1.0);
      y=-(2*j-(dimy-1))/(dimy-1.0);
    }

    if(th1>th0)  th=th1*sqrt(sqr(x)+sqr(y));
    else         th=th0*sqrt(sqr(x)+sqr(y));

    if(th>=th0 && th<=th1 && th>=tth0 && th<=tth1) {
      phi=180/pi*atan2(y,x); // -(phi1-phi0)/nphi;
      if(phi0>=0 && phi1>=0)  normalize_angle(phi);
      if(phi0<=phi && phi<=phi1) {
        xth=(nth-1)*(th-th0)/(th1-th0);         ith=int(floor(xth));
        xphi=(nphi-1)*(phi-phi0)/(phi1-phi0);   iphi=int(floor(xphi));
        if(ith<0) ith=0;
        if(iphi<0) iphi=0;
        ij=ith*nphi+iphi;
        ii=intensity[ij] + (xth-ith)*(intensity[ij+nphi]-intensity[ij])
                         + (xphi-iphi)*(intensity[ij+1]-intensity[ij]);
        if(ii>1) ii=1;   if(ii<0) ii=0;   ii=(ncolors-1)*ii;
        ccc=int(floor(ii));
        if(ccc<1) ccc=1;
        gdImageSetPixel(plot_graph, i,j, ccc);
  } } }

  // --- output ---

  gdImageInterlace(plot_graph, 1);    // interlace when viewing

  FILE *fout;
  fout=fopen(name, "wb");
  gdImageGif(plot_graph, fout);
  fclose(fout);

  // --- ending ---

  gdImageDestroy(plot_graph);
  delete [] colors;

  if(ilog) delete [] intensity;
}

// --------------------------------------------------------------------------

void plot_hologram(char *name, numero *intensity,
          numero tth0, numero tth1, numero th0, numero th1, int nth,
          numero phi0, numero phi1, int nphi, int dim)
{
  plot_hologram(name, intensity, tth0,tth1,
                th0,th1,nth, phi0,phi1,nphi, dim, 0);
}

// **************************************************************************
// *** plot_color_hemi -> plots angular distributions of RGB colors       ***
// *** *iR, *iG, *iB contain the RGB parameters for each (th,phi)         ***
// **************************************************************************

//(nameo,iR,iG,iB, th0,th1,phi0,phi1, nth,nphi, dim,intensity);

void plot_color_hemi(char *name,  // output file
          numero *iR,             // Red   intensity
          numero *iG,             // Green intensity
          numero *iB,             // Blue  intensity
          numero th0,             // initial angle th
          numero th1,             // final angle th   (all angles in degrees)
          numero phi0,            // initial angle phi      (on input)
          numero phi1,            // final angle phi
          int    nth,             // number of th values
          int    nphi,            // number of phi values
          int    dim,             // the image will contain dim x dim pixels
          numero II)              // multiplicative factor for RGB param.
{
  int i,j,jj,ij,ith,iphi,ccc,RR,GG,BB,*colors,ll,l;  colors=new int [dim*dim];
  numero th,phi,x,y,xth,xphi,R,G,B;

  // --- init ---

  gdImagePtr plot_graph;
  plot_graph = gdImageCreate(dim,dim);
  int white=gdImageColorAllocate(plot_graph, 255, 255, 255);

  // --- plotting hologram ---

  for(l=0; l<255; l++) {
    xth=l*(nth-1)/(255-1.0);  ith=int(floor(xth));
    if(ith<0) ith=0;
    if(ith>nth-2) ith=nth-2;
    R=iR[ith]+(xth-ith)*(iR[ith+1]-iR[ith]);
    G=iG[ith]+(xth-ith)*(iG[ith+1]-iG[ith]);
    B=iB[ith]+(xth-ith)*(iB[ith+1]-iB[ith]);
    R=R*II;  G=G*II;  B=B*II;
    if(R>1) R=1;   if(R<0) R=0;   RR=int(255*R);
    if(G>1) G=1;   if(G<0) G=0;   GG=int(255*G);
    if(B>1) B=1;   if(B<0) B=0;   BB=int(255*B);
    colors[ith]=gdImageColorAllocate(plot_graph, RR,GG,BB);
  }

  for(l=1; l<2; l++)
  for(i=ll=0; i<dim; i++)
  for(jj=0; jj<dim; jj++,ll++) {  j=dim-jj-1;

    x= (2*i-(dim-1))/(dim-1.0);
    y=-(2*j-(dim-1))/(dim-1.0);

    th=90*sqrt(sqr(x)+sqr(y));

    if(th0<=th && th<=th1) {
      phi=180/pi*atan2(y,x);
      if(phi0>=0 && phi1>=0)  normalize_angle(phi);
      if(phi0<=phi && phi<=phi1) {
        xth=(nth-1)*(th-th0)/(th1-th0);         ith=int(floor(xth));
        xphi=(nphi-1)*(phi-phi0)/(phi1-phi0);   iphi=int(floor(xphi));
        if(ith<0) ith=0;
        if(iphi<0) iphi=0;
        if(nphi==0) {
          R=iR[ith]+(xth-ith)*(iR[ith+1]-iR[ith]);
          G=iG[ith]+(xth-ith)*(iG[ith+1]-iG[ith]);
          B=iB[ith]+(xth-ith)*(iB[ith+1]-iB[ith]);
        } else {
          ij=ith*nphi+iphi;
          R=iR[ij]+(xth-ith)*(iR[ij+nphi]-iR[ij])+(xphi-iphi)*(iR[ij+1]-iR[ij]);
          G=iG[ij]+(xth-ith)*(iG[ij+nphi]-iG[ij])+(xphi-iphi)*(iG[ij+1]-iG[ij]);
          B=iB[ij]+(xth-ith)*(iB[ij+nphi]-iB[ij])+(xphi-iphi)*(iB[ij+1]-iB[ij]);
        }
        R=R*II;  G=G*II;  B=B*II;
        if(R>1) R=1;   if(R<0) R=0;   RR=int(255*R);
        if(G>1) G=1;   if(G<0) G=0;   GG=int(255*G);
        if(B>1) B=1;   if(B<0) B=0;   BB=int(255*B);
        if(nphi==0) ith=int((255-1.0)*xth/(nth-1));
	//        if(l==0) colors[ll]=gdImageColorAllocate(plot_graph, RR,GG,BB);
        if(l==1) gdImageSetPixel(plot_graph, i,j,ith+1);
  } } }

  // --- output ---

  gdImageInterlace(plot_graph, 1);    // interlace when viewing

  FILE *fout;
  fout=fopen(name, "wb");
  gdImageGif(plot_graph, fout);
  fclose(fout);

  // --- ending ---

  gdImageDestroy(plot_graph);
  delete [] colors;
}

// **************************************************************************
// *** plot_density -> plots 2D rectangular density in gif format         ***
// *** *intensity contains the values to plot                             ***
// **************************************************************************

void plot_density(
          char *name,             // output file
          numero *intensity,      // points in curves
          numero x0, numero x1,   // the plotted area is contained
          numero y0, numero y1,   // in the rectangle (x0,y0)<->(x1,y1);
          numero xx0,numero xx1,  // the data in intensity are defined in
          numero yy0,numero yy1,  // (xx0,yy0)<->(xx1,yy1);
          int    nx,              // number of x's in input
          int    ny,              // number of y's in input
          int    dimx,            // the image will contain dimx x dimy
          int    dimy,            // pixels
          numero IImin,           // black is plotted below this limit
          numero IImax,           // white is plotted above this limit
          int    limits,          // set to 0/1/2 to saturate, b&w, w&b
          numero *xcurve,         // if ncurve>0, a curve is plotted ...
          numero *ycurve,         // ... on top of the density plot
          int    *ccurve,         // ... with ncurve (x,y) points of color c
          int    ncurve,
          int    IInorm,          // set to 1 to normalize data to IImax...
          int    bw,              // set to 1 to plot in black and white
          char   *dtext)          // text to include in the plot
{
  int i,ic,j,jj, ix,iy, *colors, ccc, im,id, ir,ig,ib;  colors=new int [256];
  numero mmax,mmin,val, xx,xy, x,y,ii;

  // --- init ---

  gdImagePtr plot_graph;
  plot_graph = gdImageCreate(dimx, dimy);  // physical dimen. of graph
  int white=gdImageColorAllocate(plot_graph, 255,255,255);
  int black=gdImageColorAllocate(plot_graph, 0,0,0);

  // --- colors ---

  int ncolors=256-4, colcurve;

  if(ncurve>0) {
    if(ccurve[0]==0) colcurve=black; else
    if(ccurve[0]==1) colcurve=gdImageColorAllocate(plot_graph,255,0,  0); else
    if(ccurve[0]==2) colcurve=gdImageColorAllocate(plot_graph,  0,0,255); else
    if(ccurve[0]==3) colcurve=white;
  }

  for(i=0; i<ncolors; i++) {
    ic=ncolors-1-i;
    if(bw==0) {
      im=ic%(ncolors/4);  id=ic/(ncolors/4);
      if(id==0) {ir=255;      ig=4*im;     ib=0;}    else
      if(id==1) {ir=255-4*im; ig=255;      ib=0;}    else
      if(id==2) {ir=0;        ig=255;      ib=4*im;} else
                {ir=0;        ig=255-4*im; ib=255;}
      colors[i]=gdImageColorAllocate(plot_graph,ir,ig,ib);
    }
    else if(bw==1) colors[i]=gdImageColorAllocate(plot_graph, i-1,i-1,i-1);
    else if(bw==2) colors[i]=gdImageColorAllocate(plot_graph, i-1,0,0);
    else if(bw==3) colors[i]=gdImageColorAllocate(plot_graph, 0,0,i-1);
    else if(bw==4) colors[i]=gdImageColorAllocate(plot_graph, ic,ic,ic);
    else if(bw==5) colors[i]=gdImageColorAllocate(plot_graph, 255,ic,ic);
    else if(bw==6) colors[i]=gdImageColorAllocate(plot_graph, ic,ic,255);
    else if(bw==14) colors[i]=gdImageColorAllocate(plot_graph, i-1,i-1,0);
    else if(bw==25) colors[i]=gdImageColorAllocate(plot_graph, ic,255,255);
    else if(bw==26) colors[i]=gdImageColorAllocate(plot_graph,
         127+ic/2,127+ic/2,127+ic/2);
    else if(bw==27) colors[i]=gdImageColorAllocate(plot_graph, ic,255,255);
    else if(bw==28) colors[i]=gdImageColorAllocate(plot_graph, 255,ic,255);
    else if(bw==29) colors[i]=gdImageColorAllocate(plot_graph, 255,255,ic);
    else if(bw==30) colors[i]=gdImageColorAllocate(plot_graph, ic,ic,255);
    else if(bw==31) colors[i]=gdImageColorAllocate(plot_graph, ic,255,ic);
    else if(bw==32) colors[i]=gdImageColorAllocate(plot_graph, 255,ic,ic);
    else {
      if(bw==7) {
        if(i<128) {ir=ig=0;  ib=2*(127-i);}
        else {im=i-128;
          if(im<64) {ir=4*im; ig=0; ib=0;}
          else {ir=252; ig=4*(im-63); ib=4*(im-63);}
      } }
      if(bw==50) {
        if(i<128) {ir=ig=0;  ib=2*i;}
        else {im=i-128;
          if(im<64) {ir=4*im+3; ig=0; ib=256-4*im;}
          else {ir=255; ig=ib=int(4*(im-64)*255/248.0);}
      } }
      if(bw==60) {
        if(i>=10 && i<=110) {ir=int((110-i)*2.55); ig=ib=0;} else
        if(i>=111 && i<=211) {ib=255; ir=ig=int((211-i)*2.55);} else
        {ir=ig=255; ib=0;}
      }
      else if(bw==8) {  // blue-black-red
        if(i<128) {ir=ig=0;  ib=2*(127-i);}
        else      {ig=ib=0;  ir=2*(i-128);}
      }
      else if(bw==9) {  // blue-white-red
        if(i<128) {ib=255;  ir=ig=255-2*(127-i);}
        else      {ir=255;  ig=ib=255-2*(i-128);}
      }
      else if(bw==10) {  // yellow-black-red
        if(i<128) {ib=0;  ir=ig=2*(127-i);}
        else      {ig=ib=0;  ir=2*(i-128);}
      }
      else if(bw==11) {  // yellow-white-red
        if(i<128) {ir=ig=255;  ib=255-2*(127-i);}
        else      {ir=255;     ig=ib=255-2*(i-128);}
      }
      else if(bw==12) {  // green-white-red
        if(i<128) {ig=255;  ir=ib=255-2*(127-i);}
        else      {ir=255;  ig=ib=255-2*(i-128);}
      }
      else if(bw==13) {  // green-black-red
        if(i<128) {ir=ib=0;  ig=2*(127-i);}
        else      {ig=ib=0;  ir=2*(i-128);}
      }
      else if(bw==15) {  // black-red-yellow
        if(i<128) {ib=ig=0;       ir=2*i;}
        else      {ir=255; ib=0;  ig=2*(i-128);}
      }
      else if(bw==16) {  // black-yellow-white
        if(i<128) {ib=0;          ir=ig=2*i;}
        else      {ir=ig=255;     ib=2*(i-128);}
      }
      else if(bw==20) {  // black-red-white
        if(i<128) {ig=ib=0;  ir=2*(127-i);}
        else      {ir=255;  ig=ib=2*(i-128);}
      }
      else if(bw==21) {  // black-green-white
        if(i<128) {ir=ib=0;  ig=2*(127-i);}
        else      {ig=255;  ir=ib=2*(i-128);}
      }
      else if(bw==22) {  // black-blue-white
        if(i<128) {ir=ig=0;  ib=2*(127-i);}
        else      {ib=255;  ir=ig=2*(i-128);}
      }
      else if(bw==23) {  // black-pink-white
        if(i<128) {ig=0;          ir=ib=2*i;}
        else      {ir=ib=255;     ig=2*(i-128);}
      }
      else if(bw==24) {  // black-cyan-white
        if(i<128) {ir=0;          ig=ib*i;}
        else      {ig=ib=255;     ir=2*(i-128);}
      }
      colors[i]=gdImageColorAllocate(plot_graph,ir,ig,ib);
  } }

  colors[255]=gdImageColorAllocate(plot_graph, 255,255,255);

  // --- normalizing data ---

  if(IInorm) {mmin=IImin;  mmax=IImax;}
  else {
    mmax=-infinity;  mmin=infinity;
    for(i=0; i<nx; i++)
    for(j=0; j<ny; j++) {
      x=xx0+i*(xx1-xx0)/(nx-1);
      y=yy0+j*(yy1-yy0)/(ny-1);
      if(x0<=x && x<=x1 && y0<=y && y<=y1) {
        val=intensity[i*ny+j];
        if(val>mmax && val<=IImax)  mmax=val;
        if(val<mmin && val>=IImin)  mmin=val;
  } } }

  numero *iintensity;  iintensity=new numero [nx*ny];

  for(i=0; i<nx; i++)
  for(j=0; j<ny; j++)
    if(bw==60)  iintensity[i*ny+j]=intensity[i*ny+(ny-1-j)]; else
    iintensity[i*ny+j]=(intensity[i*ny+(ny-1-j)]-mmin)/(mmax-mmin);

  // --- plotting density ---

  for(i=0; i<dimx; i++)
  for(jj=0; jj<dimy; jj++) {  j=dimy-jj-1;

    x=x0+(x1-x0)*i/(dimx-1.0);
    y=y0+(y1-y0)*j/(dimy-1.0);

    if(xx0<=x && x<=xx1 && yy0<=y && y<=yy1) {
      xx=(nx-1)*(x-xx0)/(xx1-xx0);   ix=int(floor(xx));
      xy=(ny-1)*(y-yy0)/(yy1-yy0);   iy=int(floor(xy));
      ii=iintensity[ix*ny+iy];
      if(ix<nx-1) ii+=(xx-ix)*(iintensity[(ix+1)*ny+iy]-iintensity[ix*ny+iy]);
      if(iy<ny-1) ii+=(xy-iy)*(iintensity[ix*ny+iy+1]-iintensity[ix*ny+iy]);
      // ii=iintensity[ix*ny+iy]
      //   +(xx-ix)*(iintensity[(ix+1)*ny+iy]-iintensity[ix*ny+iy])
      //   +(xy-iy)*(iintensity[ix*ny+iy+1]-iintensity[ix*ny+iy]);
      if(bw==60) {
        if(-1<=ii && ii<=0) ccc=int((ii+1)*100+10); else
        if(0<=ii && ii<=1) ccc=int(ii*100+111); else
        ccc=220;
      } else {
        if(ii>=1) {if(limits==1 || limits==3) ccc=white; if(limits==2 || limits==5) ccc=black; else ccc=colors[ncolors-1];} else
        if(ii<=1e-12) {if(limits==1 || limits==4) ccc=black; if(limits==2 || limits==6) ccc=white; else ccc=colors[0];}   else {
          ccc=colors[int(floor((ncolors-1)*ii))];
      } }
      gdImageSetPixel(plot_graph, i,j, ccc);
  } }

  for(jj=0; jj<ncurve; jj++) {
    i=int(floor((xcurve[jj]-x0)/(x1-x0)*(dimx-1)));
    j=int(floor((ycurve[jj]-y0)/(y1-y0)*(dimy-1)));
    if(i>=0 && i<dimx && j>=0 && j<dimy)
      gdImageSetPixel(plot_graph, i,dimy-1-j, colcurve);
  }

  if(dtext!=NULL)
    gdImageString(plot_graph, gdFontMediumBold,
                       30,
                       dimy-30,
                       (unsigned  char *) dtext, black);

  // --- output ---

  //  gdImageInterlace(plot_graph, 1);    // interlace when viewing

  FILE *fout;
  fout=fopen(name, "wb");
  gdImageGif(plot_graph, fout);
  fclose(fout);

  // --- ending ---

  gdImageDestroy(plot_graph);
  delete [] colors;  delete [] iintensity;
}

// --------------------------------------------------------------------------

void plot_density(
          char *name,             // output file
          numero *intensity,      // points in curves
          numero x0, numero x1,   // the area plotted goes is contained
          numero y0, numero y1,   // in the rectangle (x0,y0)<->(x1,y1);
          numero xx0,numero xx1,  // the data in intensity are defined in
          numero yy0,numero yy1,  // (xx0,yy0)<->(xx1,yy1);
          int    nx,              // number of x's in input
          int    ny,              // number of y's in input
          int    dimx,            // the image will contain dimx x dimy
          int    dimy,            // pixels
          numero IImin,           // black is plotted below this limit
          numero IImax,           // white is plotted above this limit
          int    limits,   // set to 0 to plot blue/red instead of black/red
          numero *xcurve,         // if ncurve>0, a curve is plotted ...
          numero *ycurve,         // ... on top of the density plot
          int    *ccurve,         // ... with ncurve (x,y) points of color c
          int    ncurve,
          int    IInorm)          // set to 1 to normalize data to IImax...
{
  plot_density(name,intensity,x0,x1,y0,y1,xx0,xx1,yy0,yy1,
               nx,ny,dimx,dimy,IImin,IImax,limits,xcurve,ycurve,ccurve,
               ncurve,IInorm,0,NULL);
}

// --------------------------------------------------------------------------

void plot_density(
          char *name,             // output file
          numero *intensity,      // points in curves
          numero x0, numero x1,   // the area plotted goes is contained
          numero y0, numero y1,   // in the rectangle (x0,y0)<->(x1,y1);
          numero xx0,numero xx1,  // the data in intensity are defined in
          numero yy0,numero yy1,  // (xx0,yy0)<->(xx1,yy1);
          int    nx,              // number of x's in input
          int    ny,              // number of y's in input
          int    dimx,            // the image will contain dimx x dimy
          int    dimy)            // pixels
{
  plot_density(name,intensity,x0,x1,y0,y1,xx0,xx1,yy0,yy1,
               nx,ny,dimx,dimy,-infinity,infinity,1,NULL,NULL,NULL,0,0);
}

// --------------------------------------------------------------------------

void plot_density(char *name, numero *intensity,
          numero x0, numero x1, numero y0, numero y1,
          int nx, int ny, int dimx, int dimy)
{
  plot_density(name, intensity, x0,x1,y0,y1, x0,x1,y0,y1, nx,ny, dimx,dimy);
}

// --------------------------------------------------------------------------

void plot_density(char *name, numero *intensity,
          numero x0, numero x1, numero y0, numero y1, int n, int dim)
{
  plot_density(name, intensity, x0,x1,y0,y1, x0,x1,y0,y1, n,n, dim,dim);
}

#endif  // ******************************************************************
