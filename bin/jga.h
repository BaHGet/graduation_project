#ifndef jga_jga           // ************************************************
#define jga_jga 1         // ***  jga/jga.h                               ***
                          // ***                                          ***
#include <stdio.h>        // ***  <-- standard libraries of common use    ***
#include <stdlib.h>       // ***                                          ***
#include <string.h>       // ***                                          ***
#include <time.h>         // ***                                          ***
#include <math.h>         // ***                                          ***
                          // ***                                          ***
#define numero double     // ***  <-- definition of 'numero'              ***
                          // ***                                          ***
// ***  ----------------  // ***  --------------------------------------  ***
// ***                                                                    ***
// ***  Author: F. J. Garcia de Abajo, CSIC, Spain                        ***
// ***  Routines of general use in the jga library                        ***
// ***                                                                    ***
// ***  ----------------------------------------- CONSTANTS ------------  ***
//                                                                        ***
        numero                                        //                  ***
                                                      //                  ***
           pi      = 3.141592653589793238462643,      //                  ***
                                                      //                  ***
           au_eV   = 27.2113834,         // 1 a.u. of enery in eV         ***
           a0_au   = 0.529177208,        // Bohr radius a0 in Angstroms   ***
           nm      = 10/a0_au,           // 1 nm in a.u.                  ***
           c_au    = 137.03599976;       // 1/alfa=speed of light (a.u.)  ***
                                         //                               ***
           #define infinite      1.0e20  // effective infinity            ***
           #define infinity      1.0e20  // effective infinity            ***
           #define infinite_int  10000   // effective infinite integer    ***
           #define infinitesimal 1.0e-20 // effective epsilon             ***
//                                                                        ***
// ***  ----------------------------------------- ROUTINES -------------  ***
// ***                                                                    ***
// ***  --- elementary function                                           ***
// ***                                                                    ***
// ***         sqr(x);                       x^2 for int's and numeros    ***
// ***         ABS(x);                       |x| for int's and numeros    ***
// ***         sign_1l(l)                    (-1)^l                       ***
// ***         SIGN(x)                     x<0 -> -1, x>0 -> 1, x=0 -> 0  ***
// ***                                                                    ***
// ***  --- strcmpC(str1,str2)                                            ***
// ***                                                                    ***
// ***      like strcmp(...) with additional features                     ***
// ***                                                                    ***
// ***  --- input strings, int's, and numeros disregarding comments       ***
// ***      with the same format as comments in C++; also, 'infinity'     ***
// ***      and 'zero' are understood as a large number and 0, respect.   ***
// ***      (see routine also for special role of 'all')                  ***
// ***                                                                    ***
// ***         read_name(f,name);            read 'name' from file f      ***
// ***         int read_int(f);              read int from file f         ***
// ***         numero read_numero(f);        read numero from file f      ***
// ***         int read_int(name);           read int from string name    ***
// ***         numero read_numero(name);     read numero from str. name   ***
// ***         read_name(name);              read 'name' from stdin       ***
// ***         int read_int();               read int from stdin          ***
// ***         numero read_numero();         read numero from stdin       ***
// ***                                                                    ***
// ***  --- char *ftostr(x,d);                                            ***
// ***                                                                    ***
// ***      real x -> string with d figures after the decimal point       ***
// ***                                                                    ***
// ***  --- error messages; on_error abort the program after              ***
// ***      displaying an error message; their format is                  ***
// ***                                                                    ***
// ***         on_error  (where, comment);                                ***
// ***         on_error  (where, comment1, comment2);                     ***
// ***         on_error  (where, comment1, x, comment2);                  ***
// ***         on_error  (fout, where, comment);                          ***
// ***         on_error  (fout, where, comment1, comment2);               ***
// ***         on_error  (fout, where, comment1, x, comment2);            ***
// ***                                                                    ***
// ***      and the error message displayed is                            ***
// ***                                                                    ***
// ***         *** error in where: comment                                ***
// ***         *** error in where: comment1 comment2                      ***
// ***         *** error in where: comment1 x comment2                    ***
// ***                                                                    ***
// **************************************************************************

// **************************************************************************
// *** miscellaneous                                                      ***
// **************************************************************************

FILE *foutput=stdout;        // file(s) normally used as the standard output

// **************************************************************************
// *** elementary functions for int's and numeros                         ***
// **************************************************************************

numero sqr(numero a) {return a*a;}                        // a^2
int    sqr(int    i) {return i*i;}                        // i^2

// --------------------------------------------------------------------------

numero ABS(numero a) {return (a<0)?-a:a;}                 // |a|
int    ABS(int    i) {return (i<0)?-i:i;}                 // |i|

// --------------------------------------------------------------------------

numero SIGN(numero a) {return (a<0)?-1:((a>0)?1:0);}      // sign(a)
int    SIGN(int    i) {return (i<0)?-1:((i>0)?1:0);}      // sign(i)

// --------------------------------------------------------------------------

int sign_1l(int l)  {if(l%2)  return -1;  return 1;}      // (-1)^l

// **************************************************************************
// *** comparare two strings with independence of the case                ***
// **************************************************************************

int strcmpC(char *str1, char *str2)
{
  return strcmp(str1,str2);           // implement!!!  (jga)
}

// **************************************************************************
// *** read input string/int/float and disregard comments                 ***
// **************************************************************************

char read_command[80];                           // internal variable

// --------------------------------------------------------------------------

void read_name(FILE *f, char *command)
{
  fscanf(f,"%s",command);  char c;  int i;       // comments start with /*
                                                 // and end with '*/';
  if(strlen(command)>1)                          // also, '//' like in C++
  if(command[0]=='/' && command[1]=='*') {
    c=1;
    do {
      fscanf(f,"%s",command);  i=strlen(command);
      if(i>1)  if(command[i-2]=='*')  if(command[i-1]=='/')  c=0;
    } while(c);
    read_name(f, command);
  }  else
  if(command[0]=='/' && command[1]=='/') {
    do fscanf(f,"%c",&c); while(c!='\n');
    read_name(f,command);
} }

// --------------------------------------------------------------------------

numero read_numero(char *command)
{
  if(!strcmpC(command,"infinite"))   return infinity;
  if(!strcmpC(command,"infinity"))   return infinity;
  if(!strcmpC(command,"zero"))       return 0;
                                     return atof(command);
}

// --------------------------------------------------------------------------

int read_int(char *command)
{
  if(!strcmpC(command,"infinite"))   return 3*infinite_int;
  if(!strcmpC(command,"infinity"))   return 3*infinite_int;
  if(!strcmpC(command,"all"))        return 2*infinite_int;
  if(!strcmpC(command,"zero"))       return 0;
  if(!strcmpC(command,"off"))        return 0;
  if(!strcmpC(command,"on"))         return 1;
                                     return atoi(command);
}

// --------------------------------------------------------------------------

numero read_numero(FILE *f)
{
  read_name(f,read_command);
  return read_numero(read_command);
}

// --------------------------------------------------------------------------

int read_int(FILE *f)
{
  read_name(f,read_command);
  return read_int(read_command);
}

// --------------------------------------------------------------------------

numero read_numero(void)  {return read_numero(stdout);}
int    read_int(void)     {return read_int(stdout);}

// **************************************************************************
// *** transform real number x into string with d figures after the point ***
// **************************************************************************

char ftostr_data[50];

char *ftostr(double x, int d)
{
  if(d>0)  if(x-floor(x)==0)  return ftostr(x,0);

  char *val,*v;  val=ftostr_data;  v=val;
  double dec=1e20, eps=1e-12;
  int i=0, j;

  if(x<0) {val[i]='-';  i++;  x=ABS(x);}

  while(floor(x/dec+eps)==0 && dec>1)  dec=dec/10;

  while(dec>=1) {
    j=((int) floor(x/dec+eps));
    x=x-j*dec;  val[i]='0'+j;  i++;  dec=dec/10;
  }

  if(ABS(x)>eps) {
    val[i]='.';  i++;
    while(d) {
      j=int(floor(x/dec+eps));   x=x-j*dec;
      if(d==1 && j<9) if(x/dec>=0.5) j++;
      val[i]='0'+j;  i++;  dec=dec/10;
      d--;
  } }

  val[i]='\0';

  return v;
}

// **************************************************************************
// *** error messages                                                     ***
// **************************************************************************

void on_error(FILE *fout, char *where, char *com1, numero x, char *com2)
  {fprintf(fout,"*** error in '%s': %s %g %s\n",where,com1,x,com2); exit(1);}

void on_error(FILE *fout, char *where, char *com1, char *com2)
  {fprintf(fout,"*** error in '%s': %s %s\n",where,com1,com2); exit(1);}

void on_error(FILE *fout, char *where, char *com)
  {fprintf(fout,"*** error in '%s': %s\n",where,com); exit(1);}

void on_error(char *where, char *com1, numero x, char *com2)
  {on_error(foutput, where,com1,x,com2);}

void on_error(char *where, char *com1, char *com2)
  {on_error(foutput, where,com1,com2);}

void on_error(char *where, char *com)  {on_error(foutput, where,com);}

#endif  // ******************************************************************
