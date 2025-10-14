#ifndef jga_algebra       // ************************************************
#define jga_algebra 1     // ***  jga/algebra.h                           ***
                          // ***                                          ***
#include "jga.h"          // ***                                          ***
                          // ***                                          ***
// ***  ----------------  // ***  --------------------------------------  ***
// ***                                                                    ***
// ***  Author: F. J. Garcia de Abajo, CSIC, Spain         21-iii-2006    ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Handling simple algebraic expressions to input data.              ***
// ***                                                                    ***
// **************************************************************************

class algebra_variable {
 public:
  algebra_variable *next;
  char             *name;
  numero           a;
  algebra_variable(void) {name=NULL;}
} *algebra_var, *algebra_var0=NULL;

void algebra_free(void)
{
  while(algebra_var0!=NULL) {
    if(algebra_var0->name!=NULL)  delete [] algebra_var0->name;
    algebra_var=algebra_var0->next;  delete algebra_var0;
    algebra_var0=algebra_var;
} }

int algebra_find(char *name)
{
  algebra_var=algebra_var0;

  if(algebra_var==NULL)  return 0;
  while(strcmp(algebra_var->name,name)) {
    if(algebra_var->next==NULL)  return 0;
    algebra_var=algebra_var->next;
  }
  return 1;
}

int algebra_error_flag=1;

#define algebra_error                                                  \
  if(algebra_error_flag)                                               \
    on_error("algebraic evaluation", "wrong expression: ", ex);   else \
    on_warning("algebraic evaluation", "wrong expression: ", ex);

int algebra_locate(char *ex, int i0, int i1, char c)
{
  int i=i1+1, p=0;  // This returns the position i where character c is
  int cond;         // found in string ex between i0 and i1, or -1 otherwise.
  do {
    i--;
    if(ex[i]=='(' || ex[i]=='[' || ex[i]=='{') p--;
    if(ex[i]==')' || ex[i]==']' || ex[i]=='}') p++;
    if(p<0) algebra_error
    cond=(ex[i]==c);
    if(cond) if(c=='+' || c=='-') {
      if(i>i0) if(ex[i-1]=='e' || ex[i-1]=='d')
      if(i>i0+1) if(ex[i-2]=='.' || ('0'<=ex[i-2] && ex[i-2]<='9')) cond=0;
    }
  } while((!cond || p!=0)  && i>i0);
  if(p!=0) algebra_error
  if(!cond) return -1;
  return i;
}

int algebra_locate(char *ex, char c)
  {return algebra_locate(ex,0,strlen(ex)-1,c);}

int algebra_locate_absolute(char *ex, int i0, int i1, char c)
{
  int i=i0-1;  do i++; while(ex[i]!=c && i<i1);
  if(ex[i]!=c) return -1;
  return i;
}

numero algebra_eval(char *ex, int i0, int i1)
{
  char c,d;  int i;  if(i0>i1) return 0;  c=ex[i0];  d=ex[i1];
  numero val=0;

  if(c=='*' || c==':' || c=='/' || c=='^' || c==')' || c==']' || c=='}' ||
     d=='*' || d==':' || d=='/' || d=='^' || d=='(' || d=='[' || d=='{')
    algebra_error

  if((i=algebra_locate(ex,i0,i1,'+'))>=0)
    return algebra_eval(ex,i0,i-1)+algebra_eval(ex,i+1,i1);

  if((i=algebra_locate(ex,i0,i1,'-'))>=0)
    return algebra_eval(ex,i0,i-1)-algebra_eval(ex,i+1,i1);

  if((i=algebra_locate(ex,i0,i1,'*'))>=0)
    return algebra_eval(ex,i0,i-1)*algebra_eval(ex,i+1,i1);

  if((i=algebra_locate(ex,i0,i1,':'))>=0)
    return algebra_eval(ex,i0,i-1)*algebra_eval(ex,i+1,i1);

  if((i=algebra_locate(ex,i0,i1,'/'))>=0)
    return algebra_eval(ex,i0,i-1)/algebra_eval(ex,i+1,i1);

  if((i=algebra_locate(ex,i0,i1,'^'))>=0)
    return pow(algebra_eval(ex,i0,i-1),algebra_eval(ex,i+1,i1));

  if(c=='(' && d==')')  return algebra_eval(ex,i0+1,i1-1);
  if(c=='[' && d==']')  return algebra_eval(ex,i0+1,i1-1);
  if(c=='{' && d=='}')  return algebra_eval(ex,i0+1,i1-1);

  char *aux;  aux=new char [i1-i0+2];
  for(i=i0; i<=i1; i++) aux[i-i0]=ex[i];  aux[i1-i0+1]=0;  i1=i1-i0;  i0=0;

  if(('a'<=c && c<='z') || ('A'<=c && c<='Z')) {d=aux[i1];
    if(d==')' || d==']' || d=='}') {
      if((i=algebra_locate_absolute(aux,i0,i1,'('))<0)
      if((i=algebra_locate_absolute(aux,i0,i1,'['))<0)
      if((i=algebra_locate_absolute(aux,i0,i1,'{'))<0)
        algebra_error
      val=algebra_eval(aux,i+1,i1-1);
      aux[i]=0;
      #ifdef jga_bessel
      if(!strcmp(aux,"J0"))    val=besselJJ(0,val);  else
      if(!strcmp(aux,"J1"))    val=besselJJ(1,val);  else
      if(!strcmp(aux,"J2"))    val=besselJJ(2,val);  else
      if(!strcmp(aux,"J3"))    val=besselJJ(3,val);  else
      if(!strcmp(aux,"Y0"))    val=besselYY(0,val);  else
      if(!strcmp(aux,"Y1"))    val=besselYY(1,val);  else
      if(!strcmp(aux,"Y2"))    val=besselYY(2,val);  else
      if(!strcmp(aux,"Y3"))    val=besselYY(3,val);  else
      if(!strcmp(aux,"K0"))    val=besselK(0,val);   else
      if(!strcmp(aux,"K1"))    val=besselK(1,val);   else
      if(!strcmp(aux,"K2"))    val=besselK(2,val);   else
      if(!strcmp(aux,"K3"))    val=besselK(3,val);   else
      if(!strcmp(aux,"I0"))    val=besselI(0,val);   else
      if(!strcmp(aux,"I1"))    val=besselI(1,val);   else
      if(!strcmp(aux,"I2"))    val=besselI(2,val);   else
      if(!strcmp(aux,"I3"))    val=besselI(3,val);   else
      if(!strcmp(aux,"j0"))    val=besselj(0,val);   else
      if(!strcmp(aux,"j1"))    val=besselj(1,val);   else
      if(!strcmp(aux,"j2"))    val=besselj(2,val);   else
      if(!strcmp(aux,"j3"))    val=besselj(3,val);   else
      if(!strcmp(aux,"y0"))    val=bessely(0,val);   else
      if(!strcmp(aux,"y1"))    val=bessely(1,val);   else
      if(!strcmp(aux,"y2"))    val=bessely(2,val);   else
      if(!strcmp(aux,"y3"))    val=bessely(3,val);   else
      #endif
      #ifdef jga_legendre
      if(!strcmp(aux,"P00"))   val=legendre(0,0,val);  else
      if(!strcmp(aux,"P10"))   val=legendre(1,0,val);  else
      if(!strcmp(aux,"P11"))   val=legendre(1,1,val);  else
      if(!strcmp(aux,"P20"))   val=legendre(2,0,val);  else
      if(!strcmp(aux,"P21"))   val=legendre(2,1,val);  else
      if(!strcmp(aux,"P22"))   val=legendre(2,2,val);  else
      #endif
      if(!strcmp(aux,"ln"))    val=log(val);    else
      if(!strcmp(aux,"log"))   val=log(val);    else
      if(!strcmp(aux,"exp"))   val=exp(val);    else
      if(!strcmp(aux,"sin"))   val=sin(val);    else
      if(!strcmp(aux,"cos"))   val=cos(val);    else
      if(!strcmp(aux,"tan"))   val=tan(val);    else
      if(!strcmp(aux,"asin"))  val=asin(val);   else
      if(!strcmp(aux,"acos"))  val=acos(val);   else
      if(!strcmp(aux,"atan"))  val=atan(val);   else
      if(!strcmp(aux,"abs"))   val=ABS(val);    else
      if(!strcmp(aux,"mod"))   val=ABS(val);    else
      if(!strcmp(aux,"ABS"))   val=ABS(val);    else
      if(!strcmp(aux,"sqr"))   val=sqr(val);    else
      if(!strcmp(aux,"sqrt"))  val=sqrt(val);   else
      if(!strcmp(aux,"sinh"))  val=(exp(val)-exp(-val))/2;   else
      if(!strcmp(aux,"cosh"))  val=(exp(val)+exp(-val))/2;   else
      if(algebra_error_flag)
        on_error("algebraic evaluation", "undefined function", aux); else
        on_warning("algebraic evaluation", "undefined function", aux);
    } else
    if(algebra_find(aux)) val=algebra_var->a; else
    if(!strcmp(aux,"pi")) val=pi; else
    if(!strcmp(aux,"PI")) val=pi; else
    if(!strcmp(aux,"Pi")) val=pi; else
    if(!strcmp(aux,"e"))  val=exp(1.0); else
    if(!strcmp(aux,"E"))  val=exp(1.0); else
    if(!strcmp(aux,"eV"))  val=1/au_eV; else
    if(!strcmp(aux,"nm"))  val=nm; else
    if(!strcmp(aux,"c"))   val=c_au; else
    if(!strcmp(aux,"m"))   val=1e9*nm; else
    #ifdef jga_constants
    if(!strcmp(aux,"s"))   val=1/au_s; else
    if(!strcmp(aux,"kg"))  val=1/au_kg; else
    if(!strcmp(aux,"kB"))  val=kb_au; else
    #endif
    if(!strcmp(aux,"infinite")) val=infinity; else
    if(!strcmp(aux,"infinity")) val=infinity; else
    if(!strcmp(aux,"zero"))     val=0;        else
    if(algebra_error_flag)
    on_error("algebraic evaluation", "undefined variable", aux); else
    on_warning("algebraic evaluation", "undefined variable", aux);
    delete [] aux;  return val;
  }

  val=read_numero(aux);
  delete [] aux;  return val;
}

numero algebra_eval(char *ex)
{
  int i,l=strlen(ex);  char c;

  for(i=0; i<l; i++) {
    c=ex[0];
    if(c!='+' && c!='-' && c!='*' && c!=':' && c!='/' && c!='^' && c!='.' &&
       c!='(' && c!=')' && c!='[' && c!=']' && c!='{' && c!='}' &&
       (c<'a' || c>'z') && (c<'A' || c>'Z') && (c<'0' || c>'9'))
      algebra_error
  }

  return algebra_eval(ex,0,l-1);
}

int algebra_verbose=0;

void algebra_define(char *name, char *ex)
{
  numero val=algebra_eval(ex);
  if(algebra_find(name)==0) {
    algebra_var=new algebra_variable;
    algebra_var->name=new char [strlen(name)+1];
    strcpy(algebra_var->name,name);
    algebra_var->next=algebra_var0;
    algebra_var0=algebra_var;
  }
  algebra_var->a=val;
  if(algebra_verbose)
    printf(">>> %s = %s = %g <<<\n", algebra_var->name,ex,algebra_var->a);
}

numero alread_numero(char *name)
{
  strcpy(read_command,name);
  return algebra_eval(read_command);
}

int alread_int(char *name)
{
  strcpy(read_command,name);
  return ((int) algebra_eval(read_command));
}

numero alread_numero(FILE *f)
{
  read_name(f,read_command);
  return algebra_eval(read_command);
}

int alread_int(FILE *f)
{
  read_name(f,read_command);
  return ((int) algebra_eval(read_command));
}

#endif  // ******************************************************************
