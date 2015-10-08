/*
Causal Viscous Hydro Code for Non-Central Heavy Ion Collisions

by

Ulrike Romatschke and Paul Romatschke

version 0.2 February 2009
...and Matt Luzum
version 0.3 April 2010
version 0.4 September 2010


The causal viscous hydro equations have been
derived in 

R. Baier, P. Romatschke, D.T. Son, A. Starinets, \
M. Stephanov, arXiv:0712.2451,
(JHEP 0804:100,2008).

Results from this code were used in 

   P.~Romatschke and U.~Romatschke, arXiv:0706.1522,
   (Phys. Rev. Lett.99, 172301,2007).

The setup and tests are documented in

   M.~Luzum and P.~Romatschke, arXiv:0804.4015
   (Phys.Rev.C78, 034915, 2008)

Some parts of the code, in particular the paramreader module, 
were salvaged and modified from unrelated work by one of the 
authors (PR) with Michael Strickland (MS) and in fact were 
originally written by MS.

If you use this code or some part of it, be sure to reference these
articles. Also, if you use the code package as it is, be sure to
refer to the following articles:

* Equation of State: 
M.~Laine and Y.~Schroder, Phys.\ Rev.\  D {\bf 73} (2006) 085009 [arXiv:hep-ph/0603048].

* Resonance Decay Routines: 
J.~Sollfrank, P.~Koch, U.W.~Heinz,  Z.\ Phys.\  C {\bf 52} (1991) 593.
J.~Sollfrank, P.~Koch and U.~W.~Heinz, hys.\ Lett.\  B {\bf 252} (1990) 256.

* If you use the Color Glass Condensate initial conditions, refer to 

H.~J.~Drescher, A.~Dumitru, A.~Hayashigaki and Y.~Nara, Phys.\ Rev.\  C {\bf 74} (2006) 044905 [arXiv:nucl-th/0605012].

Permission to copy this code is granted provided you keep this disclaimer.

*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>


using namespace std;

// these global vars are initialized from parameters file
// defaults set here are overridden by that file

int NUMT=8;
long int STEPS=4000,UPDATE=100,SNAPUPDATE=1000;
double AT=0.05,EPS=0.001,B=0.0;
double ETAOS=0.3;
double TSTART=0.5,TF=0.1,TINIT=1.0;
double IC;
int PTASIZE,PHIPASIZE;

//controls value of tau_Pi
double COEFF=3.0;

//controls value of lambda_1
double L1COEF=2.0;
double L2COEF=1.0;

//create freeze-out surface
int FREEZE=1;

double PTMAX, TRIEPS = 0.0, TRIANGLE = 0.0;
double QUADEPS = 0.0, QUADANGLE = 0.0;
double QUINTEPS = 0.0, QUINTANGLE = 0.0;
double SEXEPS = 0.0, SEXANGLE = 0.0;
double SEPTEPS = 0.0, SEPTANGLE = 0.0;
double BIEPS = 0.0, BIANGLE = 0.0;

//flags
int wflag=0;
int reachedTf=0;
int bflag=0;

int dof=37;

//also for convenience
long int globali;
double globalx;
//defined in loadeos
long int length;

//radius of nucleus in fm
double Rnuc=6.4;
//wood-saxon parameter in fm;
double anuc=0.54;


//for the equation of state
double *eoT4,*cs2i,*poT4,*Ti;

// these hold the current values
double ***u,**e,**pixy,**pixx,**piyy;

// these hold the updated values
double ***U,**E,**Pixy,**Pixx,**Piyy;

// these hold the values from the last UPDATE time step
double ***ulast,**elast,**pixylast,**pixxlast,**piyylast;

//overall time
double t = 0;

//these are global for convenience; used in doInc
double **dtmat;

double **vec;/*bei Pauli rhs*/

//center of lattice 
int Middle;




//for convenience 
//global definition of ut(sx,sy) -- in order not to have it calculated
//a gazillion times
double globut;

double thf[4];
double dtpixx[4],dtpixy[4],dtpiyy[4];
double mypixt,mypiyt,mypitt,mypiee;

// output files
fstream freeze_out;
fstream meta;
fstream ecces;

//splines -- for fancy freeze-out
gsl_interp_accel *wac;
gsl_spline *workspline;
//splines -- for equation of state
gsl_interp_accel *pacc,*Tacc,*cs2acc;
gsl_spline *pspline,*Tspline,*cs2spline;

//to know where to stop interpolation
double lowestE;

double eos(double mye);
double T(int sx,int sy);
double Tlast(int sx,int sy);


// initialize global arrays
void allocateMemory() {

cout << "==> Allocating memory\n";


 u = new double**[2];

 for (int i=0;i<2;i++) 
   u[i] = new double*[NUMT+2];

 for (int i=0;i<2;i++)
   for (int j=0;j<NUMT+2;j++) 
     u[i][j] = new double[NUMT+2];

	 
 e = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
    e[i] = new double[NUMT+2];
 
 pixy = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pixy[i] = new double[NUMT+2];


 pixx = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pixx[i] = new double[NUMT+2];

 piyy = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   piyy[i] = new double[NUMT+2];

//////////////////////
 ulast = new double**[2];

 for (int i=0;i<2;i++) 
   ulast[i] = new double*[NUMT+2];

 for (int i=0;i<2;i++)
   for (int j=0;j<NUMT+2;j++) 
     ulast[i][j] = new double[NUMT+2];

	 
 elast = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
    elast[i] = new double[NUMT+2];
 
 pixylast = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pixylast[i] = new double[NUMT+2];


 pixxlast = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pixxlast[i] = new double[NUMT+2];

 piyylast = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   piyylast[i] = new double[NUMT+2];
//////////////////////

 U = new double**[2];

 for (int i=0;i<2;i++) 
   U[i] = new double*[NUMT+2];

 for (int i=0;i<2;i++)
   for (int j=0;j<NUMT+2;j++) 
     U[i][j] = new double[NUMT+2];

	 
 E = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   E[i] = new double[NUMT+2];


 Pixy = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   Pixy[i] = new double[NUMT+2];


 Pixx = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   Pixx[i] = new double[NUMT+2];


 Piyy = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   Piyy[i] = new double[NUMT+2];

 
 dtmat = new double*[3];
 for (int i=0;i<3;i++) dtmat[i] = new double[3];
 
 vec = new double*[3];
 for (int i=0;i<3;i++) vec[i] = new double[1];
  

}

//enforce periodic BoundaryConditions

void enforcePBCs()
{
  
  for(int sy=1;sy<=NUMT;sy++)
    {
      u[0][0][sy]=u[0][NUMT][sy];
      u[0][NUMT+1][sy]=u[0][1][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      u[0][sx][0]=u[0][sx][NUMT];
      u[0][sx][NUMT+1]=u[0][sx][1];
    }
  for(int sy=1;sy<=NUMT;sy++)
    {
      u[1][0][sy]=u[1][NUMT][sy];
      u[1][NUMT+1][sy]=u[1][1][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      u[1][sx][0]=u[1][sx][NUMT];
      u[1][sx][NUMT+1]=u[1][sx][1];
   }
  for(int sy=1;sy<=NUMT;sy++)
    {
      e[0][sy]=e[NUMT][sy];
      e[NUMT+1][sy]=e[1][sy]; 
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      e[sx][0]=e[sx][NUMT];
      e[sx][NUMT+1]=e[sx][1];
    }

 for(int sy=1;sy<=NUMT;sy++)
    {
      pixx[0][sy]=pixx[NUMT][sy];
      pixx[NUMT+1][sy]=pixx[1][sy];     
    }

  for(int sx=0;sx<=NUMT+1;sx++)
    {
      pixx[sx][0]=pixx[sx][NUMT];
      pixx[sx][NUMT+1]=pixx[sx][1];
    }

 for(int sy=1;sy<=NUMT;sy++)
    {
      pixy[0][sy]=pixy[NUMT][sy];
      pixy[NUMT+1][sy]=pixy[1][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      pixy[sx][0]=pixy[sx][NUMT];
      pixy[sx][NUMT+1]=pixy[sx][1];
    }

 for(int sy=1;sy<=NUMT;sy++)
    {
      piyy[0][sy]=piyy[NUMT][sy];
      piyy[NUMT+1][sy]=piyy[1][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      piyy[sx][0]=piyy[sx][NUMT];
      piyy[sx][NUMT+1]=piyy[sx][1];
    }
}

//Copy fields at each UPDATE time step
void copyUPDATE()
{
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	for (int i=0;i<2;i++) 
	  ulast[i][sx][sy]=U[i][sx][sy];
	elast[sx][sy]=E[sx][sy];
	pixylast[sx][sy]=Pixy[sx][sy];
	pixxlast[sx][sy]=Pixx[sx][sy];
	piyylast[sx][sy]=Piyy[sx][sy];
      }
}

//prepares next time-step
void copyDown() 
{
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	for (int i=0;i<2;i++) 
	  u[i][sx][sy]=U[i][sx][sy];
	e[sx][sy]=E[sx][sy];
	pixy[sx][sy]=Pixy[sx][sy];
	pixx[sx][sy]=Pixx[sx][sy];
	piyy[sx][sy]=Piyy[sx][sy];
      }
  
enforcePBCs();
}


//load equation of state
void loadeos()
{
  fstream eosf;

  //ideal equation of state
  //eosf.open("idEOS.dat", ios::in);

  //qcd equation of state from Mikko Laine and York Schroeder, 
  //hep-ph/0603048
  //eosf.open("qcdEOS.dat",ios::in);
  //interpolated equation of state -- numerically more stable
  eosf.open("qcdIEOS.dat",ios::in);

  if (eosf.is_open())
    {
      eosf >> length;
      cout << "Read length of " << length << endl;

      

      eoT4 = new double[length];
      cs2i = new double[length];
      poT4 = new double[length];
      Ti = new double[length];

      for (int i=1;i<=length;i++)
	{
	  eosf >> Ti[i-1];
	  eosf >> eoT4[i-1];
	  eosf >> poT4[i-1];
	  eosf >> cs2i[i-1];
	}

      eosf.close();
    }
  else
    cout << "Could not open EOS file" << endl; 

  //interpolate
  pacc=gsl_interp_accel_alloc (); 
  Tacc=gsl_interp_accel_alloc (); 
  cs2acc=gsl_interp_accel_alloc (); 
  pspline=gsl_spline_alloc (gsl_interp_cspline, length);
  Tspline=gsl_spline_alloc (gsl_interp_cspline, length);
  cs2spline=gsl_spline_alloc (gsl_interp_cspline, length);
  
  
  double warrx[length];
  double warry[length];

  for (int i=0;i<length;i++)
    {
      warrx[i]=eoT4[i]*Ti[i]*Ti[i]*Ti[i]*Ti[i];
      warry[i]=poT4[i]*Ti[i]*Ti[i]*Ti[i]*Ti[i];
    }

  lowestE=warrx[0];

  gsl_spline_init (pspline,warrx,warry,length);
  
  for (int i=0;i<length;i++)
    {
      warry[i]=Ti[i];
    }  
  gsl_spline_init (Tspline,warrx,warry,length);
  
  for (int i=0;i<length;i++)
    {
      warry[i]=cs2i[i];
    }  
  gsl_spline_init (cs2spline,warrx,warry,length);

}

//return interpolated temperature in lattice units
double getintT(int i,double x)
{
  return (Ti[i]+x*(Ti[i+1]-Ti[i]))*AT;
}

//Wood-Saxon routine, in physical units
double WS(double x, double y, double z)
{
  
  //cout <<"R= " << R << endl;

  double temp;
  temp=x*x+y*y+z*z;
  temp=sqrt(temp);
  temp-=Rnuc;
  temp/=anuc;
  return 1/(1+exp(temp));
}

//transverse density T_A, in physical units
double TA(double x, double y)
{
  double temp=0;
  for (int i=0;i<1200;i++)
    temp+=WS(x,y,(i+0.5)*Rnuc/400.)*Rnuc/400.;

  //return result normalized to a gold nucleus
  return 2*temp*197./1175.22;
}


//number density of wounded nucleons
//from Kolb et. al, hep-ph/0103234
double getwnuc(double xx,double yy,double b)
{
  
  double mTAp=TA(xx+b/2.,yy);
  double mTAm=TA(xx-b/2.,yy);

  //return mTA*2*(1.-pow(1-mTA/197.*4.,197.));
  double temp=0;
  temp+=mTAp*(1.-exp(-mTAm*4.));
  temp+=mTAm*(1.-exp(-mTAp*4.));
  return temp;
}


//reads in initial energy density profile
//from inited.dat (generate by either initE.cpp or your
//favorite routine)
void setInitialConditions()
{

  double sig;
  //extern double randGauss(double); // generates a gaussian random number with mean 0

  //Middle
  Middle=(NUMT-1)/2+1;
  //Middle=1;

  printf("mi %i\n",Middle);

  //freeze-out temperature
  TF=TF*AT;

  cout << "TF=" << TF/AT << endl;

  cout << "TSTART=" << TSTART << endl;

  //load equation of state

  loadeos();

  //convert fm/c to lattice units
  t=TINIT*5.06842/AT;

  fstream inited;
  inited.open("data/inited.dat", ios::in);

  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	u[0][sx][sy]=0.0;
	u[1][sx][sy]=0.0;
	//e[sx][sy]=e0*getwnuc((sx-Middle)*AT/5.06842,(sy-Middle)*AT/5.06842,B)/4.29048;
	inited >> e[sx][sy];
	//double NS=ETAOS*2./(3.*t)*(e[sx][sy]+eos(e[sx][sy]))/(T(sx,sy)+0.05);
	double NS=0;
	pixx[sx][sy]=NS;
	pixy[sx][sy]=0;
	piyy[sx][sy]=NS;
      }
  
  wac=gsl_interp_accel_alloc (); 
  workspline=gsl_spline_alloc (gsl_interp_cspline, Middle);
  
  inited.close();

  enforcePBCs();
}

//gets index associated with energy density -- internal use only
long int geti(double mye)
{
  long int i;
  mye/=AT*AT*AT*AT;
  for (i=0;i<length;i++)
    {
      if (eoT4[i]*Ti[i]*Ti[i]*Ti[i]*Ti[i]>mye)
	break;
    }
  return (i-1);
}

//get fraction associated with index i-- internal use only
double getx(long int i,double mye)
{
  mye/=AT*AT*AT*AT;
  
  double temp;

  temp=(eoT4[i+1]-eoT4[i])/(Ti[i+1]-Ti[i])*(Ti[i+1]+Ti[i])/2;
  temp+=2*(eoT4[i+1]+eoT4[i]);
  temp*=(Ti[i+1]+Ti[i])/2;
  temp*=(Ti[i+1]+Ti[i])/2;
  temp*=(Ti[i+1]+Ti[i])/2;

  return (mye-eoT4[i]*Ti[i]*Ti[i]*Ti[i]*Ti[i])/temp/(Ti[i+1]-Ti[i]);
}


//-----------------------------------------------------------------------------
//declaration of functions, shells for functions (for faster running)
//and spatial derivatives

double ut(int sx,int sy)
{
  double temp=1.;
  temp+=u[0][sx][sy]*u[0][sx][sy];
  temp+=u[1][sx][sy]*u[1][sx][sy];
  return sqrt(temp);
}

double umu(int mu, int sx,int sy)
{
  if (mu<2)
    return u[mu][sx][sy];
  if (mu==2)
    return ut(sx,sy);
  else
    return 0;
}



//returns correct pi
double pi (int delta,int beta,int sx,int sy)
{

  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta

  int phi=0;
  int deltap,betap;

  if (beta==2)
    betap=0;
  if (beta==3)
    betap=3;
  if (beta==1)
    betap=2;
  if (beta==0)
    betap=1;
  if (delta==2)
    deltap=0;
  if (delta==3)
    deltap=3;
  if (delta==1)
    deltap=2;
  if (delta==0)
    deltap=1;
  

  if (deltap>betap)
    {
      phi=betap;
      betap=deltap;
      deltap=phi;
    }
  if (deltap<betap)
    {
      if (betap==1 && deltap==0)
	{
	  return u[0][sx][sy]/ut(sx,sy)*pixx[sx][sy]+u[1][sx][sy]/ut(sx,sy)*pixy[sx][sy];
	}
      if (betap==2)
	{
	  if (deltap==0)
	    return u[0][sx][sy]/ut(sx,sy)*pixy[sx][sy]+u[1][sx][sy]/ut(sx,sy)*piyy[sx][sy];
	  if (deltap==1)
	    return pixy[sx][sy];
	}
      else return 0;
    }

  if (deltap==betap)
    {
      if (deltap==0)
	return u[0][sx][sy]/ut(sx,sy)*u[0][sx][sy]/ut(sx,sy)*pixx[sx][sy]+2*u[0][sx][sy]/ut(sx,sy)*u[1][sx][sy]/ut(sx,sy)*pixy[sx][sy]+u[1][sx][sy]/ut(sx,sy)*u[1][sx][sy]/ut(sx,sy)*piyy[sx][sy];
      if (deltap==1)
	return pixx[sx][sy];
      if (deltap==2)
	return piyy[sx][sy];
      if (deltap==3)
	return (-1)/(t*t)*((-1)*u[0][sx][sy]/ut(sx,sy)*u[0][sx][sy]/ut(sx,sy)*pixx[sx][sy]-2*u[0][sx][sy]/ut(sx,sy)*u[1][sx][sy]/ut(sx,sy)*pixy[sx][sy]-u[1][sx][sy]/ut(sx,sy)*u[1][sx][sy]/ut(sx,sy)*piyy[sx][sy]+pixx[sx][sy]+piyy[sx][sy]);
    }
}


//to make things faster
double pishell (int delta,int beta,int sx,int sy)
{
  int phi=0;
  if (delta>beta)
    {
      phi=beta;
      beta=delta;
      delta=phi;
    }
  if (delta==0&&beta==0)
    return pixx[sx][sy];
  if (delta==0&&beta==1)
    return pixy[sx][sy];
  if (delta==0&&beta==2)
    return mypixt;
  if (delta==0&&beta==3)
    return 0;
  if (delta==1&&beta==1)
    return piyy[sx][sy];
  if (delta==1&&beta==2)
    return mypiyt;
  if (delta==1&&beta==3)
    return 0;
  if (delta==2&&beta==2)
    return mypitt;
  if (delta==2&&beta==3)
    return 0;
  if (delta==3&&beta==3)
    return mypiee;
  
}


//this provides dx u[i]
double dxu(int i,int sx,int sy)
{
  double temp=0;
  if(sx!=1)
    {
      if(sx==NUMT)
	{
	  temp=(umu(i,sx,sy)-umu(i,sx-1,sy));
	}
      else
	temp=(umu(i,sx+1,sy)-umu(i,sx-1,sy))/2;
    }
  else
    temp=(umu(i,sx+1,sy)-umu(i,sx,sy));
  return temp;
}

//this provides dy u[i]
double dyu(int i,int sx,int sy)
{
  double temp=0;
  if(sy!=1)
    {
      if(sy==NUMT)
	{
	  temp=(umu(i,sx,sy)-umu(i,sx,sy-1));
	}
      else
	temp=(umu(i,sx,sy+1)-umu(i,sx,sy-1))/2;
    }
  else
    temp=(umu(i,sx,sy+1)-umu(i,sx,sy));
  return temp;
}

//provides \partial_i u^\mu
double diumu(int i,int mu,int sx,int sy)
{
  double temp=0;
  if (mu!=2)
    {
      if (i==0)
	temp=dxu(mu,sx,sy);
      if (i==1)
	temp=dyu(mu,sx,sy);
    }
  if (mu==2)
    {
      if (i==0)
	temp=umu(0,sx,sy)*dxu(0,sx,sy)+umu(1,sx,sy)*dxu(1,sx,sy);
      if (i==1)
	temp=umu(0,sx,sy)*dyu(0,sx,sy)+umu(1,sx,sy)*dyu(1,sx,sy);
      temp/=umu(2,sx,sy);
    }
  return temp;
}

//this provides dx e
double dxe(int sx,int sy)
{
  double temp=0;
  if(sx!=1)
    {
      if(sx==NUMT)
	{
	  temp=(e[sx][sy]-e[sx-1][sy]);
	}
      else
	temp=(e[sx+1][sy]-e[sx-1][sy])/2;
    }
  else
    temp=(e[sx+1][sy]-e[sx][sy]);
  return temp;
}


//this provides dy e
double dye(int sx,int sy)
{
  double temp=0;
  if(sy!=1)
    {
      if(sy==NUMT)
	{
	  temp=(e[sx][sy]-e[sx][sy-1]);
	}
      else
	temp=(e[sx][sy+1]-e[sx][sy-1])/2;
    }
  else
    temp=(e[sx][sy+1]-e[sx][sy]);
  return temp;
}


//term f
//di pi^mu,alpha
double djpi(int j,int mu,int alpha,int sx,int sy)
{
  double temp=0;
  if (j==0)
    {
      if(sx!=1)
	{
	  if(sx==NUMT)
	    {
	      temp=(pi(mu,alpha,sx,sy)-pi(mu,alpha,sx-1,sy));
	    }
	  else
	    temp=(pi(mu,alpha,sx+1,sy)-pi(mu,alpha,sx-1,sy))/2;
	}
      else
	temp=(pi(mu,alpha,sx+1,sy)-pi(mu,alpha,sx,sy));
    }
  if (j==1)
    {
      if(sy!=1)
	{
	  if(sy==NUMT)
	    {
	      temp=(pi(mu,alpha,sx,sy)-pi(mu,alpha,sx,sy-1));
	    }
	  else
	    temp=(pi(mu,alpha,sx,sy+1)-pi(mu,alpha,sx,sy-1))/2;
	}
      else
	temp=(pi(mu,alpha,sx,sy+1)-pi(mu,alpha,sx,sy));
    }
  return temp;
}

//returns all non zero gamma
double gamma (int alpha, int beta, int delta)
{

  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
  if (beta==3)
    {
    if (alpha==2 &&delta==3)
      return t;
    else
      if (alpha==3 && delta==2)
	return 1/t;
      else 
	return 0;
    }
  if (beta==2 && alpha==3 && delta==3)
    {
      return 1/t;
    }
  else
    return 0;
}


//metric 2 upper indices
double g(int alpha,int beta)
{
  if (alpha!=beta)
    return 0.;
  else
    {
      if (alpha==0 || alpha==1)
	return -1.;
      if (alpha==2)
	return 1.;
      if (alpha==3)
	return (-1)/(t*t);
    }
}

//metric 2 lower indices
double gdown(int alpha,int beta)
{
  if (alpha!=beta)
    return 0.;
  else
    {
      if ((alpha==0)||(alpha==1))
	return -1.;
      if (alpha==2)
	return 1.;
      if (alpha==3)
	return -(t*t);
    }
}

//Delta^{mu kappa}
double Delta(int mu,int kappa,int sx,int sy)
{
   return g(mu,kappa)- umu(mu,sx,sy)*umu(kappa,sx,sy);
}

//----------------------------------------------------------------------
//thermodynamic functions

double eos(double mye)
{
  double temp=0;

  double phys=mye/AT/AT/AT/AT;

  if (phys>lowestE)
    {
      temp=gsl_spline_eval(pspline,phys,pacc);
      temp*=AT*AT;
      temp*=AT*AT;
    }
  else
    {
      temp=gsl_spline_eval(cs2spline,lowestE,cs2acc);
      temp*=phys;
      temp*=AT*AT;
      temp*=AT*AT;
    }

  return temp;
}


//Speed of Sound squared = dp/depsilon
double cs2(int sx,int sy)
{

  double temp;
  double phys=e[sx][sy]/AT/AT/AT/AT;

  if (phys>lowestE)
    {
      temp=gsl_spline_eval(cs2spline,phys,cs2acc);
    }
  else
    {
      temp=gsl_spline_eval(cs2spline,lowestE,cs2acc);
    }
  return temp;
}

//provides Temperature*lattice spacing
double T(int sx,int sy)
{
  if (e[sx][sy]<0)
    {
      printf("Negative e at sx=%i sy=%i\n",sx,sy);
    } 

  double temp;
  double phys=e[sx][sy]/AT/AT/AT/AT;

  if (phys>lowestE)
    {
      temp=gsl_spline_eval(Tspline,phys,Tacc);
      temp*=AT;
    }
  else
    {
      temp=sqrtl(sqrtl(e[sx][sy]/eoT4[0]));
    }
  return temp;

}

double Tlast(int sx,int sy)
{

  double temp;
  double phys=elast[sx][sy]/AT/AT/AT/AT;

  if (phys>lowestE)
    {
      temp=gsl_spline_eval(Tspline,phys,Tacc);
      temp*=AT;
    }
  else
    {
      temp=sqrtl(sqrtl(elast[sx][sy]/eoT4[0]));
    }
  return temp;

}

//Provides tau_Pi/lattice spacing;
double taupi(int sx,int sy)
{
  double temp=COEFF*2.0*ETAOS/T(sx,sy);
  //printf("E %f\n",e[sx][sy]);
  if (isnan(temp)!=0)
    {
      cout << "Error in taupi\n";
    }
  return temp;
}

//eta/taupi
double etataupi(int sx,int sy)
{
  double temp=0.5/COEFF*(e[sx][sy]+eos(e[sx][sy]));
  return temp;
}



//-----------------------------------------------------------------------
//functions that are used to fill in the update matrix Eq.(9,10) of 
//nucl-th/0610108


//convention:
//coefficients denote
//a[0]-> coeff multiplying partial_t u[0]
//a[1]-> coeff multiplying partial_t u[1]
//a[2]-> coeff multiplying partial_t e
//a[3]-> remainder


//this provides D u^mu
void Dumu(double *a,int mu,int sx,int sy)
{
  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
 if (mu==0)
   {
     a[0]=globut;
     a[1]=0;
     a[2]=0;
     a[3]=u[0][sx][sy]*dxu(0,sx,sy);
     a[3]+=u[1][sx][sy]*dyu(0,sx,sy);
}

 if (mu==1)
   {
     a[0]=0;
     a[1]=globut;
     a[2]=0;
     a[3]=u[0][sx][sy]*dxu(1,sx,sy);
     a[3]+=u[1][sx][sy]*dyu(1,sx,sy);
   }
 
 if (mu==2)
   {
     a[0]=u[0][sx][sy];
     a[1]=u[1][sx][sy];
     a[2]=0;
     a[3]=u[0][sx][sy]*dxu(2,sx,sy);
     a[3]+=u[1][sx][sy]*dyu(2,sx,sy);
   }
 if (mu==3)
   {
     a[0]=0;
     a[1]=0;
     a[2]=0;
     a[3]=0;
   }
}



//Second term Grad p:


//this provides \nabla^\mu p
void Nablap(double *a,int mu,int sx,int sy)
{
  //unconventional notation: 0=x, 1=y, 2=tau
  if (mu==0)
    {
      a[0]=0;
      a[1]=0;
      a[2]=(-1)*cs2(sx,sy)*u[0][sx][sy]*globut;
      a[3]=(-1)*cs2(sx,sy)*u[0][sx][sy]*u[0][sx][sy]*dxe(sx,sy);
      a[3]-=cs2(sx,sy)*u[0][sx][sy]*u[1][sx][sy]*dye(sx,sy);
      a[3]-=cs2(sx,sy)*dxe(sx,sy);
    }
  if (mu==1)
    {
      a[0]=0;
      a[1]=0;
      a[2]=(-1)*cs2(sx,sy)*u[1][sx][sy]*globut;
      a[3]=(-1)*cs2(sx,sy)*u[1][sx][sy]*u[0][sx][sy]*dxe(sx,sy);
      a[3]-=cs2(sx,sy)*u[1][sx][sy]*u[1][sx][sy]*dye(sx,sy);
      a[3]-=cs2(sx,sy)*dye(sx,sy);
    }
 
}


//Third term:


//Second part of third term with gammas:


//g*Delta*gammas

double gDeltagamma(int mu,int sx,int sy)
{
  double temp=0;
  double h1=0;
  double h2=0;
  for (int kappa=0;kappa<=3;kappa++)
    for (int alpha=0;alpha<=3;alpha++)
      for (int beta=0;beta<=3;beta++)
	for (int delta=0;delta<=3;delta++)
	  {
	    h1=gamma(alpha,beta,delta);
	    h2=gamma(beta,beta,delta);
	    if (h1!=0)
	      temp+=Delta(mu,kappa,sx,sy)*gdown(alpha,kappa)*(gamma(alpha,beta,delta)*pishell(delta,beta,sx,sy));
	    if (h2!=0)
	      temp+=Delta(mu,kappa,sx,sy)*gdown(alpha,kappa)*(gamma(beta,beta,delta)*pishell(alpha,delta,sx,sy));
	  }
  return temp;
}

//First part of third term:

//From eq. (13)

//d_t pi^i,alpha
//term a
//double vorticity(.....)
//{
//}
//term b
//<Grad u> upper indices!

double theta3(int sx,int sy)
{
  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
  double temp=0;
  for (int kappa=0;kappa<=3;kappa++)
    {
      {
	for (int alpha=0;alpha<=3;alpha++)
	  for (int i=0;i<=1;i++)
	    temp+=gdown(alpha,kappa)*Delta(i,kappa,sx,sy)*diumu(i,alpha,sx,sy);
      }
      temp+=gdown(3,kappa)*Delta(3,kappa,sx,sy)*gamma(3,2,3)*globut;
    }
  /*
    if ((sx==Middle+8)&&(sy==Middle))
    {
    double t2=0;
    for (int kappa=0;kappa<=3;kappa++)
    for (int alpha=0;alpha<=3;alpha++)
    for (int i=0;i<=1;i++)
    t2+=gdown(alpha,kappa)*Delta(i,kappa,sx,sy)*diumu(i,alpha,sx,sy);
    
    cout << "theta " << temp;
    cout << "\t where" << t2;
    cout << endl; 
      
    }
  */
  return temp;
}

double theta0(int sx,int sy)
{
  double theta0=Delta(2,2,sx,sy)*u[0][sx][sy]/globut-Delta(2,0,sx,sy);
  return theta0;
}

double theta1(int sx ,int sy)
{
  double theta1=-Delta(2,1,sx,sy)+Delta(2,2,sx,sy)*u[1][sx][sy]/globut;;
  return theta1;
}


//<nabla_x u_x>
void gradxux(double *a,int sx,int sy)
{
  
  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
  a[0]=2*Delta(2,0,sx,sy)-2/3.*Delta(0,0,sx,sy)*thf[0];
  a[1]=(-1)*2/3.*Delta(0,0,sx,sy)*thf[1];
  a[2]=0;
  a[3]=2*(Delta(0,0,sx,sy)*dxu(0,sx,sy)+Delta(0,1,sx,sy)*dyu(0,sx,sy));
  a[3]-=2/3.*Delta(0,0,sx,sy)*thf[3];
  /*
    if ((sx==Middle+8)&&(sy==Middle))
    {
    cout << "here gradxux" << 2./3.*Delta(0,0,sx,sy)*thf[3];
    cout << "\t other" << -2*(Delta(0,0,sx,sy)*dxu(0,sx,sy)+Delta(0,1,sx,sy)*dyu(0,sx,sy));
    cout << endl;
    }
  */

}

void gradxuy(double *a,int sx,int sy)
{
  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
  a[0]=Delta(2,1,sx,sy)-2/3.*Delta(0,1,sx,sy)*thf[0];
  a[1]=Delta(2,0,sx,sy)-2/3.*Delta(0,1,sx,sy)*thf[1];
  a[2]=0;
  a[3]=Delta(0,0,sx,sy)*dxu(1,sx,sy);
  a[3]+=Delta(0,1,sx,sy)*dyu(1,sx,sy);
  a[3]+=Delta(0,1,sx,sy)*dxu(0,sx,sy);
  a[3]+=Delta(1,1,sx,sy)*dyu(0,sx,sy);
  a[3]-=2/3.*Delta(0,1,sx,sy)*thf[3];
}

//<nabla_x u_t>
void gradxut(double *a,int sx,int sy)
{
  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
  /*a[0]=umu(0,sx,sy)/globut*(2*Delta(2,0,sx,sy)-2/3.*Delta(0,0,sx,sy)*thf[0]);
  a[0]+=umu(1,sx,sy)/globut*(Delta(2,1,sx,sy)-2/3.*Delta(0,1,sx,sy)*thf[0]); 
  a[1]=umu(0,sx,sy)/globut*(-1)*2/3.*Delta(0,0,sx,sy)*thf[1]; 
  a[1]+=umu(1,sx,sy)/globut*(Delta(2,0,sx,sy)-2/3.*Delta(0,1,sx,sy)*thf[1]); 
  a[2]=0; 
  a[3]=umu(0,sx,sy)/globut*(2*(Delta(0,0,sx,sy)*dxu(0,sx,sy)+Delta(0,1,sx,sy)*dyu(0,sx,sy))-2/3.*Delta(0,0,sx,sy)*thf[3]); 
  a[3]+=umu(1,sx,sy)/globut*(Delta(0,0,sx,sy)*dxu(1,sx,sy)+Delta(0,1,sx,sy)*dyu(1,sx,sy)+Delta(0,1,sx,sy)*dxu(0,sx,sy)+Delta(1,1,sx,sy)*dyu(0,sx,sy)-2/3.*Delta(0,1,sx,sy)*thf[3]);*/

  
  //oder
  
  double grxux[4];
  gradxux(grxux,sx,sy);
  double grxuy[4];
  gradxuy(grxuy,sx,sy);
  
  a[0]=-(umu(0,sx,sy)/globut*grxux[0]+umu(1,sx,sy)/globut*grxuy[0]);
  a[1]=-(umu(0,sx,sy)/globut*grxux[1]+umu(1,sx,sy)/globut*grxuy[1]);
  a[2]=0;
  a[3]=-(umu(0,sx,sy)/globut*grxux[3]+umu(1,sx,sy)/globut*grxuy[3]);
}

void gradyuy(double *a,int sx,int sy)
{
  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
  a[0]=(-1)*2/3.*Delta(1,1,sx,sy)*thf[0];
  a[1]=2*Delta(2,1,sx,sy)-2/3.*Delta(1,1,sx,sy)*thf[1];
  a[2]=0;
  a[3]=2*(Delta(0,1,sx,sy)*dxu(1,sx,sy)+Delta(1,1,sx,sy)*dyu(1,sx,sy));
  a[3]-=2/3.*Delta(1,1,sx,sy)*thf[3];
}

//<nabla_y u_t>
void gradyut(double *a,int sx,int sy)
{
  /*//unconventional notation: 0=x, 1=y, 2=tau, 3=eta
  a[0]=umu(1,sx,sy)/globut*(-1)*2/3.*Delta(1,1,sx,sy)*thf[0]+umu(0,sx,sy)/globut*(Delta(2,1,sx,sy)-2/3.*Delta(0,1,sx,sy)*thf[0]);
  a[1]=umu(1,sx,sy)/globut*(2*Delta(2,1,sx,sy)-2/3.*Delta(1,1,sx,sy)*thf[1]); 
  a[1]+=umu(0,sx,sy)/globut*(Delta(2,0,sx,sy)-2/3.*Delta(0,1,sx,sy)*thf[1]);
  a[2]=0;
  a[3]=umu(1,sx,sy)/globut*(2*(Delta(1,1,sx,sy)*dyu(1,sx,sy)+Delta(0,1,sx,sy)*dxu(1,sx,sy))-2/3.*Delta(1,1,sx,sy)*thf[3]);
  a[3]+=umu(0,sx,sy)/globut*(Delta(0,1,sx,sy)*dyu(1,sx,sy)+Delta(0,0,sx,sy)*dxu(1,sx,sy)+Delta(1,1,sx,sy)*dyu(0,sx,sy)+Delta(0,1,sx,sy)*dxu(0,sx,sy)-2/3.*Delta(0,1,sx,sy)*thf[3]);*/
  
  //oder
  
  double grxuy[4];
  gradxuy(grxuy,sx,sy);
  double gryuy[4];
  gradyuy(gryuy,sx,sy);
    
  a[0]=-(umu(0,sx,sy)/globut*grxuy[0]+umu(1,sx,sy)/globut*gryuy[0]);
  a[1]=-(umu(0,sx,sy)/globut*grxuy[1]+umu(1,sx,sy)/globut*gryuy[1]);
  a[2]=0;
  a[3]=-(umu(0,sx,sy)/globut*grxuy[3]+umu(1,sx,sy)/globut*gryuy[3]);
} 

void gradtut(double *a,int sx,int sy)
{
  /*a[0]=umu(0,sx,sy)/globut*umu(0,sx,sy)/globut*(2*Delta(2,0,sx,sy)-2/3.*Delta(0,0,sx,sy)*thf[0]);
  a[0]+=2*umu(0,sx,sy)/globut*umu(1,sx,sy)/globut*(Delta(2,1,sx,sy)-2/3.*Delta(0,1,sx,sy)*thf[0]);
  a[0]+=umu(1,sx,sy)/globut*umu(1,sx,sy)/globut*(-1)*2/3.*Delta(1,1,sx,sy)*thf[0];
  a[1]=umu(0,sx,sy)/globut*umu(0,sx,sy)/globut*(-1)*2/3.*Delta(0,0,sx,sy)*thf[1];
  a[1]+=2*umu(0,sx,sy)/globut*umu(1,sx,sy)/globut*(Delta(2,0,sx,sy)-2/3.*Delta(0,1,sx,sy)*thf[1]);
  a[1]+=umu(1,sx,sy)/globut*umu(1,sx,sy)/globut*(2*Delta(2,1,sx,sy)-2/3.*Delta(1,1,sx,sy)*thf[1]);
  a[2]=0; 
  a[3]=umu(0,sx,sy)/globut*umu(0,sx,sy)/globut*(2*(Delta(0,0,sx,sy)*dxu(0,sx,sy)+Delta(0,1,sx,sy)*dyu(0,sx,sy))-2/3.*Delta(0,0,sx,sy)*thf[3]);
  a[3]+=2*umu(0,sx,sy)/globut*umu(1,sx,sy)/globut*(Delta(0,0,sx,sy)*dxu(1,sx,sy)+Delta(0,1,sx,sy)*dyu(1,sx,sy)+Delta(0,1,sx,sy)*dxu(0,sx,sy)+Delta(1,1,sx,sy)*dyu(0,sx,sy)-2/3.*Delta(0,1,sx,sy)*thf[3]);
  a[3]+=umu(1,sx,sy)/globut*umu(1,sx,sy)/globut*(2*(Delta(0,1,sx,sy)*dxu(1,sx,sy)+Delta(1,1,sx,sy)*dyu(1,sx,sy))-2/3.*Delta(1,1,sx,sy)*thf[3]);*/
  
  //oder
  
  double grxux[4];
  gradxux(grxux,sx,sy);
  double grxuy[4];
  gradxuy(grxuy,sx,sy);
  double gryuy[4];
  gradyuy(gryuy,sx,sy);
  
  a[0]=umu(0,sx,sy)/globut*umu(0,sx,sy)/globut*grxux[0];
  a[0]+=2*umu(0,sx,sy)/globut*umu(1,sx,sy)/globut*grxuy[0];
  a[0]+=umu(1,sx,sy)/globut*umu(1,sx,sy)/globut*gryuy[0];
  a[1]=umu(0,sx,sy)/globut*umu(0,sx,sy)/globut*grxux[1];
  a[1]+=2*umu(0,sx,sy)/globut*umu(1,sx,sy)/globut*grxuy[1];
  a[1]+=umu(1,sx,sy)/globut*umu(1,sx,sy)/globut*gryuy[1];
  a[2]=0;
  a[3]=umu(0,sx,sy)/globut*umu(0,sx,sy)/globut*grxux[3];
  a[3]+=2*umu(0,sx,sy)/globut*umu(1,sx,sy)/globut*grxuy[3];
  a[3]+=umu(1,sx,sy)/globut*umu(1,sx,sy)/globut*gryuy[3];
} 

void gradeue(double *a,int sx,int sy)
{
  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
  a[0]=(-1)*2/3.*t*t*t*t*Delta(3,3,sx,sy)*thf[0];
  a[1]=(-1)*2/3.*t*t*t*t*Delta(3,3,sx,sy)*thf[1];
  a[2]=0;
  a[3]=2*t*t*t*t*Delta(3,3,sx,sy)*gamma(3,2,3)*globut-2/3.*t*t*t*t*Delta(3,3,sx,sy)*thf[3];
}

//extra terms

//gives D_i u_j
double diuj(int i,int j, int sx,int sy)
{
  double temp;
  if ((i==0)&&(j==0))
    temp=dxu(j,sx,sy);
  if ((i==0)&&(j==1))
    temp=dxu(j,sx,sy);
  if ((i==0)&&(j==3))
    temp=0;
  if ((i==1)&&(j==0))
    temp=dyu(j,sx,sy);
  if ((i==1)&&(j==1))
    temp=dyu(j,sx,sy);
  if ((i==1)&&(j==3))
    temp=0;
  if ((i==3)&&(j==0))
    temp=0;
  if ((i==3)&&(j==1))
    temp=0;
  if ((i==3)&&(j==3))
    temp=-ut(sx,sy)*t;

  return temp;
}


//Omega_{i j}
double prevor(double *a,int i,int j,int sx,int sy)
{
  double bi[4],bj[4];
  Dumu(bi,i,sx,sy);
  Dumu(bj,j,sx,sy);
  
  a[0]=0.5*(umu(i,sx,sy)*bj[0]-umu(j,sx,sy)*bi[0]);
  a[1]=0.5*(umu(i,sx,sy)*bj[1]-umu(j,sx,sy)*bi[1]);
  a[2]=0.5*(umu(i,sx,sy)*bj[2]-umu(j,sx,sy)*bi[2]);
  a[3]=0.5*(umu(i,sx,sy)*bj[3]-umu(j,sx,sy)*bi[3]);
  a[3]-=0.5*(diuj(i,j,sx,sy)-diuj(j,i,sx,sy));

}


//Omega_{xy}
void vorticityxy(double *a,int sx,int sy)
{
  double bx[4],by[4];
  Dumu(bx,0,sx,sy);
  Dumu(by,1,sx,sy);

  a[0]=0.5*(u[1][sx][sy]*bx[0]-u[0][sx][sy]*by[0]);
  a[1]=0.5*(u[1][sx][sy]*bx[1]-u[0][sx][sy]*by[1]);
  a[2]=0.5*(u[1][sx][sy]*bx[2]-u[0][sx][sy]*by[2]);
  a[3]=0.5*(u[1][sx][sy]*bx[3]-u[0][sx][sy]*by[3]);
  a[3]-=0.5*(dxu(1,sx,sy)-dyu(0,sx,sy));

  
}


//gives \Omega_{mu,nu}
double vorticity(double *a,int mu, int nu,int sx,int sy)
{
  if ((mu==2)||(nu==2))
    {
      if ((mu==2)&&(nu==2))
	{
	  a[0]=0;
	  a[1]=0;
	  a[2]=0;
	  a[3]=0;
	}
      if ((mu==2)&&(nu!=2))
	{
	  double bx[4],by[4];
	  prevor(bx,nu,0,sx,sy);
	  prevor(by,nu,1,sx,sy);
	  a[0]=(u[0][sx][sy]*bx[0]+u[1][sx][sy]*by[0])/ut(sx,sy);
	  a[1]=(u[0][sx][sy]*bx[1]+u[1][sx][sy]*by[1])/ut(sx,sy);
	  a[2]=(u[0][sx][sy]*bx[2]+u[1][sx][sy]*by[2])/ut(sx,sy);
	  a[3]=(u[0][sx][sy]*bx[3]+u[1][sx][sy]*by[3])/ut(sx,sy);
	}
      if ((mu!=2)&&(nu==2))
	{
	  double bx[4],by[4];
	  prevor(bx,0,mu,sx,sy);
	  prevor(by,1,mu,sx,sy);
	  a[0]=(u[0][sx][sy]*bx[0]+u[1][sx][sy]*by[0])/ut(sx,sy);
	  a[1]=(u[0][sx][sy]*bx[1]+u[1][sx][sy]*by[1])/ut(sx,sy);
	  a[2]=(u[0][sx][sy]*bx[2]+u[1][sx][sy]*by[2])/ut(sx,sy);
	  a[3]=(u[0][sx][sy]*bx[3]+u[1][sx][sy]*by[3])/ut(sx,sy);
	}
    }
  else 
    {
      prevor(a,mu,nu,sx,sy);
    }

//  if ((sx==15)&&(sy==15))
  //   {
  //  printf("\n here for %i %i",mu,nu);
  //  printf("\t made a[3]=%.12g\n",a[3]);
  //}
}




void terma(double *a,int i,int alpha, int sx,int sy)
{

  //5/2 D ln T
  //a[0]=0;
  //a[1]=0;
  //a[2]=5/8./e[sx][sy]*ut(sx,sy);
  //a[3]=5/8./e[sx][sy]*(u[0][sx][sy]*dxe(sx,sy)+u[1][sx][sy]*dye(sx,sy));

  

  //subtract 1/2 theta
  
  //a[0]-=0.5*theta0(sx,sy);
  //a[1]-=0.5*theta1(sx,sy);
  //a[3]-=0.5*theta3(sx,sy);
  
  //  -4/3 theta

  a[0]=-4./3.*theta0(sx,sy);
  a[1]=-4./3.*theta1(sx,sy);
  a[2]=0;
  a[3]=-4./3.*theta3(sx,sy);


  //multiply by pi^ialpha

  a[0]*=pishell(i,alpha,sx,sy);
  a[1]*=pishell(i,alpha,sx,sy);
  a[2]*=pishell(i,alpha,sx,sy);
  a[3]*=pishell(i,alpha,sx,sy);
  

  
  //subtract vorticity
  

  double b[4];
  //double b2[4];
  
  vorticity(b,0,1,sx,sy);
  

  
  if (i==alpha)
    {
      for (int rho=0;rho<4;rho++)
	{

	  if (alpha!=rho)
	    {
	      vorticity(b,alpha,rho,sx,sy);
  

	      a[0]-=L2COEF*pishell(rho,i,sx,sy)*b[0]*g(alpha,alpha);
	      a[1]-=L2COEF*pishell(rho,i,sx,sy)*b[1]*g(alpha,alpha);
	      a[2]-=L2COEF*pishell(rho,i,sx,sy)*b[2]*g(alpha,alpha);
	      a[3]-=L2COEF*pishell(rho,i,sx,sy)*b[3]*g(alpha,alpha);
	    }
	}
    }
  else
    {
      for (int rho=0;rho<4;rho++)
	{
	  if (alpha!=rho)
	    {
	      vorticity(b,alpha,rho,sx,sy);
	      
	      a[0]-=0.5*L2COEF*pishell(rho,i,sx,sy)*b[0]*g(alpha,alpha);
	      a[1]-=0.5*L2COEF*pishell(rho,i,sx,sy)*b[1]*g(alpha,alpha);
	      a[2]-=0.5*L2COEF*pishell(rho,i,sx,sy)*b[2]*g(alpha,alpha);
	      a[3]-=0.5*L2COEF*pishell(rho,i,sx,sy)*b[3]*g(alpha,alpha);
	    }
	}
  
      for (int rho=0;rho<4;rho++)
	{
	  if (i!=rho)
	    {
	      vorticity(b,i,rho,sx,sy);
  
	      a[0]-=0.5*L2COEF*pishell(rho,alpha,sx,sy)*b[0]*g(i,i);
	      a[1]-=0.5*L2COEF*pishell(rho,alpha,sx,sy)*b[1]*g(i,i);
	      a[2]-=0.5*L2COEF*pishell(rho,alpha,sx,sy)*b[2]*g(i,i);
	      a[3]-=0.5*L2COEF*pishell(rho,alpha,sx,sy)*b[3]*g(i,i);
	    }
	}
    }
  
}



//\eta/taupi/u^tau <\nabla^i u^\alpha>
void termb(double *a,int i,int alpha,int sx,int sy)
{
  a[0]=0;
  a[1]=0;
  a[2]=0;
  a[3]=0;
  
  double grxux[4];
  double grxuy[4];
  double gryuy[4];
  double grxut[4];
  double gryut[4];
 
  int alphap,ip;

  alphap=alpha;
  ip=i;

  gradxux(grxux,sx,sy);
  gradxuy(grxuy,sx,sy);
  gradyuy(gryuy,sx,sy);
  gradxut(grxut,sx,sy);
  gradyut(gryut,sx,sy);


  int phi=0;
  if(alphap<ip)
    {
      phi=alphap;
      alphap=ip;
      ip=phi;
    }

  for(int zaehler=0;zaehler<4;zaehler++)
    {    
      
      if(ip==0 && alphap==0)
	{
	  a[zaehler]+=etataupi(sx,sy)*grxux[zaehler]/globut;
	  //if (zaehler==0&&(sx==Middle+8)&&sy==Middle)
	  //  cout << "here tb" << a[zaehler]/e[sx][sy] << endl;
	}
      if(ip==0 && alphap==1)
	{
	  a[zaehler]+=etataupi(sx,sy)*grxuy[zaehler]/globut;
	}
      if(ip==1 && alphap==1)
	{
	  a[zaehler]+=etataupi(sx,sy)*gryuy[zaehler]/globut;
	}
      if(ip==0 && alphap==2)
	{
	  a[zaehler]-=etataupi(sx,sy)*grxut[zaehler]/globut;
	}
      if(ip==1 && alphap==2)
	{
	  a[zaehler]-=etataupi(sx,sy)*gryut[zaehler]/globut;
	}
      //printf("gr %f\n", grxux[3]);    
    }
  
  //  if ((sx==Middle+8)&&(sy==Middle)&&i==0&&alpha==0)
  // cout << "terb " << grxux[3]/globut  << endl;

}

//term c
//1/(tau_Pi u^tau) Pi^{i alpha}
double termc(int i,int alpha,int sx,int sy)
{
  double temp=0;
  temp+=pishell(i,alpha,sx,sy)/(taupi(sx,sy)*globut);
  //printf("tau %f pi %f\n",taupi(sx,sy),pi(i,alpha,sx,sy));
return temp;
}


//term d
//1/u^tau (u^i Pi^alpha_kappa+u^alpha Pi^i_kappa) D u^kappa
void  termd(double *a,int i,int alpha,int sx,int sy)
{
  a[0]=0;
  a[1]=0;
  a[2]=0;
  a[3]=0;
  double Du[4];
	
  for (int kappa=0;kappa<=3;kappa++)
    {
      Dumu(Du,kappa,sx,sy);
      for (int beta=0;beta<=3;beta++)
	for (int zaehler=0;zaehler<4;zaehler++)
	  a[zaehler]+=1/globut*(umu(i,sx,sy)*gdown(kappa,beta)*pishell(alpha,beta,sx,sy)+umu(alpha,sx,sy)*gdown(kappa,beta)*pishell(i,beta,sx,sy))*Du[zaehler];
    }
}

//term e
//vanishes
/*double terme(int i,int alpha,int sx,int sy)
{
  double temp=0;
  for (int kappa=0;kappa<=3;kappa++)
    for (int beta=0;beta<=3;beta++)
      temp+=umu(kappa,sx,sy)/globut*(gamma(i,kappa,beta)*pi(beta,alpha,sx,sy)+gamma(alpha,kappa,beta)*pi(i,beta,sx,sy));
  return temp;
  }*/



//u^j/u^tau \partial_j Pi^{i alpha}
double termf(int i,int alpha,int sx,int sy)
{
  double temp=0;
  for(int j=0;j<=1;j++)
    temp+=umu(j,sx,sy)/globut*djpi(j,i,alpha,sx,sy);
  return temp;
}

//shear-shear self coupling
double termg(int i,int alpha, int sx,int sy)
{
  double temp=0;

  for(int j=0;j<=3;j++)
    for (int k=0;k<=3;k++)
      temp+=pishell(k,j,sx,sy)*pishell(k,j,sx,sy)*gdown(j,j)*gdown(k,k);
  
  temp*=Delta(i,alpha,sx,sy);
  temp/=-3.;

  for(int j=0;j<=3;j++)
    temp+=pishell(i,j,sx,sy)*pishell(alpha,j,sx,sy)*gdown(j,j);
    
  temp*=L1COEF;
  temp/=4*M_PI*ETAOS*taupi(sx,sy)*(e[sx][sy]+eos(e[sx][sy]));
				   
  return temp;
}

//d_t pi^{i,alpha} summary
void dtpiialpha(double *a,int i,int alpha,int sx,int sy)
{
  double terd[4];
  termd(terd,i,alpha,sx,sy);
  double terb[4];
  termb(terb,i,alpha,sx,sy);
  
  double tera[4];
  terma(tera,i,alpha,sx,sy);

  //  tera[0]=0.0;
  // tera[1]=0.0;
  //tera[2]=0.0;
  //tera[3]=0.0;
  

  a[0]=terb[0]-terd[0]+tera[0];
  a[1]=terb[1]-terd[1]+tera[1];
  a[2]=terb[2]-terd[2]+tera[2];
  a[3]=terb[3]-termf(i,alpha,sx,sy)-terd[3]-termc(i,alpha,sx,sy)+tera[3];
  a[3]-=termg(i,alpha,sx,sy);

   if (isnan(a[3])!=0)
    {
      cout << "Problem in dtpiialpha" << endl;
      cout << "More specific: " << terb[3] << "\t" << terd[3] << "\t";
      cout << termf(i,alpha,sx,sy) << "\t" << termc(i,alpha,sx,sy) << endl;
    }

  /*
    if ((sx==Middle+8)&&(sy==Middle)&&(i==0)&&(alpha==0))
    {
    //cout << "tb " << terb[0]/e[sx][sy] << " terd " << terd[0]/e[sx][sy] << endl;
    cout << "dtpi " << a[3]/e[sx][sy];
    cout << "\t where " << -termc(i,alpha,sx,sy)/e[sx][sy];
    cout << "\t and " << terb[3]/e[sx][sy];
    cout << endl;
    }
  */
  //a[3]=-termf(i,alpha,sx,sy)/*-terme(i,alpha,sx,sy)*/-terd[3]-termc(i,alpha,sx,sy)+terb[3]/*-terma(i,alpha,sx,sy)*/;
  //printf("i %i a %i tf %f te %f td %f tc %f tb %f ux %f uy %f\n",i,alpha,termf(i,alpha,sx,sy),terme(i,alpha,sx,sy),terd[3],termc(i,alpha,sx,sy),terb[3],u[0][sx][sy],u[1][sx][sy]);
}

//\partial_\tau \Pi^{i \alpha}
void dtpishell(double *a,int i,int alpha,int sx,int sy)
{
  int phi;
  if (alpha<i)
    {
      phi=alpha;
      alpha=i;
      i=phi;
    }
  if ((i==0)&&(alpha==0))
    {
      a[0]=dtpixx[0];
      a[1]=dtpixx[1];
      a[2]=dtpixx[2];
      a[3]=dtpixx[3];
    }
  if ((i==0)&&(alpha==1))
    {
      a[0]=dtpixy[0];
      a[1]=dtpixy[1];
      a[2]=dtpixy[2];
      a[3]=dtpixy[3];
    }
  if ((i==0)&&(alpha==2))
    dtpiialpha(a,i,alpha,sx,sy);
  if ((i==1)&&(alpha==1))
    {
      a[0]=dtpiyy[0];
      a[1]=dtpiyy[1];
      a[2]=dtpiyy[2];
      a[3]=dtpiyy[3];
    }
  if ((i==1)&&(alpha==2))
    dtpiialpha(a,i,alpha,sx,sy);

}

//d_t pi^alpha,t 
double help(int alpha,int sx,int sy)
{
  double temp=0;
  for(int i=0;i<=1;i++)
    temp+=pishell(i,alpha,sx,sy)*umu(i,sx,sy)/globut;
  return temp;
}

//u^j/u^tau \partial_\tau \Pi^{j \alpha}
void help2(double *a,int alpha,int sx,int sy)
{
  double dtpiialp[4];
  
  
  a[0]=0;
  a[1]=0;
  a[2]=0;
  a[3]=0;
  
  for(int j=0;j<2;j++)
    {
      dtpishell(dtpiialp,j,alpha,sx,sy);     
      for(int zaehler=0;zaehler<4;zaehler++)
	a[zaehler]+=umu(j,sx,sy)/globut*dtpiialp[zaehler];
	  //printf("dt %f\n",dtpiialp[3]);
    }

      // }
  
}

//\partial_tau \Pi^{\alpha tau} 
void dtpialphat(double *a,int alpha,int sx,int sy)
{
  double hel2[4];
  help2(hel2,alpha,sx,sy);

  //qcout << "a[1]" << a[1] << endl;

  a[0]=hel2[0]+pishell(0,alpha,sx,sy)/globut-help(alpha,sx,sy)*umu(0,sx,sy)/(globut*globut);
  a[1]=hel2[1]+pishell(1,alpha,sx,sy)/globut-help(alpha,sx,sy)*umu(1,sx,sy)/(globut*globut);
  a[2]=0;
  a[3]=hel2[3];
  

}

//First part of third term

//g_{alpha kappa} \Delta^{\mu \kappa} \partial_\beta \Pi^{\alpha \beta}

void gDeltapi(double *a,int mu,int sx,int sy)
{
  double temp0=0;
  double temp1=0;
  double temp3=0;
  double dtpialpt[4];
  for(int alpha=0;alpha<=3;alpha++)
    {     
      dtpialphat(dtpialpt,alpha,sx,sy);
      for(int kappa=0;kappa<=3;kappa++) 
	{
	  temp0+=gdown(alpha,kappa)*Delta(mu,kappa,sx,sy)*dtpialpt[0];
	  temp1+=gdown(alpha,kappa)*Delta(mu,kappa,sx,sy)*dtpialpt[1];
	  temp3+=gdown(alpha,kappa)*Delta(mu,kappa,sx,sy)*(dtpialpt[3]+djpi(0,alpha,0,sx,sy)+djpi(1,alpha,1,sx,sy));
	}
    }

  a[0]=temp0;
  a[1]=temp1;
  a[2]=0;
  a[3]=temp3;
  
}



//Eq. (2)

//fourth term De

void De(double *a,int sx,int sy)
{
  a[0]=0;
  a[1]=0;
  a[2]=globut;
  a[3]=umu(0,sx,sy)*dxe(sx,sy)+umu(1,sx,sy)*dye(sx,sy);
}


//fifth term (e+p)Nabla_mu u^Mu

//Nabla_mu u^mu
void nablamuumu(double *a,int sx,int sy)
{
  a[0]=thf[0];//-Delta(0,2,sx,sy)+Delta(2,2,sx,sy)*umu(0,sx,sy)/globut;
  a[1]=thf[1];//-Delta(1,2,sx,sy)+Delta(2,2,sx,sy)*umu(1,sx,sy)/globut;
  a[2]=0;
  a[3]=thf[3];//-Delta(0,0,sx,sy)*dxu(0,sx,sy)-Delta(0,1,sx,sy)*dyu(0,sx,sy)-Delta(1,0,sx,sy)*dxu(1,sx,sy)-Delta(1,1,sx,sy)*dyu(1,sx,sy)+Delta(2,0,sx,sy)*dxu(2,sx,sy)+Delta(2,1,sx,sy)*dyu(2,sx,sy);
}

//sixth term 1/2 pi^mu,nu <grad_mu u_nu>
void pigradu(double *a,int sx,int sy)
{
  a[0]=0;
  a[1]=0;
  a[2]=0;
  a[3]=0;
  
  double grxux[4];
  double grxuy[4];
  double gryuy[4];
  double grxut[4];
  double gryut[4];
  double grtut[4];
  double greue[4];

  gradxux(grxux,sx,sy);
  gradxuy(grxuy,sx,sy);
  gradyuy(gryuy,sx,sy);
  gradxut(grxut,sx,sy);
  gradyut(gryut,sx,sy);
  gradtut(grtut,sx,sy);
  gradeue(greue,sx,sy);

  int mup,nup;
  
  for(int zaehler=0;zaehler<=3;zaehler++)
    for(int mu=0;mu<=3;mu++)
        for(int nu=0;nu<=3;nu++)
	  { 
	  mup=mu;
	  nup=nu;
	  
	  int phi=0;
	  if(mu>nu)
	    {
	      phi=mup;
	      mup=nup;
	      nup=phi;
	    }
	  
	  if(mup==0 && nup==0)
	    {
	      a[zaehler]+=1/2.*pishell(mup,nup,sx,sy)*grxux[zaehler];
	    }
	  
	  if(mup==0 && nup==1)
	    {
	      a[zaehler]+=1/2.*pishell(mup,nup,sx,sy)*grxuy[zaehler];
	    }
	  if(mup==0 && nup==2)
	    {
	      a[zaehler]+=1/2.*pishell(mup,nup,sx,sy)*grxut[zaehler];
	    }
	  if(mup==1 && nup==1)
	    {
	      a[zaehler]+=1/2.*pishell(mup,nup,sx,sy)*gryuy[zaehler];
	    }
	  if(mup==1 && nup==2)
	    {
	      a[zaehler]+=1/2.*pishell(mup,nup,sx,sy)*gryut[zaehler];
	    }
	  if(mup==2 && nup==2)
	    {
	      a[zaehler]+=1/2.*pishell(mup,nup,sx,sy)*grtut[zaehler];
	    }
	  if(mup==3 && nup==3)
	    {
	      a[zaehler]+=1/2.*pishell(mup,nup,sx,sy)*greue[zaehler];
	    }
	  }
}


//gauss-jordan elimination procedure
//provides 
//void gaussj(double **a, int n,double **b, int m)
#include "GJE.cpp"
#include "diags.cpp"



//main update routine
inline void doInc(double eps) 
{
  double Dumx[4],Dumy[4],Nabpx[4],Nabpy[4],gDeltpix[4],gDeltpiy[4],Des[4],nablau[4],pigru[4];

  //cout << "am here " << endl;

//   time_t time1,time2;

  int debug=0;

  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	//these are here to fix the internal indices once
	//instead of having to recalc every time -- faster!
	//globali=geti(e[sx][sy]);
	//globalx=getx(globali,e[sx][sy]);

	globut=ut(sx,sy);

	thf[0]=theta0(sx,sy);
	thf[1]=theta1(sx,sy);
	thf[3]=theta3(sx,sy);

	mypixt=pi(0,2,sx,sy);
	mypiyt=pi(1,2,sx,sy);
	mypitt=pi(2,2,sx,sy);
	mypiee=pi(3,3,sx,sy);

	dtpiialpha(dtpixx,0,0,sx,sy);
	dtpiialpha(dtpixy,0,1,sx,sy);
	dtpiialpha(dtpiyy,1,1,sx,sy);

	

	Dumu(Dumx,0,sx,sy); 
	Dumu(Dumy,1,sx,sy);
	Nablap(Nabpx,0,sx,sy);
	Nablap(Nabpy,1,sx,sy);
	gDeltapi(gDeltpix,0,sx,sy);
	gDeltapi(gDeltpiy,1,sx,sy);
	De(Des,sx,sy);	
	nablamuumu(nablau,sx,sy);
	pigradu(pigru,sx,sy);
	
	


	//Matrix
	
	/*	dtmat[0][0]=e[sx][sy]*(1+cs2(sx,sy))*Dumx[0]-Nabpx[0];//+gDeltagamma(0,sx,sy)+gDeltpix[0];
	  dtmat[0][1]=e[sx][sy]*(1+cs2(sx,sy))*Dumx[1]-Nabpx[1];//+//gDeltagamma(0,sx,sy)+gDeltpix[1];
	  dtmat[0][2]=e[sx][sy]*(1+cs2(sx,sy))*Dumx[2]-Nabpx[2];//+gDeltagamma(0,sx,sy)+gDeltpix[2];
	  dtmat[1][0]=e[sx][sy]*(1+cs2(sx,sy))*Dumy[0]-Nabpy[0];//+gDeltagamma(1,sx,sy)+gDeltpiy[0];
	  dtmat[1][1]=e[sx][sy]*(1+cs2(sx,sy))*Dumy[1]-Nabpy[1];//+gDeltagamma(1,sx,sy)+gDeltpiy[1];
	  dtmat[1][2]=e[sx][sy]*(1+cs2(sx,sy))*Dumy[2]-Nabpy[2];//+gDeltagamma(1,sx,sy)+gDeltpiy[2];
	  dtmat[2][0]=Des[0]+e[sx][sy]*(1+cs2(sx,sy))*nablau[0];//-pigru[0];
	  dtmat[2][1]=Des[1]+e[sx][sy]*(1+cs2(sx,sy))*nablau[1];//-pigru[1];
	  dtmat[2][2]=Des[2]+e[sx][sy]*(1+cs2(sx,sy))*nablau[2];//-pigru[2];
	*/
	
	dtmat[0][0]=(e[sx][sy]+eos(e[sx][sy]))*Dumx[0]-Nabpx[0]
	  //+gDeltagamma(0,sx,sy)
	  +gDeltpix[0];
	dtmat[0][1]=(e[sx][sy]+eos(e[sx][sy]))*Dumx[1]-Nabpx[1]
	  //+gDeltagamma(0,sx,sy)
	  +gDeltpix[1];
	dtmat[0][2]=(e[sx][sy]+eos(e[sx][sy]))*Dumx[2]-Nabpx[2]//+gDeltagamma(0,sx,sy)
	  +gDeltpix[2];
	dtmat[1][0]=(e[sx][sy]+eos(e[sx][sy]))*Dumy[0]-Nabpy[0]//+gDeltagamma(1,sx,sy)
	  +gDeltpiy[0];
	dtmat[1][1]=(e[sx][sy]+eos(e[sx][sy]))*Dumy[1]-Nabpy[1]//+gDeltagamma(1,sx,sy)
	  +gDeltpiy[1];
	dtmat[1][2]=(e[sx][sy]+eos(e[sx][sy]))*Dumy[2]-Nabpy[2]//+gDeltagamma(1,sx,sy)
	  +gDeltpiy[2];
	dtmat[2][0]=Des[0]+(e[sx][sy]+eos(e[sx][sy]))*nablau[0]-pigru[0];
	dtmat[2][1]=Des[1]+(e[sx][sy]+eos(e[sx][sy]))*nablau[1]-pigru[1];
	dtmat[2][2]=Des[2]+(e[sx][sy]+eos(e[sx][sy]))*nablau[2]-pigru[2];
	

	//Vector (written like a 3*0 Matrix)
	vec[0][0]=Nabpx[3]-(e[sx][sy]+eos(e[sx][sy]))*Dumx[3]-gDeltagamma(0,sx,sy)-gDeltpix[3];

	if (isnan(vec[0][0])!=0)
	  {
	    cout << "Problem in v00" << endl;
	    cout << "More specific 1 : " << Nabpx[3] << "\t";
	    cout << "2 : " << (e[sx][sy]+eos(e[sx][sy]))*Dumx[3] << "\t";
	    cout << "3 : " << gDeltagamma(0,sx,sy) << "\t";
	    cout << "4 : " << gDeltpix[3] << endl;
	  }
	vec[1][0]=Nabpy[3]-(e[sx][sy]+eos(e[sx][sy]))*Dumy[3]-gDeltagamma(1,sx,sy)-gDeltpiy[3];
	vec[2][0]=(Des[3]+(e[sx][sy]+eos(e[sx][sy]))*nablau[3])*(-1.0)+pigru[3];
	
	int check;
	check=gaussj(dtmat,3,vec,1);
		//printf("cs %f\n",cs2(sx,sy));  
	
	if (check==0)
	  {
	    for (int i=0;i<2;i++)
	      U[i][sx][sy]=u[i][sx][sy]+eps*vec[i][0];
	    //printf("vecx %f   vecy %f\n",vec[0][0],vec[1][0]);
	    

	    E[sx][sy]=e[sx][sy]+eps*vec[2][0];
	    
	    
	    
	    

	    Pixx[sx][sy]=pixx[sx][sy]+eps*(dtpixx[0]*vec[0][0]+dtpixx[1]*vec[1][0]+dtpixx[2]*vec[2][0]+dtpixx[3]);
	    
	    
	    Pixy[sx][sy]=pixy[sx][sy]+eps*(dtpixy[0]*vec[0][0]+dtpixy[1]*vec[1][0]+dtpixy[2]*vec[2][0]+dtpixy[3]);
	    
	    Piyy[sx][sy]=piyy[sx][sy]+eps*(dtpiyy[0]*vec[0][0]+dtpiyy[1]*vec[1][0]+dtpiyy[2]*vec[2][0]+dtpiyy[3]);
	    
	    //Pixx[sx][sy]=0.0;
	    //Pixy[sx][sy]=0.0;
	    //Piyy[sx][sy]=0.0;
	  }
	else
	  {
	    cout << "Error! at " << sx  << "  " << sy << endl;
	    snapshot(t);
	    bflag=1;
	  }
      }

  
}


//main driver routine
void Evolve() 
{
  //setting step sizes to maximum step size
  double eps = EPS;
  long int i=0;
  //for (long int i=1;i<=STEPS;i++) 
  //{
  while((reachedTf==0)&&(wflag==0)) 
    {
      i++;
      // evolve fields eps forward in time storing updated fields in captial vars
	  
	  doInc(eps); 
      
	  // measurements and data dump
	  //if ( (i>1 && (i-1)%UPDATE==0)) 
	  //{
	      if ((i-1)%UPDATE==0) 
		{
		  outputMeasurements(t);
		  
		  if (FREEZE > 1) 
		  {
		    //Copy fields to memory to compare with next UPDATE time step
		    copyUPDATE();
		  }
		  
		}
	      if ((i-1)%SNAPUPDATE==0) 
		{      
		  snapshot(t); 
		}
	  
	      //}
	  
	  //copy fields from capital vars to lowercase vars
	  copyDown();
	  
	  // increment time
	  t += eps;
					       
	  if (bflag==1)
	    break;

	  }
}




void printDivider() 
{
  int dwidth = 13;
  for (int i=0;i<8*dwidth;i++) cout << "-"; cout << endl;
  return;
}

void cleanup()
{
  gsl_spline_free (workspline);
  gsl_interp_accel_free (wac);
  gsl_spline_free (pspline);
  gsl_interp_accel_free (pacc);
  gsl_spline_free (Tspline);
  gsl_interp_accel_free (Tacc);
  gsl_spline_free (cs2spline);
  gsl_interp_accel_free (cs2acc);
}



int main() 
{
  
  extern void readParameters(const char*);
  
  printDivider();
  
  readParameters("data/params.txt");
  
  
  printDivider();

  allocateMemory();
 
  setInitialConditions();

  printf("Programmstart\n");  

  freeze_out.open("data/freezeout.dat", ios::out);
  meta.open("data/meta.dat", ios::out);
  ecces.open("data/ecc.dat", ios::out);

  
  Evolve();
 
  freeze_out.close();
  meta.close();
  ecces.close();

  cleanup();

  if (bflag==0)
    cout << "Done.\n";
  else
    cout << "Aborted because encountered nan.\n";


  return 0;
}
