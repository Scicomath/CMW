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
double B=0.0,AT=0.05,EPS=0.001;
double ETAOS=0.3;
double TSTART=0.5,TF=0.1,TINIT=1.0;
double L1COEF;
double L2COEF;
double IC=0;
int PTASIZE,PHIPASIZE;
int FREEZE=0;
double PTMAX, TRIEPS = 0.0, TRIANGLE = 0.0;
double QUADEPS = 0.0, QUADANGLE = 0.0;
double QUINTEPS = 0.0, QUINTANGLE = 0.0;
double SEXEPS = 0.0, SEXANGLE = 0.0;
double SEPTEPS = 0.0, SEPTANGLE = 0.0;
double BIEPS = 0.0, BIANGLE = 0.0;

//controls value of tau_Pi
double COEFF=3.0;

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
//double Rnuc = 6.6; //Lead
double Rnuc=6.4; //Gold

//wood-saxon parameter in fm;
//double anuc = 0.55; //Lead
double anuc=0.54;  //Gold

//inelastic nucleon-nucleon cross-section in mb
//double sigmaNN = 60; //LHC (5.5 A*TeV)
double sigmaNN = 40; //RHIC (200 A*GeV)

//normalization of thickness function
// double TAnorm = 2*208/1286.8; //Lead
double TAnorm = 2*197./1175.22; //Gold



//for the equation of state
double *eoT4,*cs2i,*poT4,*Ti;

// these hold the current values
double ***u,**e,**pixy,**pixx,**piyy;

double **cyminp;

double atuomas=Rnuc/125;
double g2mua=0.43328;
double ttuomas;

//overall time
double t = 0;

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
fstream cym,inited;

//splines -- for fancy freeze-out
gsl_interp_accel *wac;
gsl_spline *workspline;
//splines -- for equation of state
gsl_interp_accel *pacc,*Tacc,*cs2acc;
gsl_spline *pspline,*Tspline,*cs2spline;

//to know where to stop interpolation
double lowestE;

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

 //allocating memory
  cyminp = new double*[500];
  for (int i=0;i<500;i++) 
    cyminp[i] = new double[500];


}


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
//(x, y in fm)
double TA(double x, double y)
{
  double temp=0;
  for (int i=0;i<1200;i++)
    temp+=WS(x,y,(i+0.5)*Rnuc/400.)*Rnuc/400.;

  //return result normalized
  return temp*TAnorm;
  
 //return result normalized to a gold nucleus (A=197)
  //return 2*temp*197./1175.22;

 //return result normalized to a lead nucleus (A=208)
  //return 2*temp*208/1286.8;
}


//number density of wounded nucleons
//from Kolb et. al, hep-ph/0103234
double getwnuc(double xx,double yy,double b)
{
  
  double mTAp=TA(xx+b/2.,yy);
  double mTAm=TA(xx-b/2.,yy);

  //return mTA*2*(1.-pow(1-mTA/197.*4.,197.));
  double temp=0;
  temp+=mTAp*(1.-exp(-mTAm*sigmaNN*0.1));
  temp+=mTAm*(1.-exp(-mTAp*sigmaNN*0.1));
  return temp;
}


//number density of binary collisions
//from Kolb et. al, hep-ph/0103234
double getbin(double xx,double yy,double b)
{
  
  double mTAp=TA(xx+b/2.,yy);
  double mTAm=TA(xx-b/2.,yy);

  //return mTA*2*(1.-pow(1-mTA/197.*4.,197.));
  double temp=0;
  
  temp=4*mTAp*mTAm;

  return temp;
}

double getscal(double ss)
{
  long int i;
  for (i=0;i<length;i++)
    {
      if ((Ti[i]*Ti[i]*Ti[i]*(eoT4[i]+poT4[i]))>ss)
	break;
    }

  double temp=eoT4[i];
  temp*=Ti[i]*AT;
  temp*=Ti[i]*AT;
  temp*=Ti[i]*AT;
  temp*=Ti[i]*AT;

  return temp;
}

//alias Npart
double Nwounded(double b)
{
  int DIVIDER=10;
  double temp=0;
  for (double xx=-2*Rnuc;xx<2*Rnuc;xx+=Rnuc/DIVIDER)
    for (double yy=-2*Rnuc;yy<2*Rnuc;yy+=Rnuc/DIVIDER)
      {
	temp+=getwnuc(xx,yy,b);
      }
  return temp*Rnuc/DIVIDER*Rnuc/DIVIDER;
}

double Ncoll(double b)
{
  int DIVIDER=10;
  double temp=0;
  for (double xx=-2*Rnuc;xx<2*Rnuc;xx+=Rnuc/DIVIDER)
    for (double yy=-2*Rnuc;yy<2*Rnuc;yy+=Rnuc/DIVIDER)
      {
	temp+=getbin(xx,yy,b);
      }
  return temp*Rnuc/DIVIDER*Rnuc/DIVIDER;
}


double aniso()
{
  double diff=0,sum=0;
  
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	diff+=e[sx][sy]*((sy-Middle)*(sy-Middle)-(sx-Middle)*(sx-Middle));
	sum+=e[sx][sy]*((sx-Middle)*(sx-Middle)+(sy-Middle)*(sy-Middle));
      }

  return diff/sum;
}


//  calculate eccentricities epsilon_n and angles psi_n
void eps()
{
  //variables to calculate epsilon1 through epsilon7 with energy density weighting
  double num2 = 0.0, numeps21 = 0.0, numeps22 = 0.0, den = 0.0;
  double den3 = 0.0, den4 = 0.0, den5 = 0.0, den6 = 0.0, den7 = 0.0;
  double numeps31 = 0.0, numeps32 = 0.0,numeps41 = 0.0, numeps42 = 0.0;
  double numeps51 = 0.0, numeps52 = 0.0,numeps61 = 0.0, numeps62 = 0.0;
  double numeps71 = 0.0, numeps72 = 0.0, numeps11 = 0.0, numeps12 = 0.0;
  double r2numeps31 = 0.0, r2numeps32 = 0.0,r2numeps41 = 0.0, r2numeps42 = 0.0;
  double r2numeps51 = 0.0, r2numeps52 = 0.0,r2numeps61 = 0.0, r2numeps62 = 0.0;
  double r2numeps71 = 0.0, r2numeps72 = 0.0;
  //variables to calculate epsilon1 through epsilon7 with entropy density weighting
  double snum2 = 0.0, snumeps21 = 0.0, snumeps22 = 0.0, sden = 0.0, stot = 0.0;
  double sden3 = 0.0, sden4 = 0.0, sden5 = 0.0, sden6 = 0.0, sden7 = 0.0;
  double snumeps31 = 0.0, snumeps32 = 0.0, snumeps41 = 0.0, snumeps42 = 0.0;
  double snumeps51 = 0.0, snumeps52 = 0.0, snumeps61 = 0.0, snumeps62 = 0.0;
  double snumeps71 = 0.0, snumeps72 = 0.0, snumeps11 = 0.0, snumeps12 = 0.0;
  double r2snumeps31 = 0.0, r2snumeps32 = 0.0, r2snumeps41 = 0.0, r2snumeps42 = 0.0;
  double r2snumeps51 = 0.0, r2snumeps52 = 0.0, r2snumeps61 = 0.0, r2snumeps62 = 0.0;
  double r2snumeps71 = 0.0, r2snumeps72 = 0.0;
  
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	double x = (sx-Middle)*AT/5.06842; //in fm
	double y = (sy-Middle)*AT/5.06842; //in fm
	double phi = atan2(y,x);
	double r2 = x*x + y*y;
	double r3 = pow(r2,1.5);
	double r4 = r2*r2;
	double r5 = r2*r3;
	double r6 = r2*r2*r2;
	double r7 = r4*r3;
	//  Calculate eccentricities epsilon1 to epsilon7 
	//  with energy density weighting
	
	//  denominators:
	den += r2*e[sx][sy];
	den3 += r3*e[sx][sy];
	den4 += r4*e[sx][sy];
	den5 += r5*e[sx][sy];
	den6 += r6*e[sx][sy];
	den7 += r7*e[sx][sy];
	//  factors for the numerator
	num2 += (y*y - x*x)*e[sx][sy];// standard eccentricity
	numeps11 += r3*cos(phi)*e[sx][sy];
	numeps12 += r3*sin(phi)*e[sx][sy];
	numeps21 += r2*cos(2*phi)*e[sx][sy];// for participant eccentricity
	numeps22 += r2*sin(2*phi)*e[sx][sy];// for psrticipant eccentricity
	numeps31 += r3*cos(3*phi)*e[sx][sy];// for triangularity with r^n weighting
	numeps32 += r3*sin(3*phi)*e[sx][sy];// for triangularity with r^n weighting
	numeps41 += r4*cos(4*phi)*e[sx][sy];
	numeps42 += r4*sin(4*phi)*e[sx][sy];
	numeps51 += r5*cos(5*phi)*e[sx][sy];
	numeps52 += r5*sin(5*phi)*e[sx][sy];
	numeps61 += r6*cos(6*phi)*e[sx][sy];
	numeps62 += r6*sin(6*phi)*e[sx][sy];
	numeps71 += r7*cos(7*phi)*e[sx][sy];
	numeps72 += r7*sin(7*phi)*e[sx][sy];
	// r^2 weighting
	r2numeps31 += r2*cos(3*phi)*e[sx][sy];// for triangularity with r^2 weighting
	r2numeps32 += r2*sin(3*phi)*e[sx][sy];// for triangularity with r^2 weighting
	r2numeps41 += r2*cos(4*phi)*e[sx][sy];
	r2numeps42 += r2*sin(4*phi)*e[sx][sy];
	r2numeps51 += r2*cos(5*phi)*e[sx][sy];
	r2numeps52 += r2*sin(5*phi)*e[sx][sy];
	r2numeps61 += r2*cos(6*phi)*e[sx][sy];
	r2numeps62 += r2*sin(6*phi)*e[sx][sy];
	r2numeps71 += r2*cos(7*phi)*e[sx][sy];
	r2numeps72 += r2*sin(7*phi)*e[sx][sy];
	
	//  entropy density weighting calculations
	double s = 0.0;
	double p, T;
	if (e[sx][sy] < lowestE)
	{
		p = AT*AT*AT*AT*gsl_spline_eval(pspline,e[sx][sy],pacc);
		T = AT*gsl_spline_eval(Tspline,e[sx][sy],pacc);
	}
	else
	{
		p = AT*AT*AT*AT*e[sx][sy]*gsl_spline_eval(cs2spline,lowestE,cs2acc);
		T = sqrtl(sqrtl(e[sx][sy]/eoT4[0]));
	}
	s = (e[sx][sy] + p)/T;
	
	stot += s;//  total entropy
	sden += r2*s;
	sden3 += r3*s;
	sden4 += r4*s;
	sden5 += r5*s;
	sden6 += r6*s;
	sden7 += r7*s;
	snum2 += (y*y - x*x)*s;// standard eccentricity
	snumeps21 += r2*cos(2*phi)*s;// for participant eccentricity
	snumeps22 += r2*sin(2*phi)*s;// for participant eccentricity
	//  r^n weighting
	snumeps11 += r3*cos(phi)*s;
	snumeps12 += r3*sin(phi)*s;
	snumeps31 += r3*cos(3*phi)*s;
	snumeps32 += r3*sin(3*phi)*s;
	snumeps41 += r4*cos(4*phi)*s;
	snumeps42 += r4*sin(4*phi)*s;
	snumeps51 += r5*cos(5*phi)*s;
	snumeps52 += r5*sin(5*phi)*s;
	snumeps61 += r6*cos(6*phi)*s;
	snumeps62 += r6*sin(6*phi)*s;
	snumeps71 += r7*cos(7*phi)*s;
	snumeps72 += r7*sin(7*phi)*s;
	//  r^2 weighting
	r2snumeps31 += r2*cos(3*phi)*s;
	r2snumeps32 += r2*sin(3*phi)*s;
	r2snumeps41 += r2*cos(4*phi)*s;
	r2snumeps42 += r2*sin(4*phi)*s;
	r2snumeps51 += r2*cos(5*phi)*s;
	r2snumeps52 += r2*sin(5*phi)*s;
	r2snumeps61 += r2*cos(6*phi)*s;
	r2snumeps62 += r2*sin(6*phi)*s;
	r2snumeps71 += r2*cos(7*phi)*s;
	r2snumeps72 += r2*sin(7*phi)*s;
      }
    cout << "Energy Density Weighting + r^n weighting:\n";
    cout << "ecc = " << num2/den << endl;
    cout << "epsilon1 = " << sqrt(numeps11*numeps11 + numeps12*numeps12)/den3 << endl;
    cout << "psi1 = " << atan2(-numeps12,-numeps11) << endl;
    cout << "epsilon2 = " << sqrt(numeps21*numeps21 + numeps22*numeps22)/den << endl;
    cout << "psi2 = " << atan2(-numeps22,-numeps21)/2 << endl;
    cout << "epsilon3 = " << sqrt(numeps31*numeps31 + numeps32*numeps32)/den3 << endl;
    cout << "psi3 = " << atan2(-numeps32,-numeps31)/3 << endl;
    cout << "epsilon4 = " << sqrt(numeps41*numeps41 + numeps42*numeps42)/den4 << endl;
    cout << "psi4 = " << atan2(-numeps42,-numeps41)/4 << endl;
    cout << "epsilon5 = " << sqrt(numeps51*numeps51 + numeps52*numeps52)/den5 << endl;
    cout << "psi5 = " << atan2(-numeps52,-numeps51)/5 << endl;
    cout << "epsilon6 = " << sqrt(numeps61*numeps61 + numeps62*numeps62)/den6 << endl;
    cout << "psi6 = " << atan2(-numeps62,-numeps61)/6 << endl;
    cout << "epsilon7 = " << sqrt(numeps71*numeps71 + numeps72*numeps72)/den7 << endl;
    cout << "psi7 = " << atan2(-numeps72,-numeps71)/7 << endl;
    
    cout << "Entropy Density Weighting + r^n weighting:\n";
    cout << "ecc = " << snum2/sden << endl;
    cout << "epsilon1 = " << sqrt(snumeps11*snumeps11 + snumeps12*snumeps12)/den3 << endl;
    cout << "psi1 = " << atan2(-snumeps12,-snumeps11) << endl;
    cout << "epsilon2 = " << sqrt(snumeps21*snumeps21 + snumeps22*snumeps22)/sden << endl;
    cout << "psi2 = " << atan2(-snumeps22,-snumeps21)/2 << endl;
    cout << "epsilon3 = " << sqrt(snumeps31*snumeps31 + snumeps32*snumeps32)/sden3 << endl;
    cout << "psi3 = " << atan2(-snumeps32,-snumeps31)/3 << endl;
    cout << "epsilon4 = " << sqrt(snumeps41*snumeps41 + snumeps42*snumeps42)/sden4 << endl;
    cout << "psi4 = " << atan2(-snumeps42,-snumeps41)/4 << endl;
    cout << "epsilon5 = " << sqrt(snumeps51*snumeps51 + snumeps52*snumeps52)/sden5 << endl;
    cout << "psi5 = " << atan2(-snumeps52,-snumeps51)/5 << endl;
    cout << "epsilon6 = " << sqrt(snumeps61*snumeps61 + snumeps62*snumeps62)/sden6 << endl;
    cout << "psi6 = " << atan2(-snumeps62,-snumeps61)/6 << endl;
    cout << "epsilon7 = " << sqrt(snumeps71*snumeps71 + snumeps72*snumeps72)/sden7 << endl;
    cout << "psi7 = " << atan2(-snumeps72,-snumeps71)/7 << endl;
    cout << "total entropy = " << stot << endl;
    
    cout << "Energy Density Weighting + r^2 weighting:\n";
    cout << "epsilon3 = " << sqrt(r2numeps31*r2numeps31 + r2numeps32*r2numeps32)/den << endl;
    cout << "psi3 = " << atan2(-r2numeps32,-r2numeps31)/3 << endl;
    cout << "epsilon4 = " << sqrt(r2numeps41*r2numeps41 + r2numeps42*r2numeps42)/den << endl;
    cout << "psi4 = " << atan2(-r2numeps42,-r2numeps41)/4 << endl;
    cout << "epsilon5 = " << sqrt(r2numeps51*r2numeps51 + r2numeps52*r2numeps52)/den << endl;
    cout << "psi5 = " << atan2(-r2numeps52,-r2numeps51)/5 << endl;
    cout << "epsilon6 = " << sqrt(r2numeps61*r2numeps61 + r2numeps62*r2numeps62)/den << endl;
    cout << "psi6 = " << atan2(-r2numeps62,-r2numeps61)/6 << endl;
    cout << "epsilon7 = " << sqrt(r2numeps71*r2numeps71 + r2numeps72*r2numeps72)/den << endl;
    cout << "psi7 = " << atan2(-r2numeps72,-r2numeps71)/7 << endl;
    
    cout << "Entropy Density Weighting + r^2 weighting:\n";
    cout << "epsilon3 = " << sqrt(r2snumeps31*r2snumeps31 + r2snumeps32*r2snumeps32)/sden << endl;
    cout << "psi3 = " << atan2(-r2snumeps32,-r2snumeps31)/3 << endl;
    cout << "epsilon4 = " << sqrt(r2snumeps41*r2snumeps41 + r2snumeps42*r2snumeps42)/sden << endl;
    cout << "psi4 = " << atan2(-r2snumeps42,-r2snumeps41)/4 << endl;
    cout << "epsilon5 = " << sqrt(r2snumeps51*r2snumeps51 + r2snumeps52*r2snumeps52)/sden << endl;
    cout << "psi5 = " << atan2(-r2snumeps52,-r2snumeps51)/5 << endl;
    cout << "epsilon6 = " << sqrt(r2snumeps61*r2snumeps61 + r2snumeps62*r2snumeps62)/sden << endl;
    cout << "psi6 = " << atan2(-r2snumeps62,-r2snumeps61)/6 << endl;
    cout << "epsilon7 = " << sqrt(r2snumeps71*r2snumeps71 + r2snumeps72*r2snumeps72)/sden << endl;
    cout << "psi7 = " << atan2(-r2snumeps72,-r2snumeps71)/7 << endl;

}

//in inverse (fm/c)^2
double dilution()
{
  double diff=0,sum=0;
  
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	diff+=e[sx][sy]*5.06842/AT*5.06842/AT;
	sum+=e[sx][sy]*((sx-Middle)*(sx-Middle)+(sy-Middle)*(sy-Middle));
      }

  diff*=2./3.;

  return diff/sum;
}

//definition according to Bhalerao, nucl-th/0508009v3
double overlapS()
{
  double mx=0,my=0,mean=0;
  
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	mx+=e[sx][sy]*(sx-Middle)*(sx-Middle);
	my+=e[sx][sy]*(sy-Middle)*(sy-Middle);
	mean+=e[sx][sy];
      }

  return 4*M_PI*sqrt(mx*my)/mean;
}


//get Tuomas' results and smooth them
void prepareCYM(double bb)
{

  char buffer[50];


  switch ((int) round(10*TINIT/atuomas*g2mua))
    {
    case 10:
      {
	printf("Loading data for TINIT=1g2mu\n");
	sprintf(buffer,"endens00229.dat");
	ttuomas=2.315*atuomas;
	break;
      }
    case 20:
      {
	printf("Loading data for TINIT=2g2mu\n");
	sprintf(buffer,"endens00458.dat");
	ttuomas=4.605*atuomas;
	break;
      }
    case 30:
      {
	printf("Loading data for TINIT=3g2mu\n");
	sprintf(buffer,"endens00687.dat");
	ttuomas=6.895*atuomas;
	break;
      }
    case 40:
      {
	printf("Loading data for TINIT=4g2mu\n");
	sprintf(buffer,"endens00916.dat");
	ttuomas=9.185*atuomas;
	break;
      }
    case 50:
      {
	printf("Loading data for t=5g2mu TINIT=%f\n",TINIT/atuomas*g2mua);
	sprintf(buffer,"endens01145.dat");
	ttuomas=11.475*atuomas;
	break;
      }
    default:
      {
	printf("No available data for initial time %f\n",TINIT);
	exit(0);
      }
    }

  char regu[3];

  sprintf(regu,"05");

  char filenam[128];

  printf("Rnuc=%f\n",Rnuc);

  switch ((int) round(10*bb/Rnuc))
    {
    case 0: 
      {
	printf("Loading data for B=0*Rnuc\n");
	sprintf(filenam,"CYM/b00L%s/%s",regu,buffer);
	cym.open(filenam, ios::in);
	break;
      }
    case 20: 
      {
	printf("Loading data for B=2*Rnuc\n");
	sprintf(filenam,"CYM/b02L%s/%s",regu,buffer);
	cym.open(filenam, ios::in);
	break;
      }
    case 40: 
      {
	printf("Loading data for B=4*Rnuc\n");
	sprintf(filenam,"CYM/b04L%s/%s",regu,buffer);
	cym.open(filenam, ios::in);
	break;
      }
    case 60: 
      {
	printf("Loading data for B=6*Rnuc\n");
	sprintf(filenam,"CYM/b06L%s/%s",regu,buffer);
	cym.open(filenam, ios::in);
	break;
      } 
    case 80: 
      {
	printf("Loading data for B=8*Rnuc\n");
	sprintf(filenam,"CYM/b08L%s/%s",regu,buffer);
	cym.open(filenam, ios::in);	
	break;
      }
    case 100: 
      {
	printf("Loading data for B=10*Rnuc\n");
	sprintf(filenam,"CYM/b10L%s/%s",regu,buffer);
	cym.open(filenam, ios::in);	
	break;
      }
    case 120: 
      {
	printf("Loading data for B=12*Rnuc\n");
	sprintf(filenam,"CYM/b12L%s/%s",regu,buffer);
	cym.open(filenam, ios::in);	
	break;
      }
    case 140: 
      {
	printf("Loading data for B=14*Rnuc\n");
	sprintf(filenam,"CYM/b14L%s/%s",regu,buffer);
	cym.open(filenam, ios::in);	
	break;
      }
    case 160: 
      {
	printf("Loading data for B=16*Rnuc\n");
	sprintf(filenam,"CYM/b16L%s/%s",regu,buffer);
	cym.open(filenam, ios::in);	
	break;
      }
    default:
      {
	printf("No available data for impact parameter %f\n",B);
	exit(0);
      }
    }

  
  
  double temp;
  for (int i=0;i<500;i++)
    for (int j=0;j<500;j++)
      {
	cym >> temp;
	cyminp[i][j]=fabs(temp);
      }

  cym.close();

}

double smoothcym1(int sx,int sy)
{
  double temp=0;
  temp+=cyminp[sx][sy];
  temp+=cyminp[500-sx-1][sy];
  temp+=cyminp[sx][500-sy-1];
  temp+=cyminp[500-sx-1][500-sy-1];
  temp*=0.25;
  return temp;
}


double smoothcym2(int sx,int sy)
{
  double temp=0;
  int SMOOTHNESS=20;

  for (int i=-SMOOTHNESS;i<=SMOOTHNESS;i++)
    for (int j=-SMOOTHNESS;j<=SMOOTHNESS;j++)
      {
	if ((sx+i>=0)&&(sx+i<500)&&(sy+j>=0)&&(sy+j<500))
	  temp+=smoothcym1(sx+i,sy+j);
	else
	  {
	    printf("Error in smoothing: out of bounds\n");
	    printf("tried sx=%i sy=%i\n",sx,sy);
	    exit(0);
	  }
      }

  temp/=(2*SMOOTHNESS+1)*(2*SMOOTHNESS+1);

  return temp;
}

int trickindex(int s)
{
  if(s<=0)
    return (s+250);
  else
    return (s+249);
}




double getcym(double x,double y)
{
  int cx,cy;

  

  if (x>=0)
    cx=(int) (x/atuomas);
  else
    cx=(int) (x/atuomas-1);
  
  if (y>=0)
    cy=(int) (y/atuomas);
  else
    cy=(int) (y/atuomas-1);

  //using bilinear interpolation:
  //(Numerical recipes)

  double t=(x-cx*atuomas)/atuomas;
  
  double u=(y-cy*atuomas)/atuomas;
  
  int tcx,tcy;

  tcx=trickindex(cx);
  tcy=trickindex(cy);

  double y1,y2,y3,y4;

  y1=smoothcym2(tcx,tcy);
  y2=smoothcym2(tcx+1,tcy);
  y3=smoothcym2(tcx+1,tcy+1);
  y4=smoothcym2(tcx,tcy+1);

   double result;
  result=(1-t)*(1-u)*y1+t*(1-u)*y2+t*u*y3+(1-t)*u*y4;

  return result;
  
}

void cycleb()
{

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

  //locate starting temperature in EOS
  long int i;
  for (i=0;i<length;i++)
    {
      if (Ti[i]>TSTART)
	break;
    }
  
  cout << "found temperature" << Ti[i] << endl;
  //energy density at that temperature in physical units
  double e0;

  e0=(eoT4[i])*Ti[i]*AT;
  e0*=Ti[i]*AT;
  e0*=Ti[i]*AT;
  e0*=Ti[i]*AT;

  //printf("Impact=%f\t Nwound=%f\t Npart=%f \n",B,Nwounded(B),Npart(B));

  

  double normcym;
  normcym=AT/5.06842;
  normcym*=AT/5.06842;
  normcym*=AT/5.06842;
  normcym*=AT/5.06842;
  //setting g^2=4
  
  for (int k=0;k<8;k++)
    { 
      double b=(double) 2*k;

      prepareCYM(b);
      normcym/=4*atuomas*atuomas*atuomas*ttuomas;


      for (int sx=1;sx<=NUMT;sx++)
	for (int sy=1;sy<=NUMT;sy++)
	  {
	    //wounded nucleon scaling
	    //e[sx][sy]=e0*getwnuc((sx-Middle)*AT/5.06842,(sy-Middle)*AT/5.06842,b)/4.29048;
	    //participant scaling
	    //e[sx][sy]=e0*getbin((sx-Middle)*AT/5.06842,(sy-Middle)*AT/5.06842,b)/18.4151;
	    e[sx][sy]=normcym*getcym((sx-Middle)*AT/5.06842,(sy-Middle)*AT/5.06842);
	  }
  

      //printf("aniso %f\n",aniso());
      
      printf("{%f,%f},\t",b,aniso());
    }

  printf("\n");
  
  /*
    printf("Number of wounded nucs:\n");
  for (double bb=0;bb<15;bb+=1)
    {
      double temp=Nwounded(bb);
      //printf("%f\t %f\n",bb,temp);
      printf("{%f,%f},\t",bb,temp);
    }
  printf("\n");

  printf("Number of participating nucs:\n");
  for (double bb=0;bb<15;bb+=1)
    {
      double temp=Nwounded(bb);
      //printf("%f\t %f\n",bb,temp);
      printf("{%f,%f},\t",bb,temp);
    }
  printf("\n");
  */

  wac=gsl_interp_accel_alloc (); 
  workspline=gsl_spline_alloc (gsl_interp_cspline, Middle);
 
}


void setInitialConditions()
{
  
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

  //locate starting temperature in EOS
  long int i;
  for (i=0;i<length;i++)
    {
      if (Ti[i]>TSTART)
	break;
    }
  
  cout << "found temperature" << Ti[i] << endl;
  //energy density at that temperature in physical units
  double e0;

  e0=(eoT4[i])*Ti[i]*AT;
  e0*=Ti[i]*AT;
  e0*=Ti[i]*AT;
  e0*=Ti[i]*AT;

  //entropy density at that temperature in physical units
  double s0;

  s0=(eoT4[i]+poT4[i]);
  s0*=Ti[i];
  s0*=Ti[i];
  s0*=Ti[i];
  
  //printf("Impact=%f\t Nwound=%f\t Npart=%f \n",B,Nwounded(B),Npart(B));

  double normcym;

  if (IC==-3)
    {
      prepareCYM(B);
      normcym=AT/5.06842;
      normcym*=AT/5.06842;
      normcym*=AT/5.06842;
      normcym*=AT/5.06842;
      //setting g^2=4
      normcym/=4*atuomas*atuomas*atuomas*ttuomas ;
    }

  if (IC<0)
    {
      if (IC==-1)
	printf("Glauber: ENpart scaling\n");
      if (IC==-2)
	printf("Glauber: ENcoll scaling\n");
      if (IC==-3)
	printf("CGC: Lappi's CYM\n");
      if (IC==-4)
	printf("Glauber: Dilution\n");
      if (IC==-5)
	printf("Woods-Saxon IC");
      if (IC==-6)
	printf("Triangular Gaussian\n");
      if (IC==-7)
	printf("Triangularly modified Glauber: ENcoll scaling\n");
      if (IC==-21)
	printf("Glauber: ENpart scaling -- for LHC\n");
      if (IC==-22)
	printf("Glauber: ENcoll scaling -- for LHC\n");
      if (IC==-27)
	printf("Triangularly modified Glauber: ENcoll scaling -- for LHC\n");
    }
  else
    {
      printf("Glauber Entropy scaling: Ncoll admixture %f percent\n",IC*100);
    }

  double b=B;
  double diluter=dilution()*2./3.;
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	if (IC<0)
	  {
	    //wounded nucleon scaling
	    if (IC==-1 || IC==-21)
	      e[sx][sy]=e0*getwnuc((sx-Middle)*AT/5.06842,(sy-Middle)*AT/5.06842,b)/4.29048;
	    //participant scaling
	    if (IC==-2 || IC==-22)
	      e[sx][sy]=e0*getbin((sx-Middle)*AT/5.06842,(sy-Middle)*AT/5.06842,b)/18.4151;
	    if (IC==-3)
	      e[sx][sy]=normcym*getcym((sx-Middle)*AT/5.06842,(sy-Middle)*AT/5.06842);
	    if (IC==-5)
	      {
		double x = (sx-Middle);
		double y = (sy-Middle);
		double scale = 5.0;
		double R = 50.0/scale;
		double a = 10.0/scale;
		double no = 0.25*scale;
		double nor = no*no*no*no;
	        e[sx][sy]=nor*(1+exp(-R/a))/(1+exp((sqrt(x*x+y*y)-R)/a));
	      }
	    if (IC==-6)
	    {
	      //Gaussien, stretched to various eccentricities (see arXiv:1007.5469)
		    double x = (sx-Middle)*AT/5.06842; //in fm
		    double y = (sy-Middle)*AT/5.06842; //in fm
		    double phi = atan2(y,x);
		    double Rgauss = 3.0; //in fm
		    double stretch = (1 + BIEPS*cos(2*(phi-BIANGLE)) + TRIEPS*cos(3*(phi-TRIANGLE))
		    			+ QUADEPS*cos(4*(phi-QUADANGLE)) + QUINTEPS*cos(5*(phi-QUINTANGLE))
		    			+ SEXEPS*cos(6*(phi-SEXANGLE)) + SEPTEPS*cos(7*(phi-SEPTANGLE)));
// 		    cout << "eps's\n" << BIEPS << "\n" << TRIEPS << "\n" << QUADEPS << "\n" << QUINTEPS
// 		    		<< "\n" << SEXEPS << endl << SEPTEPS << endl;
		    e[sx][sy] = e0*exp(-(x*x+y*y)*stretch/(2*Rgauss*Rgauss));
//		    double s0 = 0;
//		    if (e0 < lowestE)
//		    {
//			    double p0 = AT*AT*AT*AT*gsl_spline_eval(pspline,phys,pacc);
//			    s0 = (e0+p0)/TSTART;
//		    }
//		    else
//		    {
//			    double p0 = AT*AT*AT*AT*e0*gsl_spline_eval(cs2spline,lowestE,cs2acc);
//			    s0 = (e0+p0)/TSTART;
//		    }
//
//		    double x = (sx-Middle)*AT;
//		    double y = (sy-Middle)*AT;
//		    double phi = atan2(y,x);
//		    double s = s0*exp(-(x*x+y*y)/(2/Rnuc)*(1+0.3*cos(2*phi)));


	    }
	    if (IC==-7 || IC==-27)
	    {
	      //  energy density scaling with binary collisions, 
	      //  stretched to various eccentricities (see arXiv:1007.5469)
		    double x = (sx-Middle)*AT/5.06842;
		    double y = (sy-Middle)*AT/5.06842;
		    double phi = atan2(y,x);
// 		    double stretch = sqrt(1+TRIEPS*cos(3*(phi-TRIANGLE)))*sqrt(1+QUADEPS*cos(4*(phi-QUADANGLE)));
		    double stretch = sqrt(1 + BIEPS*cos(2*(phi-BIANGLE)) + TRIEPS*cos(3*(phi-TRIANGLE))
		    			+ QUADEPS*cos(4*(phi-QUADANGLE)) + QUINTEPS*cos(5*(phi-QUINTANGLE))
		    			+ SEXEPS*cos(6*(phi-SEXANGLE)) + SEPTEPS*cos(7*(phi-SEPTANGLE)));
		    e[sx][sy]=e0*getbin(x*stretch,y*stretch,b)/18.4151;
		    
	    }

	  }
	else
	  {
	    e[sx][sy]=getscal((getwnuc((sx-Middle)*AT/5.06842,(sy-Middle)*AT/5.06842,b)/4.29048*(1-IC)+getbin((sx-Middle)*AT/5.06842,(sy-Middle)*AT/5.06842,b)/18.4151*IC)*s0);
	  }
	
	
	
      }

  cout << "Npart " << Nwounded(0) << endl;

  cout << "ED at center " << e[Middle][Middle] << endl;


  cout << "overlap " << overlapS()*AT*AT << endl;

  cout << "ecc" << aniso() << endl;
  cout << "dilution" << dilution() << endl;

  wac=gsl_interp_accel_alloc (); 
  workspline=gsl_spline_alloc (gsl_interp_cspline, Middle);
  
 
}

void outputed()
{
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      inited << e[sx][sy] <<"\t";

  inited<<endl;
}


int main() 
{
  
  extern void readParameters(const char*);
  
  readParameters("data/params.txt");  
  
  allocateMemory();
  
  if (int(IC)/10 == -2)
  {
    cout << "Lead-Lead collision\n";
    Rnuc = 6.6; //Lead
    anuc = 0.55; //Lead
    sigmaNN = 60; //LHC (5.5 A*TeV)
    TAnorm = 2*208/1286.8; //Lead
  }
 
  setInitialConditions();
  
  //cycleb();

  if ((IC == -6) || (IC == -7) || (IC == -27))
  {
  //  calculate all epsilon_n and psi_n
  eps();
  }
  
  
  inited.open("data/inited.dat", ios::out);

  outputed();

  inited.close();
  
  printf("Done!\n");

  return 0;
}
