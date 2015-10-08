//
//
#include <iostream>
#include <iomanip> 
#include <cstdlib>
#include <string>
#include <algorithm>
#include <fstream>
#include "FTDipole.h"
#include "Random.h"
#include "PyRand.h"
#include "MCSaturation.h"

using namespace std;

// CGC test
void dndy();       // dN/dy with average Qs(b)
void dndyLocal();  // dN/dy with local Qs(x,y,b)
void bDepQs(double ecm);// impact parameter dependence of Ncoll, Npart, Qs(b)...
void dndyCentralityAtMidRap(double e,double lambda,double kglue,int m,int et);
void dndyLocal2(double ecm, double b, double lambda,double kglue,
       	int mode,int et, double hmass);
void CentralityBin(double ecm);
void viscinit(double ecm, double lambda, double kglue, int model,int Et);
void eccentricity(double ecm, double lambda, double kglue, int model,int Et);
void eccentricity2(double ecm, double lambda, double kglue, int model,int Et);
void eccentricity3(double ecm, double lambda, double kglue, int model,int Et);
void eccentricity4(double ecm, double lambda, double kglue, int model,int Et);
void MCAverageQs(double ecm);

int main(int argc, char *argv[]) {

    //double kglue = 0.7;
    //double kglue = 0.495;
    //double kglue = 1.23;     // fix alpha
    double kglue = 2.0;     // fix alpha
    //double kglue = 1.97;

    double lambda=0.288;

     double ecm=200;
//    double ecm=5500;
//    double b = 2.4;
    int Et = 0; int model=7;
    int nev=10; int cent = 1; double hmass = 0.0; double x=1e-2;

    // hydro parameters
    double tau0=0.6;
    //double kappa=0.057;  // CGC  alpha=1
    //double kappa=0.063;  // CGC  alpha=1
    double kappa=0.054;  // Npart  alpha=1


    for(int i=1; i<argc; i++) {
	if(!strcmp(argv[i],"-c")) cent = atoi(argv[i+1]);
	if(!strcmp(argv[i],"-n")) nev = atoi(argv[i+1]);
	if(!strcmp(argv[i],"-e")) ecm = atof(argv[i+1]);
	if(!strcmp(argv[i],"-m")) model = atoi(argv[i+1]);
	if(!strcmp(argv[i],"-x")) x = atof(argv[i+1]);
//	if(!strcmp(argv[i],"-b")) b = atof(argv[i+1]);
	if(!strcmp(argv[i],"-l")) lambda = atof(argv[i+1]);
	if(!strcmp(argv[i],"-k")) kglue = atof(argv[i+1]);
	if(!strcmp(argv[i],"-mass")) hmass = atof(argv[i+1]);
	if(!strcmp(argv[i],"-et")) Et = atoi(argv[i+1]);
	if(!strcmp(argv[i],"-ka")) kappa = atof(argv[i+1]);
	if(!strcmp(argv[i],"-tau0")) tau0 = atof(argv[i+1]);
    }

    // Set Random number generator.
    int randomSeed=13921;
    //int randomSeed=33929;
    Random* rand = new PyRand(randomSeed);
    Random::setRandom(rand);


    viscinit(ecm, lambda, kglue, model,Et);


   return 0;
}


#include "Regge96.h"
#include "OverLap.h"
#include "Saturation.h"
#include "KLNModel.h"
#include "CentralityCut.h"
#include "UnintegPartonDist.h"
#include "KLNfunc.h"
#include "KLNGamfunc.h"

// xxx everything between here and the next triple x added by Matt Luzum
// to calculate CGC initial conditions for viscous hydro code

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>

// these global vars are initialized from parameters file
// defaults set here are overridden by that file

int NUMT=8;
long int STEPS=4000,UPDATE=100,SNAPUPDATE=1000;
double B=0.0,AT=0.05,EPS=0.001;
double ETAOS=0.3;
double TSTART=0.5,TF=0.1,TINIT=1.0;
double IC=0;
int PTASIZE,PHIPASIZE;
int FREEZE=0;
double L1COEF=0.0;
double L2COEF=0.0;
double PTMAX=4.2, TRIEPS = 0.0, TRIANGLE, QUADEPS = 0.0, QUADANGLE = 0.0;
double BIEPS = 0.0, BIANGLE;
double QUINTEPS = 0.0, QUINTANGLE = 0.0, SEXEPS = 0.0, SEXANGLE = 0.0;
double SEPTEPS = 0.0, SEPTANGLE = 0.0;

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
double Rnuc=6.4;
//wood-saxon parameter in fm;
double anuc=0.54;


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
  
  double mx = 0.0, my = 0.0, mn = 0.0;


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

  lowestE=warry[0];

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

void viscinit(double ecm, double lambda, double kglue, int mode,int Et)
{
  // generates initial conditions for use with UVH2+1 viscous hydro code
  extern void readParameters(const char*);


  readParameters("data/params.txt");

  ofstream ofs,ofs2,ofs3,ofs4;

    ofs.open("data/initedold.dat");  //file to be read when running vh2
    ofs2.open("data/initgpold.dat");  //file that is more useful for visualizing with gnuplot
    ofs3.open("data/inited.dat"); 
    ofs4.open("data/initgp.dat");
  //locate starting temperature in EOS
  long int i;

  loadeos();


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


//  Most of he following code is taken from the "eccentricity" subroutine

    //....Estimate NN cross sections.
    double sig = hadronxsec::totalXsection(ecm,0);
    double sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
    double sigin = sig - sigel;
//     double sigin = 70;

    double massnumber = 197; //Gold

  if (int(IC)/10 == -2)
  {
    cout << "Lead-Lead collision\n";
    ecm = 5500; //LHC (5.5 A*TeV)
//     Rnuc = 6.6; //Lead
//     anuc = 0.55; //Lead
    sigin = 60; //LHC (5.5 A*TeV)
    massnumber = 208;  //Lead
  }


    cout << "# ecm = " << ecm;
    cout << "  sig= " << sig << " sigel= " << sigel << " sigin= " << sigin
	  << endl;
    cout << "# lambda= " << lambda << " mode= " << mode
	 << " K= " << kglue << endl;

    OverLap* overlap = new OverLap(massnumber,0.0,sigin);
//     OverLap* overlap = new OverLap(197,0.0,sigin);
//     OverLap* overlap = new OverLap(208,0.0,sigin);  //Lead

    //UnintegPartonDist* wf = new KLNfunc();
    UnintegPartonDist* wf = new KLNGamfunc(lambda);
    //UnintegPartonDist* wf = new FTDipole();
    //KLNModel2* kln = new KLNModel2(overlap,ecm,mode,wf);
    KLNModel* kln = new KLNModel(overlap,ecm,mode,wf);
    kln->setOptFixAlpha(1);
    kln->setNoramlization(1.0);
    kln->setLambdaSaturation(lambda);
    kln->setdEtdy(Et);
    kln->setdEtdy(0);
    kln->loadHist();
    kln->setKGlue(kglue);
	kln->setOptTA(2);


    const int maxx=(NUMT-1)/2;
    const int maxy=(NUMT-1)/2;
//    double dx = 0.3;

    double rapidity=0.0;

//		double b = 1.; //b0 + ib*db;

//	ofs << "# b= " << b << endl;

    double npartd = overlap->getNPartDensity(B);
//    double npart = overlap->getNPart(b);
//    double ncoll = overlap->getNBC(b);
    double qs2 = kln->SaturationScaleX(0.01,npartd/2.0);

    double hbarc = 0.197326963;
        
    double norm = 1.2/0.6 * kln->getdNdy(0.,0.,0.,rapidity);
cout << "e0 = " << e0 <<"\n";
cout << "norm = " << norm << "\n";
    for(int ix=-maxx;ix<=maxx;ix++) 
    for(int iy=-maxy;iy<=maxy;iy++) {
	double x = AT*hbarc*ix;
	double y = AT*hbarc*iy;
//	double fac=4.0;
//	if(ix==0 && iy==0) fac=1.0;
//	if(ix==0 && iy!=0) fac=2.0;
//	if(ix!=0 && iy==0) fac=2.0;

	double phi = atan2(y,x);
	double stretch = 1.0;
	if (IC == -6 || IC == -7 || IC == -27) 
	  stretch = sqrt(1 + BIEPS*cos(2*(phi-BIANGLE)) + TRIEPS*cos(3*(phi-TRIANGLE))
	    			+ QUADEPS*cos(4*(phi-QUADANGLE)) + QUINTEPS*cos(5*(phi-QUINTANGLE))
	    			+ SEXEPS*cos(6*(phi-SEXANGLE)) + SEPTEPS*cos(7*(phi-SEPTANGLE)));

	else stretch = 1.0;
        double dndy = 1.2/0.6 * kln->getdNdy(B,x*stretch,y*stretch,rapidity);
	double dndynorm = dndy * e0  / norm;
//	double en =  pow(hbarc,2)*1.71822665 * pow(dndy/2.,4./3.) /(pow(qs2,0.5) * TINIT);
	double en = pow(dndy,4./3.)*e0 / pow(norm,4./3.);
	
//	cout << dndy << "\t";

//	if(dndy == 0.) dndy = 1.0e-13;

	ofs //<< x << setw(13) << y
	    << dndynorm << "\t";
//	    << endl;
	ofs2  << AT*ix <<  "\t" << AT*iy << "\t" << dndynorm << "\n";
	ofs3 //<< x << setw(13) << y
	    << en << "\t";
//	    << endl;
	ofs4  << AT*ix <<  "\t" << AT*iy << "\t" << en << "\n";
	
	  if ((IC == -6) || (IC = -7) || IC == -27)
	      {
		    //  Calculate eccentricities epsilon1 to epsilon7 
		    //  with energy density weighting
		    double r2 = x*x + y*y;
		    double r3 = pow(r2,1.5);
		    double r4 = r2*r2;
		    double r5 = r2*r3;
		    double r6 = r2*r2*r2;
		    double r7 = r4*r3;
		    
		    //  denominators:
		    den += r2*en;
		    den3 += r3*en;
		    den4 += r4*en;
		    den5 += r5*en;
		    den6 += r6*en;
		    den7 += r7*en;
		    //  factors for the numerator
		    num2 += (y*y - x*x)*en;// standard eccentricity
		    numeps11 += r3*cos(phi)*en;
		    numeps12 += r3*sin(phi)*en;
		    numeps21 += r2*cos(2*phi)*en;// for participant eccentricity
		    numeps22 += r2*sin(2*phi)*en;// for psrticipant eccentricity
		    numeps31 += r3*cos(3*phi)*en;// for triangularity with r^n weighting
		    numeps32 += r3*sin(3*phi)*en;// for triangularity with r^n weighting
		    numeps41 += r4*cos(4*phi)*en;
		    numeps42 += r4*sin(4*phi)*en;
		    numeps51 += r5*cos(5*phi)*en;
		    numeps52 += r5*sin(5*phi)*en;
		    numeps61 += r6*cos(6*phi)*en;
		    numeps62 += r6*sin(6*phi)*en;
		    numeps71 += r7*cos(7*phi)*en;
		    numeps72 += r7*sin(7*phi)*en;
		    // r^2 weighting
		    r2numeps31 += r2*cos(3*phi)*en;// for triangularity with r^2 weighting
		    r2numeps32 += r2*sin(3*phi)*en;// for triangularity with r^2 weighting
		    r2numeps41 += r2*cos(4*phi)*en;
		    r2numeps42 += r2*sin(4*phi)*en;
		    r2numeps51 += r2*cos(5*phi)*en;
		    r2numeps52 += r2*sin(5*phi)*en;
		    r2numeps61 += r2*cos(6*phi)*en;
		    r2numeps62 += r2*sin(6*phi)*en;
		    r2numeps71 += r2*cos(7*phi)*en;
		    r2numeps72 += r2*sin(7*phi)*en;
		    
		    //  entropy density weighting calculations
		    double s = 0.0;
		    double p, T;
		    if (en < lowestE)
		    {
			    p = AT*AT*AT*AT*gsl_spline_eval(pspline,en,pacc);
			    T = AT*gsl_spline_eval(Tspline,en,pacc);
		    }
		    else
		    {
			    p = AT*AT*AT*AT*en*gsl_spline_eval(cs2spline,lowestE,cs2acc);
			    T = sqrtl(sqrtl(en/eoT4[0]));
		    }
		    s = (en + p)/T;
		    
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
	  mx += x*x*en;
	  my += y*y*en;
	  mn += en;
	
}

  if ((IC == -6) || (IC = -7) || (IC = -27))
  {  
    
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
  
  cout << "overlap " << 4*M_PI*sqrt(mx*my)/mn << endl;

}

// xxx End inserted code


void eccentricity(double ecm, double lambda, double kglue, int mode,int Et)
{

    ofstream ofs,ofs2;
    ofs.open("dens.dat");
    ofs2.open("ecc.dat");

    //....Estimate NN cross sections.
    double sig = hadronxsec::totalXsection(ecm,0);
    double sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
    double sigin = sig - sigel;
    cout << "# ecm = " << ecm;
    cout << "  sig= " << sig << " sigel= " << sigel << " sigin= " << sigin
	  << endl;
    cout << "# lambda= " << lambda << " mode= " << mode
	 << " K= " << kglue << endl;
    cout << "# b    eps_npart    eps_cgc    Qs2(b)    Qsmax     Npart" << endl;

    OverLap* overlap = new OverLap(197,0.0,sigin);
    /*
    cout << 0.0 << "  " <<  overlap->getNPart(0.0) << endl;
    cout << 0.2*6.65 << "  " <<  overlap->getNPart(0.2*6.65) << endl;
    cout << 0.4*6.65 << "  " <<  overlap->getNPart(0.4*6.65) << endl;
    cout << 0.6*6.65 << "  " <<  overlap->getNPart(0.6*6.65) << endl;
    cout << 0.8*6.65 << "  " <<  overlap->getNPart(0.8*6.65) << endl;
    cout << 1.0*6.65 << "  " <<  overlap->getNPart(1.0*6.65) << endl;
    cout << 1.2*6.65 << "  " <<  overlap->getNPart(1.2*6.65) << endl;
    cout << 1.4*6.65 << "  " <<  overlap->getNPart(1.4*6.65) << endl;
    cout << 1.6*6.65 << "  " <<  overlap->getNPart(1.6*6.65) << endl;
    cin.get();
    */


    //UnintegPartonDist* wf = new KLNfunc();
    UnintegPartonDist* wf = new KLNGamfunc(lambda);
    //UnintegPartonDist* wf = new FTDipole();
    //KLNModel2* kln = new KLNModel2(overlap,ecm,mode,wf);
    KLNModel* kln = new KLNModel(overlap,ecm,mode,wf);
    kln->setOptFixAlpha(1);
    kln->setNoramlization(1.0);
    kln->setLambdaSaturation(lambda);
    kln->setdEtdy(Et);
    kln->setdEtdy(0);
    kln->loadHist();
    kln->setKGlue(kglue);
	kln->setOptTA(2);

    double npart0 = overlap->getNPartDensity(0.0,0.0,0.0);
    double ncoll0 = overlap->getTAA(0.0,0.0,0.0);

    const int maxx=40;
    const int maxy=40;
    double dx = 0.3;

    double rapidity=0.0;
    double bmax=14;
    int maxb=14;
    //maxb=1;
    double b0=0.0;
    double db = bmax/maxb;

    //for (int ib=0; ib<=maxb;ib++) {
    for (int ib=maxb; ib>=0;ib--) {
		double b = b0 + ib*db;

	ofs << "# b= " << b << endl;

	/*
    double epsx1=0.0, epsy1=0.0;
    if(abs(rapidity)>0) {
    for(int ix=0;ix<maxx;ix++) 
    for(int iy=0;iy<maxy;iy++) {
	double x = dx*ix;
	double y = dx*iy;
	double dndy = 1.2/0.6 * kln->getdNdy(b,x,y,rapidity);
	epsx1 += x*dndy;
	epsy1 += y*dndy;
    }
    }
    */

    double epsx=0.0;
    double epsy=0.0;
    double epsxNC=0.0;
    double epsyNC=0.0;
    double epsx1=0.0;
    double epsy1=0.0;
    double epsx2=0.0;
    double epsy2=0.0;
    double epsx3=0.0;
    double epsy3=0.0;
    double epsxCGC=0.0;
    double epsyCGC=0.0;
    double qsmax=0.0;
    for(int ix=0;ix<maxx;ix++) 
    for(int iy=0;iy<maxy;iy++) {
	double x = dx*ix;
	double y = dx*iy;
	double fac=4.0;
	if(ix==0 && iy==0) fac=1.0;
	if(ix==0 && iy!=0) fac=2.0;
	if(ix!=0 && iy==0) fac=2.0;

	// Npart or Binary Collison scaling
	double npart1 =  overlap->getNPart1(x,y,b);
	double npart2 =  overlap->getNPart1(x,y,-b);
	double T1 = overlap->getThickness(x+b/2.0,y); 
	double T2 = overlap->getThickness(x-b/2.0,y); 
	double probA = overlap->getProb(x,y,b);
	double probB = overlap->getProb(x,y,-b);
	double P1=0.0, P2=0.0;
	if(probA > 1e-6 && probB > 1e-6) {
	    P1 = overlap->getThickness(x+b/2.0,y)/probA; 
	    P2 = overlap->getThickness(x-b/2.0,y)/probB; 
	}
	double npmin = min(npart1,npart2);
	//double npmax = max(npart1,npart2);
	double Tmin = min(T1,T2);
	//double Tmax = max(T1,T2);
	double Pmin = min(P1,P2);
	//double Pmax = max(P1,P2);
	//double np = npmin*npmin*(npmax - 2.0/3.0*npmin);

	epsx1 += x*x*fac*npmin;
	epsy1 += y*y*fac*npmin;
	epsx2 += x*x*fac*Tmin;
	epsy2 += y*y*fac*Tmin;
	epsx3 += x*x*fac*Pmin;
	epsy3 += y*y*fac*Pmin;

	//double np = npart1+npart2;
	double npart = 45/npart0 * overlap->getNPartDensity(x,y,b);
	double ncoll = 45/ncoll0 * overlap->getTAA(x,y,b);
	epsx += x*x*fac*npart;
	epsy += y*y*fac*npart;
	epsxNC += x*x*fac*ncoll;
	epsyNC += y*y*fac*ncoll;

	//epsx1 += x*x*pow(np,1.5);
	//epsy1 += y*y*pow(np,1.5);

	// CGC
	//double dndy = 1.0;
	double dndy = 1.2/0.6 * kln->getdNdy(b,x,y,rapidity);
	epsxCGC += x*x*fac*dndy;
	epsyCGC += y*y*fac*dndy;

	// very simple CGC
        //double qs2 = kln->SaturationScaleX(0.01,npart);
	//double dndy2 = npart*log((qs2+0.04)/0.04);
	//epsx2 += x*x*fac*dndy2;
	//epsy2 += y*y*fac*dndy2;

	if(ix==0 && iy==0) qsmax = kln->getQs2();

	//kln->setNPartDensity(npart);
	//Qs2 = kln->getSaturationScale();

	if(y==0)
	ofs << x << setw(13) << y
	    << setw(13) << npart
	    << setw(13) << ncoll
	    << setw(13) << dndy
	    << endl;

	//if(iy==0) cout << x << setw(13) << npart << setw(13) << dndy << endl;
    }

    double npartd = overlap->getNPartDensity(b);
    double npart = overlap->getNPart(b);
    double ncoll = overlap->getNBC(b);
    double qs2 = kln->SaturationScaleX(0.01,npartd/2.0);

       cout << b
	    << setprecision(5) << setw(14) << (epsy-epsx)/(epsx+epsy)
	   //<< setprecision(4) << setw(14) << (epsy1-epsx1)/(epsx1+epsy1)
	   //<< setprecision(4) << setw(14) << (epsy2-epsx2)/(epsx2+epsy2)
	   //<< setprecision(4) << setw(14) << (epsy3-epsx3)/(epsx3+epsy3)
	    << setprecision(5) << setw(14) << (epsyNC-epsxNC)/(epsxNC+epsyNC)
	    << setprecision(5) << setw(14) << (epsyCGC-epsxCGC)/(epsxCGC+epsyCGC)
	   << setprecision(4) << setw(10) << qs2
	   << setprecision(4) << setw(10) << qsmax
	   << endl;

       ofs2 << "    {" << b << ", " << npart << ", " << ncoll << ", "
	    << (epsy-epsx)/(epsx+epsy)  << ",  "
	    << (epsyNC-epsxNC)/(epsxNC+epsyNC)  << ",  "
	    << (epsyCGC-epsxCGC)/(epsxCGC+epsyCGC)  << "},"
	   << endl;
    }
       ofs2 << "    };" << endl;

       ofs2.close();
}

// use Gauss integral to get dN/dy from local Qs(x,y,b).
void dndyLocal2(double ecm, double b, double lambda=0.25,double kglue=1.0,
       	int mode=11, int Et=0, double hmass=0.0)
{
    //....Estimate NN cross sections.
    double sig = hadronxsec::totalXsection(ecm,0);
    double sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
    double sigin = sig - sigel;
    cout << "# ecm = " << ecm;
    cout << "  sig= " << sig << " sigel= " << sigel << " sigin= " << sigin
	  << endl;
    cout << "# lambda= " << lambda << " model= " << mode << endl;

    OverLap* overlap = new OverLap(197,b,sigin);

    Saturation* kln=0;

    if(mode <10) {
    //UnintegPartonDist* wf= new KLNfunc();
    UnintegPartonDist* wf = new KLNGamfunc(lambda);
    //UnintegPartonDist* wf = new FTDipole();
	kln = new KLNModel(overlap,ecm,mode,wf);
    } else {
	cerr << " invalid mode = " << mode << endl;
	exit(1);
    }

    kln->setOptFixAlpha(1);
    kln->setNoramlization(1.0);
    kln->setLambdaSaturation(lambda);
    if(hmass > 0.0) {
	kln->setdNdeta(1);   // dN/deta or dN/dy
	kln->setMHadron(hmass);
    }
    kln->setdEtdy(Et);
    kln->loadHist();
    kln->setKGlue(kglue);

    int format = 1;
    double ymin=-5.0;
    double ymax= 5.0;
    //const int    ny = 12;
    //const int    ny = 23;
    //double dy = (ymax-ymin)/(ny+1);
    //double dy = 0.5;
    double dy = 1.0;
    int ny = (int)((ymax-ymin)/dy)+1;
    cout << "# np= " << ny << endl;
    //double dndy[ny*2];
    for (int i=0; i<ny; i++) {
	double y = i*dy + ymin;
	double dndy2 = kln->getdNdyInteg(b,y);
	//double dndy2 = 0.0;
	double dndy = 0.0;
	if(mode <= 10) dndy = kln->getdNdy(b,y);
	//dndy[i]=dndy2;
	if(format == 1) {
	    cout << setw(6) <<  y << setw(10) << dndy2
		<< setw(10) << dndy << endl;
	} else {
	    cout << "{" <<  y << ",  "  <<  dndy2
		<< ",  " << dndy << "}," << endl;
	}
	   


    }

    kln->printHist();

    delete kln;
    delete overlap;

}

// centraity dependence of multiplicity at y=0.
void dndyCentralityAtMidRap(double ecm, double lambda,double kglue,int mode,
	int Et)
{
    ofstream ofs;
    ofs.open("cent.dat");
    double b = 0.0;

    cout << "#"<<endl;
    cout << "# This results are produced by dndyCentralityAtMidRap"
	<< endl;
    cout << "#"<<endl;

    //....Estimate NN cross sections.
    double sig = hadronxsec::totalXsection(ecm,0);
    double sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
    double sigin = sig - sigel;
    cout << "# ecm = " << ecm;
    cout << "  sig= " << sig << " sigel= " << sigel << " sigin= " << sigin
	  << endl;
    cout << "# lambda= " << lambda << endl;

    OverLap* overlap2 = new OverLap(197,b,sigin);
    //OverLap* overlap = new OverLap(197,b,39.5277);
    //OverLap* overlap = new OverLap(197,b,41.975);
    OverLap* overlap = overlap2;

    KLNModel *kln=0;
    ///UnintegPartonDist* wf= new KLNfunc();
    UnintegPartonDist* wf = new KLNGamfunc(lambda);
    //UnintegPartonDist* wf = new FTDipole();
	kln = new KLNModel(overlap,ecm,mode,wf);
     kln->setOptFixAlpha(0);
	kln->setNoramlization(1.0);
	kln->setLambdaSaturation(lambda);
	kln->setdNdeta(0);
	kln->setdEtdy(Et);
	kln->setKGlue(kglue);
	kln->setOptTA(2);

    double delta=0.79;
    double bmax = 12.0;
    int nbmax = 12;
    double db = bmax/nbmax;
    double h = 0.0;

    // loop over impact parameter.
    for(int ib=0;ib<nbmax;ib++) {
	b = ib*db;
	double npart = overlap2->getNPart(b);
	double npartdens = overlap->getNPartDensity(b);
	kln->setNPartDensity(npartdens);
	double Qs2 = kln->getSaturationScale();

	//kln->setModel(2);  // KL model.
	//double dndy2local = kln->getdNdyInteg(b,h);
	//double dndy2 = kln->getdNdy(b,h);

	kln->setModel(mode);  // KLN model.
	double dndy3local = kln->getdNdyInteg(b,h);
	double dndy3 = kln->getdNdy(b,h);

	//double QsAve=kln->getAverageQs(b);
	double fac=pow(ecm/200,lambda);
	//Qs2 = kln->SaturationScale0(npartdens/2.0);
	Qs2 = kln->SaturationScaleX(0.01,npartdens/2.0);

	cout << setprecision(5) <<setw(6) << npart
	     << setprecision(4) <<setw(10)<<  2*dndy3/npart
	     << setprecision(4) <<setw(10)<<  2*dndy3local/npart
	     << setprecision(3) <<setw(4) << b
	     << setprecision(4) <<setw(9)<< 2*fac*log(Qs2*fac/0.2/0.2)
	     << setprecision(4) <<setw(9)<< 2*fac*pow(npart,(1-delta)/(3*delta))
	     //<< setprecision(4) <<setw(9)<<  kln->getQs2()
	     << setprecision(4) <<setw(9)<<  Qs2
	     << setprecision(4) <<setw(9)<<  npartdens
	     << endl;

	ofs << setprecision(5) <<setw(6) << "{" << npart << ","
	     << setprecision(4) <<setw(9)<<  2*dndy3local/npart << ","
	     << setprecision(3) <<setw(4) << b << ","
	     << setprecision(4) <<setw(9)<< 2*fac*log(Qs2*fac/0.2/0.2)
	     << "},"
	     << endl;
    }

    delete kln;
    delete overlap;

    ofs << "    };" << endl;
    ofs.close();

}

void dndyLocal()
{
    double b=2.0;
    //double b=3.6;
    //double b=4.0;
    //double b=4.5;
    //double b=5.5;
    //double b=6.3;
    //double b=7.2;
    //double b=8.5;
    //double b=10.0;

    double sig=42.0;
    double ecm=200;
    int mode=3;

    OverLap* overlap = new OverLap(197,b,sig);
    UnintegPartonDist* wf= new KLNfunc();
    KLNModel *kln = new KLNModel(overlap,ecm,mode,wf);

    kln->setNoramlization(1.0);
    kln->setLambdaSaturation(0.25);
    kln->setdNdeta(0);

    const int maxx=40;
    const int maxy=40;
    double dx = 0.3;
    double** Qs2 = new double* [maxx];
    for(int ix=0;ix<maxx;ix++)
	Qs2[ix] = new double[maxy];


    double ymin=-3.0;
    double ymax= 3.0;
    int    ny = 11;
    double dy = (ymax-ymin)/ny;

    for(int ix=0;ix<maxx;ix++) 
    for(int iy=0;iy<maxy;iy++) {
	double x = dx*ix;
	double y = dx*iy;
	double npart = overlap->getNPartDensity(x,y,b);
	kln->setNPartDensity(npart);
	Qs2[ix][iy] = kln->getSaturationScale();
	if(Qs2[ix][iy] < 0.0) {
	    cout << "qs= " << Qs2[ix][iy] << endl;
	   Qs2[ix][iy]=0.0;
	}
    }


    for (int i=0; i<ny; i++) {
	double y = i*dy + ymin;
	double dNdy=0.0, dNdeta=0.0;
	for(int ix=0;ix<maxx;ix++) 
	for(int iy=0;iy<maxy;iy++) {

	    kln->setSaturationScale(Qs2[ix][iy]);

	    int fac = 1;
	    if(ix == 0 && iy != 0 ) fac=2;
	    else if(ix != 0 && iy == 0 ) fac=2;
	    else if(ix != 0 && iy != 0 ) fac=4;

	    double dndy = kln->getdNdy(y);
	    dNdy += dndy*dx*dx*fac;

	    kln->setdNdeta(1);
	    double dndeta = kln->getdNdy(y);
	    kln->setdNdeta(0);
	    dNdeta += dndeta*dx*dx*fac;

	    //cout << " Qs= " << Qs2[ix][iy] << endl;
	    //cout << "ix= " << ix << " iy= " << iy << " dndy= " << dndy << endl;
	    //cin.get();

	}

	cout << y << setw(15) << dNdy << setw(15) << dNdeta << endl;
    }

    //double qs2= kln->getQs2();
    //cout << "# Qs2= " << qs2 << endl;


    for(int ix=0;ix<maxx;ix++) delete [] Qs2[ix];
    delete [] Qs2;
}

void dndy()
{
    double b=2.0;
    //double b=4.5;
    //double b=6.3;
    double sig=42.0;
    double ecm=200;
    int mode=2;

    OverLap* overlap = new OverLap(197,b,sig);
    UnintegPartonDist* wf= new KLNfunc();
    KLNModel *kln = new KLNModel(overlap,ecm,mode,wf);
    kln->setNoramlization(1.0);
    kln->setLambdaSaturation(0.3);
    kln->setdNdeta(0);

    double npart = overlap->getNPart(b);
    double npart_dens = overlap->getNPartDensity(b);
    kln->setNParticipants(npart);
    kln->setNPartDensity(npart_dens);
    double qs2 = kln->getSaturationScale();
    kln->setSaturationScale(qs2);


    double ymin=-6.5;
    double ymax=5.5;
    int    ny = 30;
    double dy = (ymax-ymin)/ny;
    for (int i=0; i<ny; i++) {
	double y = i*dy + ymin;
	double dndy = kln->getdNdy(y);
	kln->setdNdeta(1);
	double dndeta = kln->getdNdy(y);
	kln->setdNdeta(0);
	cout << y << setw(15) << dndy << setw(15) << dndeta << endl;
    }
    //double qs2= kln->getQs2();
    cout << "# Qs2= " << qs2 << endl;
}

void bDepQs(double ecm)
{
    //....Estimate NN cross sections.
    double sig = hadronxsec::totalXsection(ecm,0);
    double sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
    double sigin = sig - sigel;
    int mode=2;
    double xq = sqrt(2.0)/ecm;

    OverLap* overlap = new OverLap(197,0.0,sigin);
    UnintegPartonDist* wf= new KLNfunc();
    KLNModel *kln = new KLNModel(overlap,ecm,mode,wf);
    kln->setNoramlization(1.0);
    kln->setLambdaSaturation(0.3);
    kln->setdNdeta(0);

    // Calculate total AA cross section to get the centrality.
    int nb=200;
    double db=0.1;
    double* sigAA = new double [nb];
    double* sigAAsum = new double [nb];
    sigAAsum[0]=0.0;
    for(int i=0;i<nb;i++) {
	double b=i*db;
	double nbc= overlap->getNBC(b);
	sigAA[i]=2*M_PI*b*0.1*(1-exp(-nbc*sigin/10));
	sigAAsum[i] = sigAA[i];
	if(i > 0) sigAAsum[i] += sigAAsum[i-1];
    }
    cout << "# cross section " << sigAAsum[199]/100 << " barn" <<endl;

    //double bmax=12.0;
    cout << "# x= " << xq << endl;
    cout << "# b (fm)   Ncoll    Npart      rho       Qs^2(GeV^2)";
    cout <<  " sigma_AA   centrality(%)" << endl;

    for(int i=0;i<nb;i++) {
	double b = i*db;
	double npart = overlap->getNPart(b);
	double npartdens = overlap->getNPartDensity(b);
	kln->setNParticipants(npart);
	kln->setNPartDensity(npartdens);
	double qs2, xg;
	    qs2 = kln->getSaturationScale();
	    xg = kln->xG(qs2);

	cout << setprecision(3) << setw(4) << b
             << setprecision(7) << setw(10) << overlap->getNBC(b)*sigin/10
	     << setprecision(7) << setw(10) << npart
	     << setprecision(7) << setw(10) << npartdens
	     //<< setprecision(7) << setw(10) << qs2
	     //<< setw(10) << xg
	     //<< setprecision(7) << setw(12) << sigAA[i]
	     << setprecision(7) << setw(12) << sigAAsum[i]/sigAAsum[199]*100
	     << endl;
    }


}

void CentralityBin(double ecm)
{
    //....Estimate NN cross sections.
    double sig = hadronxsec::totalXsection(ecm,0);
    double sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
    double sigin = sig - sigel;

    OverLap* overlap = new OverLap(197,0.0,sigin);
    CentralityCut* cen = new CentralityCut(overlap);

    //overlap->setRadius(6.38);
    //overlap->setDiffuseness(0.53);

    cout << "# rad= " << overlap->getRadius() << endl;
    cout << "# difuseness= " << overlap->getDiffuseness() << endl;
    double sigAA = cen->getTotalCross();

    cout << "# cross section " << sigAA/100 << " barn" <<endl;

    //double bm = cen->print(0,100,0.0);
    // PHENIX
    double b2 = cen->print(0,10,0.0);
    double b3 = cen->print(10,20,b2);
    double b4 = cen->print(20,30,b3);
    double b5 = cen->print(30,40,b4);
    double b6 = cen->print(40,50,b5);
    //double b7 = cen->print(50,60,b6);
    //double b8 = cen->print(60,80,b7);
    //double b9 = cen->print(60,92,b7);

    // PHOBOS
    /*
    double b1 = cen->print(0,6,0.0);
    double b2 = cen->print(6,15,b1);
    double b3 = cen->print(15,25,b2);
    double b4 = cen->print(25,35,b3);
    double b5 = cen->print(35,45,b4);
    */

}

void eccentricity2(double ecm, double lambda, double kglue, int mode,int Et)
{

    ofstream ofs,ofs2;
    ofs.open("dens.dat");
    ofs2.open("ecc.dat");

    //....Estimate NN cross sections.
    double sig = hadronxsec::totalXsection(ecm,0);
    double sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
    double sigin = sig - sigel;
    cout << "# ecm = " << ecm;
    cout << "  sig= " << sig << " sigel= " << sigel << " sigin= " << sigin
	  << endl;
    cout << "# lambda= " << lambda << " mode= " << mode
	 << " K= " << kglue << endl;
    cout << "b    eps_npart    eps_cgc    Qs2(b)    Qsmax     Npart" << endl;

    OverLap* overlap = new OverLap(197,0.0,sigin);

    //UnintegPartonDist* wf = new KLNfunc();
    UnintegPartonDist* wf = new KLNGamfunc(lambda);
    //UnintegPartonDist* wf = new FTDipole();
    //KLNModel2* kln = new KLNModel2(overlap,ecm,mode,wf);
    KLNModel* kln = new KLNModel(overlap,ecm,mode,wf);
    kln->setOptFixAlpha(0);
    kln->setNoramlization(1.0);
    kln->setLambdaSaturation(lambda);
    kln->setdEtdy(Et);
    kln->setdEtdy(1);
    kln->loadHist();
    kln->setKGlue(kglue);

    double npart0 = overlap->getNPartDensity(0.0,0.0,0.0);
    double ncoll0 = overlap->getTAA(0.0,0.0,0.0);

    const int maxx=40;
    const int maxy=40;
    double dx = 0.3;

    double rapidity=1.5;
    double bmax=14;
    int maxb=14;
    //maxb=1;
    double b0=0.0;
    double db = bmax/maxb;

    //for (int ib=0; ib<=maxb;ib++) {
    for (int ib=maxb; ib>=0;ib--) {
		double b = b0 + ib*db;

	ofs << "# b= " << b << endl;

    double epsx1=0.0, epsy1=0.0;
    if(abs(rapidity)>0) {
    for(int iy=0;iy<maxy;iy++)
    for(int ix=0;ix<maxx;ix++) { 
	double x = dx*ix;
	double y = dx*iy;
	double dndy  = kln->getdNdy(b,x,y,rapidity);
	double dndy2 = kln->getdNdy(b,-x,y,rapidity);

	    epsx1 += x*dndy;
	    epsy1 += y*dndy;
	if(ix==0 && iy !=0) {
	    epsx1 += x*dndy;
	    epsy1 += -y*dndy;
	} else if(ix !=0 && iy == 0 ) {
	    epsx1 += -x*dndy2;
	    epsy1 += y*dndy2;
	} else if(ix !=0 && iy != 0 ) {
	    epsx1 += x*dndy;
	    epsy1 += -y*dndy;
	    epsx1 += -x*dndy2;
	    epsy1 += y*dndy2;
	    epsx1 += -x*dndy2;
	    epsy1 += -y*dndy2;
	}

	/*
	cout << x << "  " << y << " dndy= "  <<dndy << endl;
	cout << " (x,-y) " << kln->getdNdy(b,x,-y,rapidity) << endl;
	cout << " (-x,y) " << kln->getdNdy(b,-x, y,rapidity) << endl;
	cout << " (-x,-y) " << kln->getdNdy(b,-x, -y,rapidity) << endl;
	cin.get();
	*/

    }
    }

    double epsx=0.0;
    double epsy=0.0;
    double epsxNC=0.0;
    double epsyNC=0.0;
    double epsxCGC=0.0;
    double epsyCGC=0.0;
    double epsxCGC1=0.0;
    double epsyCGC1=0.0;
    double qsmax=0.0;
    for(int ix=0;ix<maxx;ix++) 
    for(int iy=0;iy<maxy;iy++) {
	double x = dx*ix;
	double y = dx*iy;

	// Npart or Binary Collison scaling
	//double npart1 =  overlap->getNPart1(x,y,b);
	//double npart2 =  overlap->getNPart1(x,y,-b);
	//double npmin = min(npart1,npart2);
	//double npmax = max(npart1,npart2);
	//double np = npmin*npmin*(npmax - 2.0/3.0*npmin);
	//double np = npart1+npart2;
	double npart = 45/npart0 * overlap->getNPartDensity(x,y,b);
	double ncoll = 45/ncoll0 * overlap->getTAA(x,y,b);
	epsx += x*x*npart;
	epsy += y*y*npart;
	epsxNC += x*x*ncoll;
	epsyNC += y*y*ncoll;

	// CGC
	double dndy = 1.2/0.6 * kln->getdNdy(b,x,y,rapidity);
	epsxCGC += x*x*dndy;
	epsyCGC += y*y*dndy;
	if(dndy>3.0) {
	    epsxCGC1 += x*x*dndy;
	    epsyCGC1 += y*y*dndy;
	}

	// very simple CGC
        //double qs2 = kln->SaturationScaleX(0.01,npart);
	//double dndy2 = npart*log((qs2+0.04)/0.04);
	//epsx2 += x*x*dndy2;
	//epsy2 += y*y*dndy2;

	if(ix==0 && iy==0) qsmax = kln->getQs2();

	//kln->setNPartDensity(npart);
	//Qs2 = kln->getSaturationScale();

	if(y==0)
	ofs << x << setw(13) << y
	    << setw(13) << npart
	    << setw(13) << ncoll
	    << setw(13) << dndy
	    << endl;

	//if(iy==0) cout << x << setw(13) << npart << setw(13) << dndy << endl;
    }

    double npartd = overlap->getNPartDensity(b);
    double npart = overlap->getNPart(b);
    double ncoll = overlap->getNBC(b);
    double qs2 = kln->SaturationScaleX(0.01,npartd/2.0);

       cout << b
	    << setprecision(5) << setw(14) << (epsy-epsx)/(epsx+epsy)
	    << setprecision(5) << setw(14) << (epsyNC-epsxNC)/(epsxNC+epsyNC)
	    << setprecision(5) << setw(14) << (epsyCGC-epsxCGC)/(epsxCGC+epsyCGC)
	   << setprecision(4) << setw(10) << qs2
	   << setprecision(4) << setw(10) << qsmax
	   << endl;

       ofs2 << "    {" << b << ", " << npart << ", " << ncoll << ", "
	    << (epsy-epsx)/(epsx+epsy)  << ",  "
	    << (epsyNC-epsxNC)/(epsxNC+epsyNC)  << ",  "
	    << (epsyCGC-epsxCGC)/(epsxCGC+epsyCGC)  << "},"
	   << endl;
    }
       ofs2 << "    };" << endl;

       ofs2.close();
}


void eccentricity3(double ecm, double lambda, double kglue, int mode,int Et)
{

    ofstream ofs,ofs2;
    ofs.open("dens.dat");
    ofs2.open("ecc.dat");

    //....Estimate NN cross sections.
    double sig = hadronxsec::totalXsection(ecm,0);
    double sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
    double sigin = sig - sigel;
    cout << "# ecm = " << ecm;
    cout << "  sig= " << sig << " sigel= " << sigel << " sigin= " << sigin
	  << endl;
    cout << "# lambda= " << lambda << " mode= " << mode
	 << " K= " << kglue << endl;

    cout << "b    eps_npart    eps_cgc    Qs2(b)    Qsmax     Npart" << endl;

    OverLap* overlap = new OverLap(197,0.0,sigin);

    //UnintegPartonDist* wf = new KLNfunc();
    UnintegPartonDist* wf = new KLNGamfunc(lambda);
    //UnintegPartonDist* wf = new FTDipole();
    KLNModel* kln = new KLNModel(overlap,ecm,mode,wf);
    kln->setOptFixAlpha(0);
    kln->setNoramlization(1.0);
    kln->setLambdaSaturation(lambda);
    kln->setdEtdy(Et);
    kln->setdEtdy(1);
    kln->loadHist();
    kln->setKGlue(kglue);

    double nnn = overlap->getNPart1(0.0,0.0,0.0);
    double qs2 = kln->SaturationScaleX(0.01,nnn);
    cout << "Qs= " << qs2 << " q= " << kln->SaturationScaleX(0.01,1.53)
	<< endl;

    double npart0 = overlap->getNPartDensity(0.0,0.0,0.0);
    double ncoll0 = overlap->getTAA(0.0,0.0,0.0);

    const int maxx=40;
    const int maxy=40;
    double dx = 0.3;

    double rapidity=2.0;
    double bmax=14;
    int maxb=14;
    bmax=10;
    maxb=10;
    double b0=0.0;
    double db = bmax/maxb;

    //for (int ib=0; ib<=maxb;ib++) {
    for (int ib=maxb; ib>=0;ib--) {
		double b = b0 + ib*db;

	ofs << "# b= " << b << endl;

    double epsx1=0.0;

	/*
	if(abs(rapidity)>0) {
    for(int ix=-maxx-1;ix<maxx;ix++)
    for(int iy=0;iy<maxy;iy++) {
	double x = dx*ix;
	double y = dx*iy;
	double dndy  = kln->getdNdy(b,x,y,rapidity);
	double fac=2.0;
	if(iy==0) fac=1.0;
	epsx1 += x*dndy*fac;
	vol += dndy*fac;

	if(dndy>0) {
	cout << x << "  " << y << " dndy= "  <<dndy << endl;
	cout << " (x,-y) " << kln->getdNdy(b,x,-y,rapidity) << endl;
	cout << " (-x,y) " << kln->getdNdy(b,-x, y,rapidity) << endl;
	cout << " (-x,-y) " << kln->getdNdy(b,-x, -y,rapidity) << endl;
	cin.get();
	}

    }
    epsx1 /= vol;
	}
    cout << " <x> = " << epsx1 << " <y>= " << epsy1 << " vol=" << vol*dx*dx;

    double z[38],zw[38];
    double zini=-10.0;
    double zfin=10.0;
    OverLap::Gauss38(zini,zfin,z,zw);
    double dNdy = 0.0;
    double dNdyx = 0.0;
    double dNdyy = 0.0;
    for(int ix=0;ix<38;ix++) {
	for(int iy=0;iy<38;iy++) {
	  double dndy = kln->getdNdy(b,z[ix],z[iy],rapidity);
	    dNdy += dndy*zw[ix]*zw[iy];
	    dNdyx += dndy*z[ix]*zw[ix]*zw[iy];
	    dNdyy += dndy*z[iy]*zw[ix]*zw[iy];
	}
    }
    dNdyx /= dNdy;
    dNdyy /= dNdy;
    cout << " <x> = " << dNdyx << " <y>= " << dNdyy << endl;

    //epsx1=dNdyx;
    //epsy1=dNdyy;
    */

    double epsx=0.0;
    double epsy=0.0;
    double epsxNC=0.0;
    double epsyNC=0.0;
    double epsxCGC=0.0;
    double epsyCGC=0.0;
    double qsmax=0.0;
    double volCGC=0.0;
    double x1CGC=0.0;
    double x3CGC=0.0;
    for(int ix=-maxx-1;ix<maxx;ix++) 
    for(int iy=-maxy-1;iy<maxy;iy++) {
	double x = dx*ix;
	double y = dx*iy;

	// Npart or Binary Collison scaling
	//double npart1 =  overlap->getNPart1(x,y,b);
	//double npart2 =  overlap->getNPart1(x,y,-b);
	//double npmin = min(npart1,npart2);
	//double npmax = max(npart1,npart2);
	//double np = npmin*npmin*(npmax - 2.0/3.0*npmin);
	//double np = npart1+npart2;
	double npart = 45/npart0 * overlap->getNPartDensity(x,y,b);
	double ncoll = 45/ncoll0 * overlap->getTAA(x,y,b);
	epsx += x*x*npart;
	epsy += y*y*npart;
	epsxNC += x*x*ncoll;
	epsyNC += y*y*ncoll;

	// CGC
	double dndy = kln->getdNdy(b,x,y,rapidity);
	epsxCGC += (x-epsx1)*(x-epsx1)*dndy;
	epsyCGC += y*y*dndy;
	volCGC += dndy;
	x1CGC += x*dndy;
	x3CGC += x*x*x*dndy;

	if(ix==0 && iy==0) qsmax = kln->getQs2();

	if(y==0)
	ofs << x << setw(13) << y
	    << setw(13) << npart
	    << setw(13) << ncoll
	    << setw(13) << dndy
	    << endl;

    }

    //double npartd = overlap->getNPartDensity(b);
    double npart = overlap->getNPart(b);
    double ncoll = overlap->getNBC(b);
    //double qs2 = kln->SaturationScaleX(0.01,npartd/2.0);

       cout << b
	    << setprecision(4) << setw(10) << (epsy-epsx)/(epsx+epsy)
	    << setprecision(4) << setw(10) << (epsyNC-epsxNC)/(epsxNC+epsyNC)
	    << setprecision(4) << setw(10) << (epsyCGC-epsxCGC)/(epsxCGC+epsyCGC)
	   << setprecision(4) << setw(10) << epsxCGC/volCGC
	   << setprecision(4) << setw(10) << epsyCGC/volCGC
	   << setprecision(4) << setw(10) << x1CGC/volCGC
	   << setprecision(4) << setw(10) << x3CGC/volCGC
	   << endl;

       ofs2 << "    {" << b << ", " << npart << ", " << ncoll << ", "
	    << (epsy-epsx)/(epsx+epsy)  << ",  "
	    << (epsyNC-epsxNC)/(epsxNC+epsyNC)  << ",  "
	    << (epsyCGC-epsxCGC)/(epsxCGC+epsyCGC)  << "},"
	   << endl;
    }
       ofs2 << "    };" << endl;

       ofs2.close();
}


void eccentricity4(double ecm, double lambda, double kglue, int mode,int Et)
{
    ofstream ofs,ofs2;
    ofs.open("dens.dat");
    ofs2.open("ecc.dat");

    //....Estimate NN cross sections.
    double sig = hadronxsec::totalXsection(ecm,0);
    double sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
    double sigin = sig - sigel;
    cout << "# ecm = " << ecm;
    cout << "  sig= " << sig << " sigel= " << sigel << " sigin= " << sigin
	  << endl;
    cout << "# lambda= " << lambda << " mode= " << mode
	 << " K= " << kglue << endl;

    cout << "b    eps_npart    eps_cgc    Qs2(b)    Qsmax     Npart" << endl;

    OverLap* overlap = new OverLap(197,0.0,sigin);

    UnintegPartonDist* wf = new KLNfunc();
    //UnintegPartonDist* wf = new KLNGamfunc(lambda);
    //UnintegPartonDist* wf = new FTDipole();
    KLNModel* kln = new KLNModel(overlap,ecm,mode,wf);
    kln->setOptFixAlpha(0);
    kln->setNoramlization(1.0);
    kln->setLambdaSaturation(lambda);
    kln->setdEtdy(Et);
    kln->setdEtdy(0);
    kln->loadHist();
    kln->setKGlue(kglue);
	kln->setOptTA(2);

    double nnn = overlap->getNPart1(0.0,0.0,0.0);
    double qs2 = kln->SaturationScaleX(0.01,nnn);
    cout << "# Qs= " << qs2 << " q= " << kln->SaturationScaleX(0.01,1.53) << endl;

    double npart0 = overlap->getNPartDensity(0.0,0.0,0.0);
    double ncoll0 = overlap->getTAA(0.0,0.0,0.0);

    const int maxx=40;
    const int maxy=40;
    double dx = 0.3;

    double bmax=2.0;
    int maxb=4;
    double b0=0.0;
    double db = bmax/maxb;
    double b=9.0;

    //for (int ib=0; ib<=maxb;ib++) {
    for (int ib=maxb; ib>=0;ib--) {
		double rapidity = b0 + ib*db;

	ofs << "# y= " << rapidity << endl;

    double epsx1=0.0, epsy1=0.0;
	double vol=0.0;

	if(abs(rapidity)>0) {
    for(int ix=-maxx-1;ix<maxx;ix++)
    for(int iy=0;iy<maxy;iy++) {
	double x = dx*ix;
	double y = dx*iy;
	double dndy  = kln->getdNdy(b,x,y,rapidity);
	double fac=2.0;
	if(iy==0) fac=1.0;
	epsx1 += x*dndy*fac;
	vol += dndy*fac;


    }
    epsx1 /= vol;
	}
    cout << "# <x> = " << epsx1 << " <y>= " << epsy1;

    double z[38],zw[38];
    double zini=-10.0;
    double zfin=10.0;
    OverLap::Gauss38(zini,zfin,z,zw);
    double dNdy = 0.0;
    double dNdyx = 0.0;
    double dNdyy = 0.0;
    for(int ix=0;ix<38;ix++) {
	for(int iy=0;iy<38;iy++) {
	  double dndy = kln->getdNdy(b,z[ix],z[iy],rapidity);
	    dNdy += dndy*zw[ix]*zw[iy];
	    dNdyx += dndy*z[ix]*zw[ix]*zw[iy];
	    dNdyy += dndy*z[iy]*zw[ix]*zw[iy];
	}
    }
    dNdyx /= dNdy;
    dNdyy /= dNdy;
    cout << "  <x> = " << dNdyx << " <y>= " << dNdyy << endl;

    //epsx1=dNdyx;
    //epsy1=dNdyy;

    double epsx=0.0;
    double epsy=0.0;
    double epsxNC=0.0;
    double epsyNC=0.0;
    double epsxCGC=0.0;
    double epsyCGC=0.0;
    double qsmax=0.0;
    double volCGC=0.0;
    double x1CGC=0.0;
    double x3CGC=0.0;
    for(int ix=-maxx-1;ix<maxx;ix++) 
    for(int iy=-maxy-1;iy<maxy;iy++) {
	double x = dx*ix;
	double y = dx*iy;

	// Npart or Binary Collison scaling
	double npart = 45/npart0 * overlap->getNPartDensity(x,y,b);
	double ncoll = 45/ncoll0 * overlap->getTAA(x,y,b);
	epsx += x*x*npart;
	epsy += y*y*npart;
	epsxNC += x*x*ncoll;
	epsyNC += y*y*ncoll;

	// CGC
	double dndy = kln->getdNdy(b,x,y,rapidity);
	double xx=x-epsx1;
	epsxCGC += xx*xx*dndy;
	epsyCGC += y*y*dndy;
	volCGC += dndy;
	x1CGC += x*dndy;
	x3CGC += xx*xx*xx*dndy;

	if(ix==0 && iy==0) qsmax = kln->getQs2();

	if(y==0)
	ofs << x << setw(13) << y
	    << setw(13) << npart
	    << setw(13) << ncoll
	    << setw(13) << dndy
	    << endl;

    }

    double npartd = overlap->getNPartDensity(b);
    double npart = overlap->getNPart(b);
    double ncoll = overlap->getNBC(b);
    double qs2 = kln->SaturationScaleX(0.01,npartd/2.0);

       cout << rapidity
	    //<< setprecision(4) << setw(10) << (epsy-epsx)/(epsx+epsy)
	    //<< setprecision(4) << setw(10) << (epsyNC-epsxNC)/(epsxNC+epsyNC)
	    << setprecision(4) << setw(12) << (epsyCGC-epsxCGC)/(epsxCGC+epsyCGC)
	   << setprecision(4) << setw(12) << epsyCGC/volCGC
	   << setprecision(4) << setw(12) << epsxCGC/volCGC
	   << setprecision(4) << setw(12) << x1CGC/volCGC
	   << setprecision(4) << setw(12) << x3CGC/volCGC
	   << endl;

       ofs2 << "    {" << b << ", " << npart << ", " << ncoll << ", "
	    << (epsy-epsx)/(epsx+epsy)  << ",  "
	    << (epsyNC-epsxNC)/(epsxNC+epsyNC)  << ",  "
	    << (epsyCGC-epsxCGC)/(epsxCGC+epsyCGC)  << "},"
	   << endl;
    }
       ofs2 << "    };" << endl;

       ofs2.close();
}




void MCAverageQs(double ecm)
{
    const int maxx=40;
    const int maxy=40;
    double dx = 0.3;
    MCSaturation* mcs = new MCSaturation(ecm,maxx,maxy,dx,dx);

    double qs1[maxx],qs2[maxx];
    double ta1[maxx],ta2[maxx];
    for(int ix=0;ix<maxx;ix++) {
	qs1[ix]=qs2[ix]=0.0;
	ta1[ix]=ta2[ix]=0.0;
    }

    double b=6.0;
    int nevent=20000;
    for(int iev=0;iev<nevent;iev++) {
	mcs->generate(b);
	for(int ix=0;ix<maxx;ix++) {
	    double y= 0.0;
	    double x= ix*dx;
	    qs1[ix] += mcs->getQs1(x,y);
	    qs2[ix] += mcs->getQs2(x,y);
	    ta1[ix] += mcs->getTA1(x,y);
	    ta2[ix] += mcs->getTA2(x,y);
	}
    }

    for(int ix=0;ix<maxx;ix++) {
	double y=0.0;
	double x=ix*dx;
	cout << x
	    << setw(12) << qs1[ix]/nevent
	    << setw(12) << qs2[ix]/nevent
	    << setw(12) << ta1[ix]/nevent
	    << setw(12) << ta2[ix]/nevent
	    << setw(12) << mcs->getNpart1(x,y,b)
	    << setw(12) << mcs->getNpart2(x,y,b)
	    << setw(12) << mcs->getTA1(x,y,b)
	    << setw(12) << mcs->getTA2(x,y,b)
	    << endl;
    }
}

