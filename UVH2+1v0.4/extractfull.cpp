#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <convert.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

using namespace std;

fstream pispec,ptarr;
fstream Kspec,pspec;
fstream respiv0,respiv2,respiv4;
fstream resKv0,resKv2,resKv4;
fstream respv0,respv2,respv4;
fstream resHv0,resHv2,resHv4;
fstream vallH, vallpis, vallKs, vallps;
// fstream vsinH, vsinpis, vsinKs, vsinps;

// int PTASIZE;
// int PHIPASIZE;

double * ptarray;
double * phiparray;

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
double PTMAX=4.2, TRIEPS, TRIANGLE = 0.0;
double QUADEPS = 0.0, QUADANGLE = 0.0, QUINTEPS, QUINTANGLE;
double SEXEPS, SEXANGLE, SEPTEPS, SEPTANGLE;
double BIEPS = 0.0, BIANGLE = 0.0;

//controls value of tau_Pi
double COEFF=3.0;

void doall()
{
  //pions
  

  ptarr >> PTASIZE;
  ptarr >> PHIPASIZE;

  printf("Found %i pt values and %i phip values\n",PTASIZE,PHIPASIZE);

  ptarray=new double[PTASIZE];
  phiparray=new double[4*PHIPASIZE];
  
  for (int i=0;i<PTASIZE;i++)
    {
      ptarr >> ptarray[i];
      cout << "ptarray[" << i << "] = " << ptarray[i] << endl;
    }

  for (int k=0;k<4*PHIPASIZE;k++) phiparray[k] = k*M_PI/2/PHIPASIZE;

  double resbuff[PTASIZE][4*PHIPASIZE];
  double CHbuff[PTASIZE][4*PHIPASIZE];

  double workhorsearr[4*PHIPASIZE+1];
  double workhorse[4*PHIPASIZE+1];
//   double workhorseangle[4*PHIPASIZE+1];
//   double workhorsequadangle[4*PHIPASIZE+1];

  gsl_fft_real_wavetable * real;
  gsl_fft_halfcomplex_wavetable * hc;
  gsl_fft_real_workspace * work;

  gsl_spline * workspline1=gsl_spline_alloc (gsl_interp_cspline_periodic, 4*PHIPASIZE+1);
  gsl_interp_accel * workacc1=gsl_interp_accel_alloc ();

  //pions


  for (int k=0;k<4*PHIPASIZE;k++)
    for (int j=0;j<PTASIZE;j++)
      {
	workhorsearr[k]=phiparray[k];
	pispec >> resbuff[j][k];
	CHbuff[j][k]=resbuff[j][k];
	//cout << "rb: " << resbuff[j][k] << endl;
      }

//   for (int j=0;j<PTASIZE;j++)
//     for (int k=0;k<PHIPASIZE;k++)
//       {
// 	workhorsearr[4*PHIPASIZE-k-1]=-phiparray[k]+2*M_PI;
// 	resbuff[j][4*PHIPASIZE-k-1]=resbuff[j][k];
//       }
//   for (int j=0;j<PTASIZE;j++)
//     for (int k=0;k<PHIPASIZE;k++)
//       {
// 	workhorsearr[2*PHIPASIZE-k-1]=M_PI-phiparray[k];
// 	resbuff[j][2*PHIPASIZE-k-1]=resbuff[j][k];
// 	resbuff[j][2*PHIPASIZE+k]=resbuff[j][4*PHIPASIZE-k-1];
// 	workhorsearr[2*PHIPASIZE+k]=M_PI+phiparray[k];
//       }

  

  work = gsl_fft_real_workspace_alloc (4*PHIPASIZE);
  real = gsl_fft_real_wavetable_alloc (4*PHIPASIZE);
  hc = gsl_fft_halfcomplex_wavetable_alloc (4*PHIPASIZE);    
  //first set
  
  

  for (int j=0;j<PTASIZE;j++)
    {
      for (int k=0;k<4*PHIPASIZE;k++)
	{
	  //workhorse[k]=10*cos(2*(2*M_PI/4/PHIPASIZE*k-M_PI));
	  //workhorse[k]=10;
	  //workhorsearr[k]=2*M_PI/4/PHIPASIZE*k-M_PI;
	  workhorse[k]=resbuff[j][k];
	}

      
      workhorsearr[4*PHIPASIZE]=phiparray[0]+2*M_PI;
      workhorse[4*PHIPASIZE]=workhorse[0];
      
      
      gsl_spline_init (workspline1, workhorsearr, workhorse, 4*PHIPASIZE+1);
      

      if (j==0)
       for (int k=0;k<4*PHIPASIZE;k++)
	{
	  //printf("{%f,%f},",workhorsearr[k],workhorse[k]);
	  //printf("{%f,%f},",workhorsearr[k],gsl_spline_eval(workspline1,workhorsearr[k],workacc1));
	}

      for (int k=0;k<4*PHIPASIZE;k++)
      {
	double angle = 2*M_PI/4/PHIPASIZE*k + TRIANGLE;
	double ang;
	if (angle < 2*M_PI) ang = angle;
	else ang = angle - 2*M_PI;
	double angle4 = 2*M_PI/4/PHIPASIZE*k + QUADANGLE;
	double ang4;
	if (angle4 < 2*M_PI) ang4 = angle4;
	else ang4 = angle4 - 2*M_PI;
	workhorse[k]=gsl_spline_eval(workspline1,2*M_PI/4/PHIPASIZE*k,workacc1);
// 	workhorseangle[k]=gsl_spline_eval(workspline1,ang,workacc1);
// 	workhorsequadangle[k]=gsl_spline_eval(workspline1,ang4,workacc1);
      }

      //printf("pt %f wh %f\n",ptarray[j],workhorse[0]);
      gsl_fft_real_transform (workhorse,1 ,4*PHIPASIZE, real, work);
//       gsl_fft_real_transform (workhorseangle,1 ,4*PHIPASIZE, real, work);
//       gsl_fft_real_transform (workhorsequadangle,1 ,4*PHIPASIZE, real, work);
      

      respiv0 << ptarray[j] << "\t";
      respiv0 << workhorse[0]/4/PHIPASIZE;
      respiv0 << "\n";
      respiv2 << ptarray[j]<< "\t";
      respiv2 << workhorse[3]/workhorse[0];
      respiv2 << "\n";
      respiv4 << ptarray[j]<< "\t";
      respiv4 << workhorse[7]/workhorse[0];
      respiv4 << "\n";
      vallpis << ptarray[j] << "\t" //p_t
	      << workhorse[0]/4/PHIPASIZE << "\t" //v_0
	      << workhorse[1]/workhorse[0] << "\t" //<cos(phi)>
	      << -workhorse[2]/workhorse[0] << "\t" //<sin(phi)>
	      << workhorse[3]/workhorse[0] << "\t" //<cos(2 phi)>
	      << -workhorse[4]/workhorse[0] << "\t" //<sin(2 phi)>
	      << workhorse[5]/workhorse[0] << "\t" //<cos(3 phi)>
	      << -workhorse[6]/workhorse[0] << "\t" //<sin(3 phi)>
	      << workhorse[7]/workhorse[0] << "\t" //<cos(4 phi)>
	      << -workhorse[8]/workhorse[0] << "\t" //<sin(4 phi)>
	      << workhorse[9]/workhorse[0] << "\t" //<cos(5 phi)>
	      << -workhorse[10]/workhorse[0] << "\t" //<sin(5 phi)>
	      << workhorse[11]/workhorse[0] << "\t" //<cos(6 phi)>
	      << -workhorse[12]/workhorse[0] << "\t" //<sin(6 phi)>
	      << workhorse[13]/workhorse[0] << "\t" //<cos(7 phi)>
	      << -workhorse[14]/workhorse[0] << "\n"; //<sin(7 phi)>
	      
	      
// 	      << workhorseangle[5]/workhorseangle[0] << "\t"
// 	      << workhorsequadangle[7]/workhorsequadangle[0] << "\n";
//       vsinpis << ptarray[j] << "\t"
// 	      << workhorse[2]/workhorse[0] << "\t"
// 	      << workhorse[4]/workhorse[0] << "\t"
// 	      << workhorseangle[6]/workhorseangle[0] << "\t"
// 	      << workhorsequadangle[8]/workhorsequadangle[0] << "\n";
      
    }

  //kaons

  for (int k=0;k<4*PHIPASIZE;k++)
    for (int j=0;j<PTASIZE;j++)
      {
	workhorsearr[k]=phiparray[k];
	Kspec >> resbuff[j][k];
	CHbuff[j][k]+=resbuff[j][k];
      }

//   for (int j=0;j<PTASIZE;j++)
//     for (int k=0;k<PHIPASIZE;k++)
//       {
// 	workhorsearr[4*PHIPASIZE-k-1]=-phiparray[k]+2*M_PI;
// 	resbuff[j][4*PHIPASIZE-k-1]=resbuff[j][k];
//       }
//   for (int j=0;j<PTASIZE;j++)
//     for (int k=0;k<PHIPASIZE;k++)
//       {
// 	workhorsearr[2*PHIPASIZE-k-1]=M_PI-phiparray[k];
// 	resbuff[j][2*PHIPASIZE-k-1]=resbuff[j][k];
// 	resbuff[j][2*PHIPASIZE+k]=resbuff[j][4*PHIPASIZE-k-1];
// 	workhorsearr[2*PHIPASIZE+k]=M_PI+phiparray[k];
//       }

  

  for (int j=0;j<PTASIZE;j++)
    {
      for (int k=0;k<4*PHIPASIZE;k++)
	workhorse[k]=resbuff[j][k];
      
      workhorsearr[4*PHIPASIZE]=phiparray[0]+2*M_PI;
      workhorse[4*PHIPASIZE]=workhorse[0];
      
      
      gsl_spline_init (workspline1, workhorsearr, workhorse, 4*PHIPASIZE+1);

      for (int k=0;k<4*PHIPASIZE;k++)
      {
	double angle = 2*M_PI/4/PHIPASIZE*k + TRIANGLE;
	double ang;
	if (angle < 2*M_PI) ang = angle;
	else ang = angle - 2*M_PI;
	workhorse[k]=gsl_spline_eval(workspline1,ang,workacc1);
      }

      gsl_fft_real_transform (workhorse,1 ,4*PHIPASIZE, real, work);
      
      for (int k=0;k<4*PHIPASIZE;k++)
	{
	  if (j==0)
	    printf("%i wh %f\n",k,workhorse[k]/4/PHIPASIZE);
	}

      resKv0 << ptarray[j] << "\t";
      resKv0 << workhorse[0]/4/PHIPASIZE;
      resKv0 << "\n";
      resKv2 << ptarray[j]<< "\t";
      resKv2 << workhorse[3]/workhorse[0];
      resKv2 << "\n";
      resKv4 << ptarray[j]<< "\t";
      resKv4 << workhorse[7]/workhorse[0];
      resKv4 << "\n";
      vallKs << ptarray[j] << "\t"
	      << workhorse[0]/4/PHIPASIZE << "\t"
	      << workhorse[1]/workhorse[0] << "\t"
	      << -workhorse[2]/workhorse[0] << "\t"
	      << workhorse[3]/workhorse[0] << "\t"
	      << -workhorse[4]/workhorse[0] << "\t"
	      << workhorse[5]/workhorse[0] << "\t"
	      << -workhorse[6]/workhorse[0] << "\t"
	      << workhorse[7]/workhorse[0] << "\t"
	      << -workhorse[8]/workhorse[0] << "\t"
	      << workhorse[9]/workhorse[0] << "\t"
	      << -workhorse[10]/workhorse[0] << "\t"
	      << workhorse[11]/workhorse[0] << "\t"
	      << -workhorse[12]/workhorse[0] << "\t"
	      << workhorse[13]/workhorse[0] << "\t"
	      << -workhorse[14]/workhorse[0] << "\n";
    }


  //protons

  for (int k=0;k<4*PHIPASIZE;k++)
    for (int j=0;j<PTASIZE;j++)
      {
	workhorsearr[k]=phiparray[k];
	pspec >> resbuff[j][k];
	CHbuff[j][k]+=resbuff[j][k];
      }
  
//   for (int j=0;j<PTASIZE;j++)
//     for (int k=0;k<PHIPASIZE;k++)
//       {
// 	workhorsearr[4*PHIPASIZE-k-1]=-phiparray[k]+2*M_PI;
// 	resbuff[j][4*PHIPASIZE-k-1]=resbuff[j][k];
//       }
//   for (int j=0;j<PTASIZE;j++)
//     for (int k=0;k<PHIPASIZE;k++)
//       {
// 	workhorsearr[2*PHIPASIZE-k-1]=M_PI-phiparray[k];
// 	resbuff[j][2*PHIPASIZE-k-1]=resbuff[j][k];
// 	resbuff[j][2*PHIPASIZE+k]=resbuff[j][4*PHIPASIZE-k-1];
// 	workhorsearr[2*PHIPASIZE+k]=M_PI+phiparray[k];
//       }

  

  for (int j=0;j<PTASIZE;j++)
    {
      for (int k=0;k<4*PHIPASIZE;k++)
	workhorse[k]=resbuff[j][k];
      
      workhorsearr[4*PHIPASIZE]=phiparray[0]+2*M_PI;
      workhorse[4*PHIPASIZE]=workhorse[0];
      
      
      gsl_spline_init (workspline1, workhorsearr, workhorse, 4*PHIPASIZE+1);

      for (int k=0;k<4*PHIPASIZE;k++)
      {
	double angle = 2*M_PI/4/PHIPASIZE*k + TRIANGLE;
	double ang;
	if (angle < 2*M_PI) ang = angle;
	else ang = angle - 2*M_PI;
	workhorse[k]=gsl_spline_eval(workspline1,ang,workacc1);
      }


      gsl_fft_real_transform (workhorse,1 ,4*PHIPASIZE, real, work);
      
      respv0 << ptarray[j] << "\t";
      respv0 << workhorse[0]/4/PHIPASIZE;
      respv0 << "\n";
      respv2 << ptarray[j]<< "\t";
      respv2 << workhorse[3]/workhorse[0];
      respv2 << "\n";
      respv4 << ptarray[j]<< "\t";
      respv4 << workhorse[7]/workhorse[0];
      respv4 << "\n";
      vallps << ptarray[j] << "\t"
	      << workhorse[0]/4/PHIPASIZE << "\t"
	      << workhorse[1]/workhorse[0] << "\t"
	      << -workhorse[2]/workhorse[0] << "\t"
	      << workhorse[3]/workhorse[0] << "\t"
	      << -workhorse[4]/workhorse[0] << "\t"
	      << workhorse[5]/workhorse[0] << "\t"
	      << -workhorse[6]/workhorse[0] << "\t"
	      << workhorse[7]/workhorse[0] << "\t"
	      << -workhorse[8]/workhorse[0] << "\t"
	      << workhorse[9]/workhorse[0] << "\t"
	      << -workhorse[10]/workhorse[0] << "\t"
	      << workhorse[11]/workhorse[0] << "\t"
	      << -workhorse[12]/workhorse[0] << "\t"
	      << workhorse[13]/workhorse[0] << "\t"
	      << -workhorse[14]/workhorse[0] << "\n";
      
    }

  //"charged hadron"
  //which I'll simply define 
  //as sum pi+,pi-,K+,K-,p,pbar divided by two

//   for (int j=0;j<4*PTASIZE;j++)
//     for (int k=0;k<PHIPASIZE;k++)
//       {
// 	workhorsearr[4*PHIPASIZE-k-1]=-phiparray[k]+2*M_PI;
// 	CHbuff[j][4*PHIPASIZE-k-1]=CHbuff[j][k];
//       }
//   for (int j=0;j<PTASIZE;j++)
//     for (int k=0;k<PHIPASIZE;k++)
//       {
// 	workhorsearr[2*PHIPASIZE-k-1]=M_PI-phiparray[k];
// 	CHbuff[j][2*PHIPASIZE-k-1]=CHbuff[j][k];
// 	CHbuff[j][2*PHIPASIZE+k]=CHbuff[j][4*PHIPASIZE-k-1];
// 	workhorsearr[2*PHIPASIZE+k]=M_PI+phiparray[k];
//       }

  

  for (int j=0;j<PTASIZE;j++)
    {
      for (int k=0;k<4*PHIPASIZE;k++)
	workhorse[k]=CHbuff[j][k];
      
      workhorsearr[4*PHIPASIZE]=phiparray[0]+2*M_PI;
      workhorse[4*PHIPASIZE]=workhorse[0];
      
      
      gsl_spline_init (workspline1, workhorsearr, workhorse, 4*PHIPASIZE+1);

      for (int k=0;k<4*PHIPASIZE;k++)
      {
	double angle = 2*M_PI/4/PHIPASIZE*k + TRIANGLE;
	double ang;
	if (angle < 2*M_PI) ang = angle;
	else ang = angle - 2*M_PI;
	workhorse[k]=gsl_spline_eval(workspline1,ang,workacc1);
      }


      gsl_fft_real_transform (workhorse,1 ,4*PHIPASIZE, real, work);
      
      for (int k=0;k<4*PHIPASIZE;k++)
	{
	  if (j==0)
	    printf("%i wh %f\n",k,workhorse[k]/4/PHIPASIZE);
	}

      resHv0 << ptarray[j] << "\t";
      resHv0 << workhorse[0]/4/PHIPASIZE;
      resHv0 << "\n";
      resHv2 << ptarray[j]<< "\t";
      resHv2 << workhorse[3]/workhorse[0];
      resHv2 << "\n";
      resHv4 << ptarray[j]<< "\t";
      resHv4 << workhorse[7]/workhorse[0];
      resHv4 << "\n";
      vallH << ptarray[j] << "\t"
	      << workhorse[0]/4/PHIPASIZE << "\t"
	      << workhorse[1]/workhorse[0] << "\t"
	      << -workhorse[2]/workhorse[0] << "\t"
	      << workhorse[3]/workhorse[0] << "\t"
	      << -workhorse[4]/workhorse[0] << "\t"
	      << workhorse[5]/workhorse[0] << "\t"
	      << -workhorse[6]/workhorse[0] << "\t"
	      << workhorse[7]/workhorse[0] << "\t"
	      << -workhorse[8]/workhorse[0] << "\t"
	      << workhorse[9]/workhorse[0] << "\t"
	      << -workhorse[10]/workhorse[0] << "\t"
	      << workhorse[11]/workhorse[0] << "\t"
	      << -workhorse[12]/workhorse[0] << "\t"
	      << workhorse[13]/workhorse[0] << "\t"
	      << -workhorse[14]/workhorse[0] << "\n";
      
    }




  gsl_spline_free (workspline1);
  gsl_interp_accel_free (workacc1);
  

  gsl_fft_halfcomplex_wavetable_free (hc);
  gsl_fft_real_wavetable_free (real);
  gsl_fft_real_workspace_free (work);

  
}

int main()
{
  
  extern void readParameters(const char*);
  
  readParameters("data/params.txt");  

  pispec.open("data/results/spec_211.dat", ios::in);
  Kspec.open("data/results/spec_321.dat", ios::in);
  pspec.open("data/results/spec_2212.dat", ios::in);
  ptarr.open("data/ptarr.dat", ios::in);
  respiv0.open("data/results/pisv0.dat", ios::out);
  respiv2.open("data/results/pisv2.dat", ios::out);
  respiv4.open("data/results/pisv4.dat", ios::out);
  resKv0.open("data/results/Ksv0.dat", ios::out);
  resKv2.open("data/results/Ksv2.dat", ios::out);
  resKv4.open("data/results/Ksv4.dat", ios::out);
  respv0.open("data/results/psv0.dat", ios::out);
  respv2.open("data/results/psv2.dat", ios::out);
  respv4.open("data/results/psv4.dat", ios::out);
  resHv0.open("data/results/Hv0.dat", ios::out);
  resHv2.open("data/results/Hv2.dat", ios::out);
  resHv4.open("data/results/Hv4.dat", ios::out);
  vallpis.open("data/results/pisvall.dat", ios::out);
  vallKs.open("data/results/Ksvall.dat", ios::out);
  vallps.open("data/results/psvall.dat", ios::out);
  vallH.open("data/results/Hvall.dat", ios::out);
//   vsinpis.open("data/results/pisvsin.dat", ios::out);
//   vsinKs.open("data/results/Ksvsin.dat", ios::out);
//   vsinps.open("data/results/psvsin.dat", ios::out);
//   vsinH.open("data/results/Hvsin.dat", ios::out);
  


  doall();

  pispec.close();
  Kspec.close();
  pspec.close();
  ptarr.close();

  respiv0.close();
  respiv2.close();
  respiv4.close();

  resKv0.close();
  resKv2.close();
  resKv4.close();

  respv0.close();
  respv2.close();
  respv4.close();

  resHv0.close();
  resHv2.close();
  resHv4.close();
  
  vallpis.close();
  vallKs.close();
  vallps.close();
  vallH.close();

  delete [] ptarray;
  delete [] phiparray;

  printf("Done!\n");

  return 0;
}
