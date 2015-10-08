#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <convert.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_multifit.h>

// #include	"reso.h"
// #include	"functions.h"

using namespace std;

int method = 1;

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

//controls value of tau_Pi
double COEFF=3.0;

fstream pttab,ptfile;
const int maxline = 128; // maximum line length used in the buffer for reading
fstream massfile,namesfile,gsfile;
fstream v0file, v2file, v4file, vallfile;//, vsinfile;

int main() 
{
  
  extern void readParameters(const char*);
  
  readParameters("data/params.txt");  
  
  double gsfact=1;
  double tempmass=0;
  double oldmass=0;
  char buffer[maxline];
  int particle=0;
  int i = 0;
  double resbuff[PTASIZE][4*PHIPASIZE]; 
  double ptbuff[PTASIZE];
  double phipbuff[4*PHIPASIZE];
  double weightbuff[PHIPASIZE];
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

  for (int k=0;k<4*PHIPASIZE;k++) phipbuff[k] = k*M_PI/2/PHIPASIZE;
  
  cout << "phipbuff[0] = " << phipbuff[0] << endl;
  cout << "phipbuff[PHIPASIZE-1] = " << phipbuff[PHIPASIZE-1] << endl;
  
  
  massfile.open("pasim.dat", ios::in);
  namesfile.open("pasinames.dat", ios::in);
  gsfile.open("gslist.dat",ios::in);

  pttab.open("data/phipspectra.dat", ios::in);
  ptfile.open("data/ptarr.dat", ios::in);
  
  while (!massfile.eof())
  {
    massfile >> tempmass;
    gsfile >> gsfact;
    namesfile.getline(buffer,maxline,'\n');
    if (method)
    {
// 	  printf("Integrating for %s with mass %f and spin gs=%f\n",buffer,tempmass,gsfact);
	  if (tempmass < 1.0)
	  {
	    
	    for (int k=0;k<4*PHIPASIZE;k++)
	    {
	      for (int j=0;j<PTASIZE;j++)
		{
		  pttab >> resbuff[j][k];
// 		  if (k<PHIPASIZE) workhorsearr[k]=phipbuff[k];
		  workhorsearr[k]=phipbuff[k];
//		  CHbuff[j][k]=resbuff[j][k];
		}
// 	      pttab << "\n";
	    }
	      int j=0;
// 	      cout << "i = " << particle << endl;
	      switch (particle) 
	      {
		case 0:
		  cout << "Case 0 = " << buffer << endl;
		  v0file.open("data/results/preresopisv0.dat", ios::out);
		  v2file.open("data/results/preresopisv2.dat", ios::out);
		  v4file.open("data/results/preresopisv4.dat", ios::out);
		  vallfile.open("data/results/preresopisvall.dat", ios::out);
// 		  vsinfile.open("data/results/preresopisvsin.dat", ios::out);
		  break;
		case 3:
		  cout << "Case 3 = " << buffer << endl;
		  v0file.open("data/results/preresoKsv0.dat", ios::out);
		  v2file.open("data/results/preresoKsv2.dat", ios::out);
		  v4file.open("data/results/preresoKsv4.dat", ios::out);
		  vallfile.open("data/results/preresoKsvall.dat", ios::out);
// 		  vsinfile.open("data/results/preresoKsvsin.dat", ios::out);
		  break;
		case 16:
		  cout << "Case 16 = " << buffer << endl;
		  v0file.open("data/results/preresopsv0.dat", ios::out);
		  v2file.open("data/results/preresopsv2.dat", ios::out);
		  v4file.open("data/results/preresopsv4.dat", ios::out);
		  vallfile.open("data/results/preresopsvall.dat", ios::out);
// 		  vsinfile.open("data/results/preresopsvsin.dat", ios::out);
		  break;
		default:
		  v0file.open("/dev/null", ios::out);
		  v2file.open("/dev/null", ios::out);
		  v4file.open("/dev/null", ios::out);
		  vallfile.open("/dev/null", ios::out);
// 		  vsinfile.open("/dev/null", ios::out);
	      } 
	      //phip-table
	    for (double pt=0.01;pt<PTMAX;pt+=PTMAX/PTASIZE)
	      {
		for(int k=0;k<PHIPASIZE;k++)
		  {
		    ptbuff[j]=pt;
// 		    workhorsearr[4*PHIPASIZE-k-1]=-phipbuff[k]+2*M_PI;
//		    resbuff[j][4*PHIPASIZE-k-1]=resbuff[j][k];
		  }
		j++;
	      }
// 	    for (int j=0;j<PTASIZE;j++)
// 	      for (int k=0;k<PHIPASIZE;k++)
// 	      {
// 		workhorsearr[2*PHIPASIZE-k-1]=M_PI-phipbuff[k];
// //		resbuff[j][2*PHIPASIZE-k-1]=resbuff[j][k];
// //		resbuff[j][2*PHIPASIZE+k]=resbuff[j][4*PHIPASIZE-k-1];
// 		workhorsearr[2*PHIPASIZE+k]=M_PI+phipbuff[k];
// 	      }
	    work = gsl_fft_real_workspace_alloc (4*PHIPASIZE);
	    real = gsl_fft_real_wavetable_alloc (4*PHIPASIZE);
	    hc = gsl_fft_halfcomplex_wavetable_alloc (4*PHIPASIZE);   
	    
	    for (int j=0;j<PTASIZE;j++)
	    {
	      for (int k=0;k<4*PHIPASIZE;k++)
		{
		  //workhorse[k]=10*cos(2*(2*M_PI/4/PHIPASIZE*k-M_PI));
		  //workhorse[k]=10;
		  //workhorsearr[k]=2*M_PI/4/PHIPASIZE*k-M_PI;
		  workhorse[k]=resbuff[j][k];
		}

      
	      workhorsearr[4*PHIPASIZE]=phipbuff[0]+2*M_PI;
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
// 		  workhorseangle[k]=gsl_spline_eval(workspline1,ang,workacc1);
// 		  workhorsequadangle[k]=gsl_spline_eval(workspline1,ang4,workacc1);
		  if (j==0)
		  {
		    //printf("%f wh %f\n",2*M_PI/4/PHIPASIZE*k-M_PI,workhorse[k]);
		    //printf("{%f,%f},",2*M_PI/4/PHIPASIZE*k,workhorse[k]);
		  }
		}
	      //printf("pt %f wh %f\n",ptbuff[j],workhorse[0]);
	      gsl_fft_real_transform (workhorse,1 ,4*PHIPASIZE, real, work);
// 	      gsl_fft_real_transform (workhorseangle,1 ,4*PHIPASIZE, real, work);
// 	      gsl_fft_real_transform (workhorsequadangle,1 ,4*PHIPASIZE, real, work);
	    
	      v0file << ptbuff[j] << "\t";
	      v0file << workhorse[0]/4/PHIPASIZE;
	      v0file << "\n";
	      v2file << ptbuff[j]<< "\t";
	      v2file << workhorse[3]/workhorse[0];
	      v2file << "\n";
	      v4file << ptbuff[j]<< "\t";
	      v4file << workhorse[7]/workhorse[0];
	      v4file << "\n";
	      vallfile << ptbuff[j] << "\t"
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
	    v0file.close();
	    v2file.close();
	    v4file.close();
	    vallfile.close();
	  }
	  else
	    {
	      //don't do anything, just repeat last result
	    }
	  



	  oldmass=tempmass;
      
	  particle++;
	}
      else
	{
// 	printf("testIntegrating for %s with mass %f and spin gs=%f\n",buffer,tempmass,gsfact);
	  if (tempmass < 1.0)
	    {	
	    
	    for (int k=0;k<PHIPASIZE;k++)
	    {
	      for (int j=0;j<PTASIZE;j++)
		{
		  pttab >> resbuff[j][k];
// 		  pttab << resbuff[j][k] << "\t";
		}
// 	      pttab << "\n";
	    }
	    
	      int j=0;
// 	      cout << "i = " << particle << endl;
	      switch (particle) 
	      {
		case 0:
		  cout << "Case 0 = " << buffer << endl;
		  v0file.open("data/results/GLpreresopisv0.dat", ios::out);
		  v2file.open("data/results/GLpreresopisv2.dat", ios::out);
		  v4file.open("data/results/GLpreresopisv4.dat", ios::out);
		  break;
		case 3:
		  cout << "Case 3 = " << buffer << endl;
		  v0file.open("data/results/GLpreresoKsv0.dat", ios::out);
		  v2file.open("data/results/GLpreresoKsv2.dat", ios::out);
		  v4file.open("data/results/GLpreresoKsv4.dat", ios::out);
		  break;
		case 16:
		  cout << "Case 16 = " << buffer << endl;
		  v0file.open("data/results/GLpreresopsv0.dat", ios::out);
		  v2file.open("data/results/GLpreresopsv2.dat", ios::out);
		  v4file.open("data/results/GLpreresopsv4.dat", ios::out);
		  break;
		default:
		  v0file.open("/dev/null", ios::out);
		  v2file.open("/dev/null", ios::out);
		  v4file.open("/dev/null", ios::out);
	      }  
	      for (double pt=0.01;pt<PTMAX;pt+=PTMAX/PTASIZE)
		{
		  double sumv0 = 0.0;
		  double sumv2 = 0.0;
		  double sumv4 = 0.0;
		  v0file << pt << "\t";
		  v2file << pt << "\t";
		  v4file << pt << "\t";
		  for(int k=0;k<PHIPASIZE;k++)
		    {
// 		      ptbuff[j]=pt;
// 		      if (j==99) cout << resbuff[j][k] << endl;
		      sumv0 += weightbuff[k] * resbuff[j][k];
		      sumv2 += weightbuff[k] * resbuff[j][k] * cos(2 * phipbuff[k]);
		      sumv4 += weightbuff[k] * resbuff[j][k] * cos(4 * phipbuff[k]);
		    }
		  sumv2/=sumv0;
		  sumv4/=sumv0;
		  sumv0/=2.0;
		    
		  v0file << sumv0 << "\n";
		  v2file << sumv2 << "\n";
		  v4file << sumv4 << "\n";

		  j++;
		}

		v0file.close();
		v2file.close();
		v4file.close();

	      
	    }
	  else
	    {
	      //don't do anything, just repeat last result
	    }
	  


	  oldmass=tempmass;
      
	  particle++;
    }
  }


  
  massfile.close();
  namesfile.close();
  gsfile.close();
  pttab.close();
  ptfile.close();
  
}
