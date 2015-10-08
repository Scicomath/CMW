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

int PTASIZE;
int PHIPASIZE;

double * ptarray;
double * phiparray;


void doall()
{
  //pions
  

  ptarr >> PTASIZE;
  ptarr >> PHIPASIZE;

  printf("Found %i pt values and %i phip values\n",PTASIZE,PHIPASIZE);

  ptarray=new double[PTASIZE];
  phiparray=new double[PHIPASIZE];
  
  for (int i=0;i<PTASIZE;i++)
    {
      ptarr >> ptarray[i];
    }

  switch (PHIPASIZE) {
  case 2:
    for(int i=0;i<1;i++){
      phiparray[1-i] = 0.25*M_PI*(gaus2x[i] + 1.0);
      phiparray[i] = 0.25*M_PI*(1.0 - gaus2x[i]);
    }
    break;
  case 4:
    for(int i=0;i<2;i++){
      phiparray[3-i] = 0.25*M_PI*(gaus4x[i] + 1.0);
      phiparray[i] = 0.25*M_PI*(1.0 - gaus4x[i]);
    }
    break;
  case 8:
    for(int i=0;i<4;i++){
      phiparray[7-i] = 0.25*M_PI*(gaus8x[i] + 1.0);
      phiparray[i] = 0.25*M_PI*(1.0 - gaus8x[i]);
    }
    break;
  case 10:
    for(int i=0;i<5;i++){
      phiparray[9-i] = 0.25*M_PI*(gaus10x[i] + 1.0);
      phiparray[i] = 0.25*M_PI*(1.0 - gaus10x[i]);
    }
    break;
  case 12:
    for(int i=0;i<6;i++){
      phiparray[11-i] = 0.25*M_PI*(gaus12x[i] + 1.0);
      phiparray[i] = 0.25*M_PI*(1.0 - gaus12x[i]);
    }
    break;
  case 16:
    for(int i=0;i<8;i++){
      phiparray[15-i] = 0.25*M_PI*(gaus16x[i] + 1.0);
      phiparray[i] = 0.25*M_PI*(1.0 - gaus16x[i]);
    }
    break;
  case 20:
    for(int i=0;i<10;i++){
      phiparray[19-i] = 0.25*M_PI*(gaus20x[i] + 1.0);
      phiparray[i] = 0.25*M_PI*(1.0 - gaus20x[i]);
    }
    break;
  case 48:
    for(int i=0;i<24;i++){
      phiparray[47-i] = 0.25*M_PI*(gaus48x[i] + 1.0);
      phiparray[i] = 0.25*M_PI*(1.0 - gaus48x[i]);
    }
    break;
  default:
    printf(" No abscissas for nPhi = %i !\n",PHIPASIZE);
    printf(" GOOD BYE AND HAVE A NICE DAY! \n");
    exit(0);
  }

  double resbuff[PTASIZE][4*PHIPASIZE];
  double CHbuff[PTASIZE][4*PHIPASIZE];

  double workhorsearr[4*PHIPASIZE+1];
  double workhorse[4*PHIPASIZE+1];

  gsl_fft_real_wavetable * real;
  gsl_fft_halfcomplex_wavetable * hc;
  gsl_fft_real_workspace * work;

  gsl_spline * workspline1=gsl_spline_alloc (gsl_interp_cspline_periodic, 4*PHIPASIZE+1);
  gsl_interp_accel * workacc1=gsl_interp_accel_alloc ();

  //pions


  for (int k=0;k<PHIPASIZE;k++)
    for (int j=0;j<PTASIZE;j++)
      {
	workhorsearr[k]=phiparray[k];
	pispec >> resbuff[j][k];
	CHbuff[j][k]=resbuff[j][k];
	//cout << "rb: " << resbuff[j][k] << endl;
      }

  for (int j=0;j<PTASIZE;j++)
    for (int k=0;k<PHIPASIZE;k++)
      {
	workhorsearr[4*PHIPASIZE-k-1]=-phiparray[k]+2*M_PI;
	resbuff[j][4*PHIPASIZE-k-1]=resbuff[j][k];
      }
  for (int j=0;j<PTASIZE;j++)
    for (int k=0;k<PHIPASIZE;k++)
      {
	workhorsearr[2*PHIPASIZE-k-1]=M_PI-phiparray[k];
	resbuff[j][2*PHIPASIZE-k-1]=resbuff[j][k];
	resbuff[j][2*PHIPASIZE+k]=resbuff[j][4*PHIPASIZE-k-1];
	workhorsearr[2*PHIPASIZE+k]=M_PI+phiparray[k];
      }

  

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
	  workhorse[k]=gsl_spline_eval(workspline1,2*M_PI/4/PHIPASIZE*k,workacc1);
	  if (j==0)
	  {
	    //printf("%f wh %f\n",2*M_PI/4/PHIPASIZE*k-M_PI,workhorse[k]);
	    //printf("{%f,%f},",2*M_PI/4/PHIPASIZE*k,workhorse[k]);
	  }
	}
      //printf("pt %f wh %f\n",ptarray[j],workhorse[0]);
      gsl_fft_real_transform (workhorse,1 ,4*PHIPASIZE, real, work);
      

      respiv0 << ptarray[j] << "\t";
      respiv0 << workhorse[0]/4/PHIPASIZE;
      respiv0 << "\n";
      respiv2 << ptarray[j]<< "\t";
      respiv2 << workhorse[3]/workhorse[0];
      respiv2 << "\n";
      respiv4 << ptarray[j]<< "\t";
      respiv4 << workhorse[7]/workhorse[0];
      respiv4 << "\n";
      
    }

  //kaons

  for (int k=0;k<PHIPASIZE;k++)
    for (int j=0;j<PTASIZE;j++)
      {
	workhorsearr[k]=phiparray[k];
	Kspec >> resbuff[j][k];
	CHbuff[j][k]+=resbuff[j][k];
      }

  for (int j=0;j<PTASIZE;j++)
    for (int k=0;k<PHIPASIZE;k++)
      {
	workhorsearr[4*PHIPASIZE-k-1]=-phiparray[k]+2*M_PI;
	resbuff[j][4*PHIPASIZE-k-1]=resbuff[j][k];
      }
  for (int j=0;j<PTASIZE;j++)
    for (int k=0;k<PHIPASIZE;k++)
      {
	workhorsearr[2*PHIPASIZE-k-1]=M_PI-phiparray[k];
	resbuff[j][2*PHIPASIZE-k-1]=resbuff[j][k];
	resbuff[j][2*PHIPASIZE+k]=resbuff[j][4*PHIPASIZE-k-1];
	workhorsearr[2*PHIPASIZE+k]=M_PI+phiparray[k];
      }

  

  for (int j=0;j<PTASIZE;j++)
    {
      for (int k=0;k<4*PHIPASIZE;k++)
	workhorse[k]=resbuff[j][k];
      
      workhorsearr[4*PHIPASIZE]=phiparray[0]+2*M_PI;
      workhorse[4*PHIPASIZE]=workhorse[0];
      
      
      gsl_spline_init (workspline1, workhorsearr, workhorse, 4*PHIPASIZE+1);

      for (int k=0;k<4*PHIPASIZE;k++)
	  workhorse[k]=gsl_spline_eval(workspline1,2*M_PI/4/PHIPASIZE*k,workacc1);


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
      
    }


  //protons

  for (int k=0;k<PHIPASIZE;k++)
    for (int j=0;j<PTASIZE;j++)
      {
	workhorsearr[k]=phiparray[k];
	pspec >> resbuff[j][k];
	CHbuff[j][k]+=resbuff[j][k];
      }
  
  for (int j=0;j<PTASIZE;j++)
    for (int k=0;k<PHIPASIZE;k++)
      {
	workhorsearr[4*PHIPASIZE-k-1]=-phiparray[k]+2*M_PI;
	resbuff[j][4*PHIPASIZE-k-1]=resbuff[j][k];
      }
  for (int j=0;j<PTASIZE;j++)
    for (int k=0;k<PHIPASIZE;k++)
      {
	workhorsearr[2*PHIPASIZE-k-1]=M_PI-phiparray[k];
	resbuff[j][2*PHIPASIZE-k-1]=resbuff[j][k];
	resbuff[j][2*PHIPASIZE+k]=resbuff[j][4*PHIPASIZE-k-1];
	workhorsearr[2*PHIPASIZE+k]=M_PI+phiparray[k];
      }

  

  for (int j=0;j<PTASIZE;j++)
    {
      for (int k=0;k<4*PHIPASIZE;k++)
	workhorse[k]=resbuff[j][k];
      
      workhorsearr[4*PHIPASIZE]=phiparray[0]+2*M_PI;
      workhorse[4*PHIPASIZE]=workhorse[0];
      
      
      gsl_spline_init (workspline1, workhorsearr, workhorse, 4*PHIPASIZE+1);

      for (int k=0;k<4*PHIPASIZE;k++)
	  workhorse[k]=gsl_spline_eval(workspline1,2*M_PI/4/PHIPASIZE*k,workacc1);


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
      
    }

  //"charged hadron"
  //which I'll simply define 
  //as sum pi+,pi-,K+,K-,p,pbar divided by two

  for (int j=0;j<PTASIZE;j++)
    for (int k=0;k<PHIPASIZE;k++)
      {
	workhorsearr[4*PHIPASIZE-k-1]=-phiparray[k]+2*M_PI;
	CHbuff[j][4*PHIPASIZE-k-1]=CHbuff[j][k];
      }
  for (int j=0;j<PTASIZE;j++)
    for (int k=0;k<PHIPASIZE;k++)
      {
	workhorsearr[2*PHIPASIZE-k-1]=M_PI-phiparray[k];
	CHbuff[j][2*PHIPASIZE-k-1]=CHbuff[j][k];
	CHbuff[j][2*PHIPASIZE+k]=CHbuff[j][4*PHIPASIZE-k-1];
	workhorsearr[2*PHIPASIZE+k]=M_PI+phiparray[k];
      }

  

  for (int j=0;j<PTASIZE;j++)
    {
      for (int k=0;k<4*PHIPASIZE;k++)
	workhorse[k]=CHbuff[j][k];
      
      workhorsearr[4*PHIPASIZE]=phiparray[0]+2*M_PI;
      workhorse[4*PHIPASIZE]=workhorse[0];
      
      
      gsl_spline_init (workspline1, workhorsearr, workhorse, 4*PHIPASIZE+1);

      for (int k=0;k<4*PHIPASIZE;k++)
	  workhorse[k]=gsl_spline_eval(workspline1,2*M_PI/4/PHIPASIZE*k,workacc1);


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
      
    }




  gsl_spline_free (workspline1);
  gsl_interp_accel_free (workacc1);
  

  gsl_fft_halfcomplex_wavetable_free (hc);
  gsl_fft_real_wavetable_free (real);
  gsl_fft_real_workspace_free (work);

  
}

int main()
{

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

  delete [] ptarray;
  delete [] phiparray;

  printf("Done!\n");

  return 0;
}
