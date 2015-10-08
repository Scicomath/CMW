/*
This is a variant of the convert routine that does not
split the freezeout surfact into to regions, but proceeds 
with the first parametrization. It may be numerically
more stable for some cases, but still isn't great.
Best is to run freezeout with both convert routines
and check that the results are similar.
*/

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


int probon=0;

using namespace std;


int NUMT=8;
long int STEPS=4000,UPDATE=100,SNAPUPDATE=1000;
double AT=0.05,EPS=0.001,B=0.0;
double ETAOS=0.3;
double TSTART=0.5,TF=0.1,TINIT=1.0;
double IC;
int PTASIZE=100;
int PHIPASIZE=4;
const double INTACC=1.e-2;
const int INTSPACE=2000;

//controls value of tau_Pi
double COEFF=3.0;


double L1COEF=2.0;
double L2COEF=0.0;
int FREEZE=1;
double PTMAX=4.2, TRIEPS, TRIANGLE;

double *xp,*yp,*ux,*uy,*pixx,*pixy,*piyy,*phi;
double *taus;


//splines
gsl_spline **xspline;
gsl_interp_accel **xacc;
gsl_spline **yspline;
gsl_interp_accel **yacc;
gsl_spline **uxspline;
gsl_interp_accel **uxacc;
gsl_spline **uyspline;
gsl_interp_accel **uyacc;
gsl_spline **pixxspline;
gsl_interp_accel **pixxacc;
gsl_spline **pixyspline;
gsl_interp_accel **pixyacc;
gsl_spline **piyyspline;
gsl_interp_accel **piyyacc;
gsl_interp_accel *tsacc; 
gsl_spline *tsspline;
gsl_interp_accel **dtxacc; 
gsl_spline **dtxspline;
gsl_interp_accel **dtyacc; 
gsl_spline **dtyspline;


//splines for last 10%
gsl_spline **tspline;
gsl_interp_accel **tacc;
gsl_spline **ouxspline;
gsl_interp_accel **ouxacc;
gsl_spline **ouyspline;
gsl_interp_accel **ouyacc;
gsl_spline **opixxspline;
gsl_interp_accel **opixxacc;
gsl_spline **opixyspline;
gsl_interp_accel **opixyacc;
gsl_spline **opiyyspline;
gsl_interp_accel **opiyyacc;
gsl_spline **dxtspline;
gsl_interp_accel **dxtacc;


struct dummy_params
{
  gsl_spline * myspline;
  gsl_interp_accel * myacc;
};




gsl_integration_workspace * w = gsl_integration_workspace_alloc (INTSPACE);

//const int LIMITS=1024;
const int LIMITS=5200;
const int MAXL=1024;
const int POLYORD=10;
const int maxline = 128; // maximum line length used in the buffer for reading

//const int PTASIZE=100;  // how many grid points in PT
//const int PHIPASIZE=12;   // how many grid points in PHIP

long int numset[LIMITS];
int totalnum;
long int subnum;
int middle,numpoints;

double stepper;

int * counterarr;
double *boundarr;


// input files
fstream freeze_out,dummy;

fstream massfile,namesfile,gsfile;

// output files
fstream pttab,ptfile;

//Find out how many sets
void countsets()
{
  char buffer[1024];
  int limiter=0;
  totalnum=0;
  subnum=-1;
  //determine how many sets
  while (!freeze_out.eof())
    {
      freeze_out.getline(buffer,1024,'\n');
      subnum++;
      if (buffer[0]=='T')
	{
	  subnum--;
	  //printf("line= %s\n",buffer);
	  //printf("found at %i\n",limiter);
	  numset[totalnum]=limiter;
	  //printf("got %i\n",numset[totalnum]);
	  totalnum++;
	  if (totalnum>LIMITS)
	    {
	      printf("More than LIMITS sets. Aborting\n");
	      break;
	    }
	}
      limiter++;
    }
  printf("Found %i data lines \n",subnum);

  //totalnum--;

  numpoints=2*(totalnum);

  //make it odd
  if (numpoints%2==0)
    numpoints++;

  middle=(numpoints-1)/2;

  if (middle<5)
    {
      numpoints*=4;
      numpoints++;
      middle=(numpoints-1)/2;
    }


}

//allocate Memory
void allocMem()
{
  xp=new double[subnum];
  yp=new double[subnum];
  phi=new double[subnum];
  ux=new double[subnum];
  uy=new double[subnum];
  pixx=new double[subnum];
  pixy=new double[subnum];
  piyy=new double[subnum];
  taus=new double[totalnum];

  xspline=new gsl_spline*[totalnum];
  xacc=new gsl_interp_accel*[totalnum];
  yspline=new gsl_spline*[totalnum];
  yacc=new gsl_interp_accel*[totalnum];
  uxspline=new gsl_spline*[totalnum];
  uxacc=new gsl_interp_accel*[totalnum];
  uyspline=new gsl_spline*[totalnum];
  uyacc=new gsl_interp_accel*[totalnum];
  pixxspline=new gsl_spline*[totalnum];
  pixxacc=new gsl_interp_accel*[totalnum];
  pixyspline=new gsl_spline*[totalnum];
  pixyacc=new gsl_interp_accel*[totalnum];
  piyyspline=new gsl_spline*[totalnum];
  piyyacc=new gsl_interp_accel*[totalnum];

  dtxspline=new gsl_spline*[totalnum];
  dtxacc=new gsl_interp_accel*[totalnum];
  dtyspline=new gsl_spline*[totalnum];
  dtyacc=new gsl_interp_accel*[totalnum];
  //dtxacc=gsl_interp_accel_alloc ();
  //dtyspline=gsl_spline_alloc (gsl_interp_cspline, 3);
  // dtyacc=gsl_interp_accel_alloc ();


  tspline=new gsl_spline*[numpoints];
  tacc=new gsl_interp_accel*[numpoints];
  ouxspline=new gsl_spline*[numpoints];
  ouxacc=new gsl_interp_accel*[numpoints];
  ouyspline=new gsl_spline*[numpoints];
  ouyacc=new gsl_interp_accel*[numpoints];
  opixxspline=new gsl_spline*[numpoints];
  opixxacc=new gsl_interp_accel*[numpoints];
  opixyspline=new gsl_spline*[numpoints];
  opixyacc=new gsl_interp_accel*[numpoints];
  opiyyspline=new gsl_spline*[numpoints];
  opiyyacc=new gsl_interp_accel*[numpoints];
  dxtspline=new gsl_spline*[numpoints];
  dxtacc=new gsl_interp_accel*[numpoints];


  counterarr=new int[numpoints];
  boundarr=new double[numpoints];

  tsspline=gsl_spline_alloc (gsl_interp_cspline, totalnum-1);
  tsacc=gsl_interp_accel_alloc ();

  //dtxspline=gsl_spline_alloc (gsl_interp_cspline, 3);
  //dtxacc=gsl_interp_accel_alloc ();
  //dtyspline=gsl_spline_alloc (gsl_interp_cspline, 3);
  // dtyacc=gsl_interp_accel_alloc ();
}


void freeMem()
{
  
  for (int i=0;i<totalnum;i++)
    {
      gsl_spline_free (xspline[i]);
      gsl_interp_accel_free (xacc[i]);
      gsl_spline_free (yspline[i]);
      gsl_interp_accel_free (yacc[i]);
      gsl_spline_free (uxspline[i]);
      gsl_interp_accel_free (uxacc[i]);
      gsl_spline_free (uyspline[i]);
      gsl_interp_accel_free (uyacc[i]);
      gsl_spline_free (pixxspline[i]);
      gsl_interp_accel_free (pixxacc[i]);
      gsl_spline_free (pixyspline[i]);
      gsl_interp_accel_free (pixyacc[i]);
      gsl_spline_free (piyyspline[i]);
      gsl_interp_accel_free (piyyacc[i]);
      gsl_spline_free (dtxspline[i]);
      gsl_interp_accel_free (dtxacc[i]);
      gsl_spline_free (dtyspline[i]);
      gsl_interp_accel_free (dtyacc[i]);
    }

  for (int i=0;i<numpoints;i++)
    {
      if (counterarr[i]>2)
	{
	  gsl_spline_free (tspline[i]);
	  gsl_interp_accel_free (tacc[i]);
	  gsl_spline_free (ouxspline[i]);
	  gsl_interp_accel_free (ouxacc[i]);
	  gsl_spline_free (ouyspline[i]);
	  gsl_interp_accel_free (ouyacc[i]);
	  gsl_spline_free (opixxspline[i]);
	  gsl_interp_accel_free (opixxacc[i]);
	  gsl_spline_free (opixyspline[i]);
	  gsl_interp_accel_free (opixyacc[i]);
	  gsl_spline_free (opiyyspline[i]);
	  gsl_interp_accel_free (opiyyacc[i]);
	  gsl_spline_free (dxtspline[i]);
	  gsl_interp_accel_free (dxtacc[i]);
	}
    }


  delete [] counterarr;
  delete [] boundarr;

  gsl_spline_free(tsspline);
  gsl_interp_accel_free(tsacc);

  gsl_integration_workspace_free (w);
}

//Get data
int readsets()
{
  double temp;

  //freeze_out.seekg(0);
  int dumsize=0;

  long int pos=0;
  long int setpos=0;
  char d[128];

  char obda;

  while((!dummy.eof())&&(setpos<totalnum))
    {
      if (pos!=(numset[setpos]-setpos))
	{
	  /*
	  dummy >> temp;
	  
	  if (pos%7==0)
	    {
	      cout << endl;
	      //length=dummy.tellg();
	      //printf("pos = %i\n",length);
	      
	    }
	  printf("%f \t",temp);
	  pos++;
	  */
	  dummy >> xp[pos];
	  dummy >> yp[pos];
	  phi[pos]=atan2(yp[pos],xp[pos]);
	  //take advantage of full 2Pi information:
	  if (phi[pos]<0) phi[pos]+=2*M_PI;
	  dummy >> ux[pos];
	  dummy >> uy[pos];
	  dummy >> pixx[pos];
	  dummy >> pixy[pos];
	  dummy >> piyy[pos];
	  // if ((yp[pos]<0)&&(xp[pos]==0)&&(setpos>switcher))
	  //  printf("{%i,%f},\n",setpos,yp[pos]);
	  //printf("%f %f %f %f %f %f %f \t %f\n",xp[pos],yp[pos],ux[pos],uy[pos],pixx[pos],pixy[pos],piyy[pos],phi[pos]);
	  pos++;
	}
      else
	{
	  dummy >> obda;
	  dummy >> obda;
	  dummy >> obda;
	  dummy >> obda;
	  dummy >> temp;
	  taus[setpos]=temp;
	  //printf("END for t= %f \n",taus[setpos]);
	  setpos++;
	  //pos+=7;
	}
    }
  cout << endl;

  return 0;
}

//Smoothing data -- important for corse lattices
void smoothdata(double *datax,double *data1y,double *data2y, int L)
{
  gsl_spline * workspline1=gsl_spline_alloc (gsl_interp_cspline_periodic, L);
  gsl_interp_accel * workacc1=gsl_interp_accel_alloc ();
  gsl_spline * workspline2=gsl_spline_alloc (gsl_interp_cspline_periodic, L);
  gsl_interp_accel * workacc2=gsl_interp_accel_alloc ();

  gsl_spline_init (workspline1, datax, data1y, L);
  gsl_spline_init (workspline2, datax, data2y, L);
  
  //for (int i=0;i<L;i++)
  //  {
  //    printf("{%f,%f},\n",datax[i],datay[i]);
  //  }

  for (int i=0;i<(L-1);i++)
    {
      double temp=i*2*M_PI/(L-1);
      data1y[i]=gsl_spline_eval(workspline1,temp,workacc1);
      data2y[i]=gsl_spline_eval(workspline2,temp,workacc2);
      //datax[i]=2*M_PI/L*i;
      //datay[i]=sin(datax[i]);
    }

  gsl_fft_real_wavetable * real;
  gsl_fft_halfcomplex_wavetable * hc;
  gsl_fft_real_workspace * work;


  work = gsl_fft_real_workspace_alloc (L-1);
  real = gsl_fft_real_wavetable_alloc (L-1);
  hc = gsl_fft_halfcomplex_wavetable_alloc (L-1);    
  //first set
  
  gsl_fft_real_transform (data1y, 1, L-1, real, work);

  //smoothing:
  //remove constant (should be centered)
  data1y[0]=0;
  
  //remove higher orders
  //if (L>15)
  //  for (int i=cut;i<L-1;i++)
  //   {
  //	data1y[i]=0;
  ////printf("%i %f",i,a[i]);
  //  }

  //remove sin and even cos parts
  for (int i=1;i<L-1;i++)
  {
    //printf("smoother ux %i %f\n",i,data1y[i]);
    if (probon==1)
      printf("%i data1y=%f",i, data1y[i]);
    if (((i-1)%4)!=0)
      data1y[i]=0;
    if (probon==1)
      printf("vs %f\n",data1y[i]);
  }

  gsl_fft_halfcomplex_inverse (data1y, 1, L-1, hc, work);
  
  //second set
  
  gsl_fft_real_transform (data2y, 1, L-1, real, work);

  data2y[0]=0;
  for (int i=1;i<L-1;i++)
  {
    //printf("smoother uy %i %f\n",i,data2y[i]);
    if (probon==1)
      printf("%i data2y=%f",i, data2y[i]);
    if (((i-2)%4)!=0)
     data2y[i]=0;
    if (probon==1)
      printf("vs %f\n",data2y[i]);
  }

  //if (L>15)
  //    for (int i=cut;i<L-1;i++)
  //  data2y[i]=0;
 
  gsl_fft_halfcomplex_inverse (data2y, 1, L-1, hc, work);


  //clean up
  gsl_fft_halfcomplex_wavetable_free (hc);
  gsl_fft_real_wavetable_free (real);
  gsl_fft_real_workspace_free (work);

  gsl_spline_free (workspline1);
  gsl_interp_accel_free (workacc1);
  gsl_spline_free (workspline2);
  gsl_interp_accel_free (workacc2);

}

//Interpolation routines

void interpolate1()
{
  const gsl_interp_type *t = gsl_interp_cspline_periodic; 

  for (int i=0;i<totalnum;i++)
    {
      int length=numset[0]+1;
      int start=0;
      if (i!=0)
	{
	  length=numset[i]-numset[i-1];
	  start=numset[i-1]-i+1;
	}
      //cout << "length=" << length;
      xacc[i] = gsl_interp_accel_alloc ();
      yacc[i] = gsl_interp_accel_alloc ();
      uxacc[i] = gsl_interp_accel_alloc ();
      uyacc[i] = gsl_interp_accel_alloc ();
      pixxacc[i] = gsl_interp_accel_alloc ();
      pixyacc[i] = gsl_interp_accel_alloc ();
      piyyacc[i] = gsl_interp_accel_alloc ();
      xspline[i] = gsl_spline_alloc (t, length);
      yspline[i] = gsl_spline_alloc (t, length);
      uxspline[i] = gsl_spline_alloc (t, length);
      uyspline[i] = gsl_spline_alloc (t, length);
      pixxspline[i] = gsl_spline_alloc (t, length);
      pixyspline[i] = gsl_spline_alloc (t, length);
      piyyspline[i] = gsl_spline_alloc (t, length);
      
      double *tempx,*tempy,*tempux,*tempuy,*temppixx,*temppixy,*temppiyy;
      tempx=new double[length];
      tempy=new double[length];
      tempux=new double[length];
      tempuy=new double[length];
      temppixx=new double[length];
      temppixy=new double[length];
      temppiyy=new double[length];
      double *tempp;
      tempp=new double[length];
      

      for (int j=0;j<length-1;j++)
	{

	  tempp[j]=phi[start+j];
	  if (j>0)
	    {
	      if (tempp[j]>tempp[j-1])
		{
		  tempx[j]=xp[start+j];
		  tempy[j]=yp[start+j];
		  tempux[j]=ux[start+j];
		  tempuy[j]=uy[start+j];
		  temppixx[j]=pixx[start+j];
		  temppixy[j]=pixy[start+j];
		  temppiyy[j]=piyy[start+j];
		}
	      else
		{
		  /*
		  tempp[j]=tempp[j-1];
		  tempp[j-1]=phi[start+j];
		  
		  tempx[j]=tempx[j-1];
		  tempx[j-1]=xp[start+j];
		  tempy[j]=tempy[j-1];
		  tempy[j-1]=yp[start+j];

		  tempux[j]=tempux[j-1];
		  tempux[j-1]=ux[start+j];
		  tempuy[j]=tempuy[j-1];
		  tempuy[j-1]=uy[start+j];

		  temppixx[j]=temppixx[j-1];
		  temppixx[j-1]=pixx[start+j];
		  temppixy[j]=temppixy[j-1];
		  temppixy[j-1]=pixy[start+j];
		  temppiyy[j]=temppiyy[j-1];
		  temppiyy[j-1]=piyy[start+j];
		  */

		  printf("WARNING: non-monotonic behaviour at %i of %i\n",i,totalnum);
		}
	    }
	  else
	    {
	      tempx[j]=xp[start+j];
	      tempy[j]=yp[start+j];
	      tempux[j]=ux[start+j];
	      tempuy[j]=uy[start+j];
	      temppixx[j]=pixx[start+j];
	      temppixy[j]=pixy[start+j];
	      temppiyy[j]=piyy[start+j];
	    }

	  
	  // if (i>switcher)
	  // printf("%i %i %f vs %f\n",i,j,ux[start+j]/uy[start+j],xp[start+j]/yp[start+j]);
	  //if ((tempx[j]==0)&&(tempy[j]<0)&&(i>switcher))
	  // printf("{%f,%f},\n",taus[i],tempy[j]);
	  
	  
	  
	  
	}
       

      /*
      for (int j=0;j<length-1;j++)
	{
	  if (i==totalnum-2)
	printf("tot2 %f %f at %f\n",tempx[j],tempy[j],tempp[j]);
	  if (i==totalnum-1)
	    printf("tot1 %i %f %f at %f\n",i,tempx[j],tempy[j],tempp[j]);
	}
      */
      
      tempx[length-1]=xp[start];
      tempy[length-1]=yp[start]; 
      tempux[length-1]=ux[start];
      tempuy[length-1]=uy[start];
      temppixx[length-1]=pixx[start];
      temppixy[length-1]=pixy[start];
      temppiyy[length-1]=piyy[start];
      tempp[length-1]=2*M_PI;
 
      
      //smoothdata(tempp,tempx,tempy,length);
      //if (i>switcher)
      //	probon=1;
      //smoothdata(tempp,tempux,tempuy,length);
      //probon=0;
      //smoothdata(tempp,tempx,tempy,length,(int) 3);

      /*
      double ddd[length];

      for (int ll=0;ll<(length-1);ll++)
	{
	  double temp=ll*2*M_PI/(length-1);
	  ddd[ll]=temp;
	}
      


      //phi-splines
      gsl_spline_init (uxspline[i], ddd, tempux, length);
      gsl_spline_init (uyspline[i], ddd, tempuy, length);
      */

      gsl_spline_init (uxspline[i], tempp, tempux, length);
      gsl_spline_init (uyspline[i], tempp, tempuy, length);
      gsl_spline_init (pixxspline[i], tempp, temppixx, length);
      gsl_spline_init (pixyspline[i], tempp, temppixy, length);
      gsl_spline_init (piyyspline[i], tempp, temppiyy, length);

      
 
 
      tempx[length-1]=tempx[0];
      tempy[length-1]=tempy[0]; 
      tempp[length-1]=2*M_PI;

      //if (i==0)
      //	{
      //	  for (int k=0;k<length;k++)
      //	    {
      //	      printf("{%f,%f},\n",tempp[k],tempx[k]);
      //	    }
      //	}
      //smoothdata(tempx,length-1);
      //smoothdata(tempy,length-1);

      gsl_spline_init (xspline[i], tempp, tempx, length);
      gsl_spline_init (yspline[i], tempp, tempy, length);

      /*      if (i==totalnum-1)
	{
	  for (double xi=0;xi<2*M_PI;xi+=0.1)
	    {
	      printf("%f %f %f\n",xi,gsl_spline_eval (xspline[i],xi,xacc[i]),gsl_spline_eval (yspline[i],xi,xacc[i]));
	    }
	}
      */

       //      if (i>switcher)
      //	printf("{%f,%f},\n",taus[i],gsl_spline_eval (yspline[i], 3*M_PI/2., yacc[i]));

      //if ((tempx[j]==0)&&(tempy[j]<0)&&(i>switcher))
	  // printf("{%f,%f},\n",taus[i],tempy[j]);


      /*   
      double yi;
      for (double xi=0;xi<M_PI;xi+=0.1)
        {
	  //yi = gsl_spline_eval (xspline[i], xi, xacc[i]);
	  yi = gsl_spline_eval (uxspline[i], xi, uxacc[i]);
	  //if (i==4000)
	  //{
	      //printf ("%g %g\n", yi,gsl_spline_eval (uxspline[i], 2*M_PI-xi, uxacc[i]));
	      //printf ("%g %g\n", xi,yi);
	      // }
	      if (yi-gsl_spline_eval (uxspline[i], 2*M_PI-xi, uxacc[i])>1.e-2)
		printf ("found at i=%i xi=%f %g %g\n", i,xi,yi,gsl_spline_eval (uxspline[i], 2*M_PI-xi, uxacc[i]));
	      
	}
      */

      delete [] tempx;
      delete [] tempy;
      delete [] tempux;
      delete [] tempuy;
      delete [] temppixx;
      delete [] temppixy;
      delete [] temppiyy;
      delete [] tempp;
    
    }

  //create dtx, dty

  double *tempx,*tempy,*tempp;
  double *datax,*datay;
  tempx=new double[totalnum];
  tempy=new double[totalnum];
  tempp=new double[101];
  datax=new double[101];
  datay=new double[101];

  double bigblockx[100][totalnum];
  double bigblocky[100][totalnum];
  

  gsl_spline * workspline1=gsl_spline_alloc (gsl_interp_cspline, totalnum);
  gsl_interp_accel * workacc1=gsl_interp_accel_alloc ();
  gsl_spline * workspline2=gsl_spline_alloc (gsl_interp_cspline, totalnum);
  gsl_interp_accel * workacc2=gsl_interp_accel_alloc ();
      


  for (int j=0;j<100;j++)
    {
      tempp[j]=2*M_PI/100*j;
      for (int i=0;i<totalnum;i++)
	{
	  tempx[i]=gsl_spline_eval(xspline[i],tempp[j],xacc[i]);
	  tempy[i]=gsl_spline_eval(yspline[i],tempp[j],yacc[i]);
	}
      
      
      gsl_spline_init (workspline1, taus, tempx, totalnum);
      gsl_spline_init (workspline2, taus, tempy, totalnum);
  
      for (int i=0;i<totalnum;i++)
	{
	  bigblockx[j][i]=gsl_spline_eval_deriv(workspline1,taus[i],workacc1);
	  bigblocky[j][i]=gsl_spline_eval_deriv(workspline2,taus[i],workacc2);
	  //bigblockx[j][i]=gsl_spline_eval(workspline1,taus[i],workacc1);
	  //bigblocky[j][i]=gsl_spline_eval(workspline2,taus[i],workacc2);
	}

      //printf("testing %f %f\n",bigblockx[j][0],gsl_spline_eval(xspline[0],tempp[j],xacc[0]));
    }
  
  for (int i=0;i<totalnum;i++)
    {
      for (int j=0;j<100;j++)
	{
	  datax[j]=bigblockx[j][i];
	  datay[j]=bigblocky[j][i];
	}

      datax[100]=datax[0];
      datay[100]=datay[0];
      tempp[100]=2*M_PI;

       
       dtxspline[i] = gsl_spline_alloc (t, 101);      
       dtxacc[i] = gsl_interp_accel_alloc ();
       dtyspline[i] = gsl_spline_alloc (t, 101);      
       dtyacc[i] = gsl_interp_accel_alloc ();


       gsl_spline_init (dtxspline[i], tempp, datax, 101);
       gsl_spline_init (dtyspline[i], tempp, datay, 101);
    }

  //for (int j=0;j<100;j++)
  // {
  ///   printf("t2 %f %f\n",gsl_spline_eval(dtxspline[0],2*M_PI/100*j,dtxacc[0]),gsl_spline_eval(xspline[0],2*M_PI/100*j,xacc[0]));
    //}
  
  gsl_spline_free (workspline1);
  gsl_spline_free (workspline2);

  gsl_interp_accel_free(workacc1);
  gsl_interp_accel_free(workacc2);
    


  delete [] tempx;
  delete [] tempy;
  delete [] tempp; 

}


//gsl_interp_accel * wac=gsl_interp_accel_alloc ();


double rootfunction(double x,void *params)
{
  struct dummy_params *point= (struct dummy_params *) params;

  gsl_spline * ws= point->myspline;
  gsl_interp_accel * wac = point->myacc;

  double temp;
  temp=gsl_spline_eval (ws, x, wac);

  return temp;
}



double getit(double x,double y,double tt,gsl_spline ** ws,gsl_interp_accel ** wsacc)
{

  double tphi=atan2(y,x);
  if (tphi<0)
    tphi+=2*M_PI;

  //int thisnum=6*(totalnum-switcher)/4;
  //if (totalnum-thisnum<0)
  // thisnum=totalnum;

  int thisnum=totalnum;
  double datax[thisnum];
  double datay[thisnum];

  
  gsl_spline * workspline=gsl_spline_alloc (gsl_interp_cspline, thisnum);
  
  for (int i=0;i<thisnum;i++)
    {
      datax[i]=taus[totalnum-thisnum+i];
      datay[i]=gsl_spline_eval (ws[totalnum-thisnum+i], tphi, wsacc[totalnum-thisnum+i]);
      //double hphi=atan2(-y,x);
      //if (probon==1)
      //	printf("inside %f %f\n",datay[i],gsl_spline_eval (ws[totalnum-thisnum+i], hphi, wsacc[totalnum-thisnum+i]));
    }

  gsl_spline_init (workspline, datax, datay, thisnum);
  gsl_interp_accel * wac=gsl_interp_accel_alloc ();

  double temp=gsl_spline_eval (workspline, tt,wac);
  

  gsl_spline_free (workspline);
  gsl_interp_accel_free(wac);

  return temp;
}






double allintegrand(double tphi,void * params)
{
  double *par= (double *) params;
  //cout << "par1= " << par[0] << "\t";
  //cout << "par2= " << par[1] << "\n";
  //par[0]=Temperature
  //par[1]=particle rest-mass
  //par[2]=p_T
  //par[3]=phi_p
  //par[4]=setnumber eqiv tau
  double px=par[2]*cos(par[3]);
  double py=par[2]*sin(par[3]);
  double mt=sqrt(par[2]*par[2]+par[1]*par[1]);
  int thisset=(int) par[4];

  double mux=gsl_spline_eval (uxspline[thisset], tphi, uxacc[thisset]);
  double muy=gsl_spline_eval (uyspline[thisset], tphi, uyacc[thisset]);

  double ut=sqrt(1+mux*mux+muy*muy);

  double dpmx=gsl_spline_eval_deriv (xspline[thisset], tphi, xacc[thisset]);
  double dpmy=gsl_spline_eval_deriv (yspline[thisset], tphi, yacc[thisset]);
  double dtmx,dtmy;
  //nice but slow?
  
  //getdtsplines(thisset,phi);

  dtmx=gsl_spline_eval(dtxspline[thisset],tphi,dtxacc[thisset]);
  dtmy=gsl_spline_eval(dtyspline[thisset],tphi,dtyacc[thisset]);
  
  //printf("hmm %f %f\n",dtmx,gsl_spline_eval(xspline[thisset],tphi,xacc[thisset]));

  //dtmx=(gsl_spline_eval (xspline[thisset+1], tphi, xacc[thisset+1])-gsl_spline_eval (xspline[thisset], tphi, xacc[thisset]))/(taus[thisset+1]-taus[thisset]);
  //dtmy=(gsl_spline_eval (yspline[thisset+1], tphi, yacc[thisset+1])-gsl_spline_eval (yspline[thisset], tphi, yacc[thisset]))/(taus[thisset+1]-taus[thisset]);
  
  //printf("dtmx=%f \t dtmy=%f\t dpmx=%f\t dpmy=%f\n",dtmx,dtmy,dpmx,dpmy);

  double d0k0=gsl_sf_bessel_K0(mt/par[0]*ut);
  double d1k0=gsl_sf_bessel_K1(mt/par[0]*ut);


  double result;

  double temp1,temp2;
  
  //printf("get %f with dtmy=%f and dtmx=%f\n",(dtmx*dpmy-dtmy*dpmx),dtmy/sin(tphi),dtmx/cos(tphi));
  //printf("{%f,%f},",tphi,gsl_spline_eval (xspline[thisset], tphi, xacc[thisset]));
  //printf("{%f,%f},",tphi,gsl_spline_eval (xspline[thisset+1], tphi, xacc[thisset+1]));

  //temp=1;
  temp1=(dtmx*dpmy-dtmy*dpmx)*mt;
  temp2=px*dpmy-py*dpmx;
  

  //if eta is non-negligible, take visc effects into account:
  if (ETAOS>0.001)
    {
 
      double mpixx=gsl_spline_eval (pixxspline[thisset], tphi, pixxacc[thisset])/(TF*TF*2);
      double mpixy=gsl_spline_eval (pixyspline[thisset], tphi, pixyacc[thisset])/(TF*TF*2);
      double mpiyy=gsl_spline_eval (piyyspline[thisset], tphi, piyyacc[thisset])/(TF*TF*2);
      
      double vx=mux/ut;
      double vy=muy/ut;
      double mpitt=vx*vx*mpixx+2*vx*vy*mpixy+vy*vy*mpiyy;
      double mpiee=mpixx+mpiyy-mpitt;
      double mpitx=vx*mpixx+vy*mpixy;
      double mpity=vx*mpixy+vy*mpiyy;
      
      double one=1+mt*mt*mpiee+px*px*mpixx+py*py*mpiyy+2*px*py*mpixy;
      double two=-2*mt*(px*mpitx+py*mpity);
      double three=mt*mt*(mpitt-mpiee);
           
      double d2k0=(gsl_sf_bessel_Kn(2,mt/par[0]*ut)+d0k0)/2;
      double d3k0=(gsl_sf_bessel_Kn(3,mt/par[0]*ut)+3*d1k0)/4;


      temp1*=(one*d1k0+two*d2k0+three*d3k0);
      temp2*=(one*d0k0+two*d1k0+three*d2k0);

      result=temp1-temp2;
    }
  else
    {
      temp1*=d1k0;
      temp2*=d0k0;
      //printf("tphi: %f temp1 %f temp2 %f\n",tphi,temp1,temp2);
      //printf("first int %f second int %f\n",firstintegrand(tphi,params),secondintegrand(tphi,params));
      result=temp1-temp2;
    }

  result*=exp((px*mux+py*muy)/par[0]);
  //printf("result %f comp %f\n",result,firstintegrand(tphi,params)-secondintegrand(tphi,params));
  return result;

}



double ointegrate1(double T,double m0,double pt, double phip, int thisset)
{
  double result,error;
  size_t neval;
  double parameters[5];

  

  parameters[0]=T;
  parameters[1]=m0;
  parameters[2]=pt;
  parameters[3]=phip;
  parameters[4]=(double) thisset;

  gsl_function F;
  F.function = &allintegrand;
  F.params = &parameters;
  
  

  gsl_error_handler_t * new_handler;

  new_handler=gsl_set_error_handler_off();
    

  int bad=3;

  

  int code=gsl_integration_qag(&F,0,2*M_PI,1e-10,INTACC,INTSPACE,3,w,&result,&error);


  while(code==GSL_EROUND)
    {
      bad++;
      //printf("Roundoff error, badness %i, set %i\n",bad-1,thisset);
      code=gsl_integration_qag(&F,0,2*M_PI,1e-10,INTACC,INTSPACE,bad,w,&result,&error);
      
      if (bad==7)
	{
	  printf("I1: Unrecoverable roundoff error detected. Aborting\n");
	  printf("result=%.12g error=%.12g pt=%f, phip=%f set=%i\n",result,error,pt,phip,thisset);
	  exit(1);
	}
    }


  /* restore original handler */
  gsl_set_error_handler (new_handler);

  /*
  F.function = &firstintegrand;
  F.params = &parameters;
  double temp;
  gsl_integration_qag(&F,0,2*M_PI,1e-10,1e-4,2000,3,w,&temp,&error);
  F.function = &secondintegrand;
  F.params = &parameters;

  gsl_integration_qag(&F,0,2*M_PI,error,1e-4,2000,3,w,&result,&error);
  */
  //printf("comp %f\n",temp-result);
  
  result*=taus[thisset];

  return result;

}

double prepareint(double T,double m0,double pt, double phip)
{
  

  double temp;
  double *tarr;
  tarr = new double[totalnum];

  for (int i=0;i<totalnum;i++)
    {
      
      tarr[i]=ointegrate1(T,m0,pt,phip,i);
      //printf("%i %.12g\n",i,tarr[i]);
    }
  
  //printf("last should be at %i %f\n",switcher-1,taus[switcher-1]);

  //gsl_spline_init (tsspline,tarrx,tarr,switcher-1);
  gsl_spline_init (tsspline,taus,tarr,totalnum-1);




  temp=gsl_spline_eval_integ (tsspline,taus[0],taus[totalnum-1] ,tsacc);

  //printf("done1: %f\n",temp);
  
  if (temp>0)
    printf("negative result %.12g at pt=%f phi=%f\n",temp,pt,phip);

  
  temp/=2*M_PI;
  temp/=2*M_PI;
  temp/=2*M_PI;
  temp*=2;
  temp*=5.06842;
  temp*=5.06842;
  temp*=-5.06842;

  

  delete [] tarr;
  
  return temp;
}

double dummyshell(double tau,void * params)
{
  double temp;
  temp=gsl_spline_eval(tsspline,tau,tsacc);
  return temp;
}

/*
double doint()
{
  double dum=1.0;
  size_t neval;

  double result,error;

  gsl_function F;
  F.function = &dummyshell;
  F.params = &dum;
  
  gsl_integration_qag(&F,taus[0],taus[switcher-1],1e-10,1e-3,2000,3,w,&result,&error);

  printf("result %f\t error %f\n",result,error);

  result/=2*M_PI;
  result/=2*M_PI;
  result/=2*M_PI;
  result*=2;
  result*=5.06842;
  result*=5.06842;
  result*=-5.06842;

  return result;
}
*/

void testing()
{

  double result1,error1;
  double result2,error2;
  size_t neval;
  double parameters[5];

  double T=0.15;
  double m0=0.13957;
  double pt=0.1;
  int thisset=1;

   parameters[0]=T;
   parameters[1]=m0;
   parameters[2]=pt;
   parameters[4]=(double) thisset;
   parameters[3]=0.1;

   //printf("why %f\n",firstintegrand(0.5,parameters));

   double mux;
   double muy;

   double mx,nx;
   double my,ny;
   int jk=0;
   double rr=0;
   double nrr=0;
   double ut=0;
   double ur=0;
   double drdt=0;

   for (double phip=0;phip<2*M_PI;phip+=0.1)
    {
      jk++;
      mx=gsl_spline_eval (xspline[thisset], phip, xacc[thisset]);
      my=gsl_spline_eval (yspline[thisset], phip, yacc[thisset]);
      nx=gsl_spline_eval (xspline[thisset+1], phip, xacc[thisset+1]);
      ny=gsl_spline_eval (yspline[thisset+1], phip, yacc[thisset+1]);
      mux=gsl_spline_eval (uxspline[thisset], phip, uxacc[thisset]);
      muy=gsl_spline_eval (uyspline[thisset], phip, uyacc[thisset]);
      //printf("rr= %f\n",sqrt(mx*mx+my*my));
      rr+=sqrt(mx*mx+my*my);
      nrr+=sqrt(nx*nx+ny*ny);
      ut+=sqrt(1+mux*mux+muy*muy);
      ur+=(mx*mux+my*muy);
      drdt+=(sqrt(nx*nx+ny*ny)-sqrt(mx*mx+my*my))/(taus[thisset+1]-taus[thisset]);
    }

   
   nrr/=jk;
   rr/=jk;
   drdt/=jk;
   ut/=jk;
   ur/=jk*rr;

   double mt=sqrt(pt*pt+m0*m0);
   //double ut=sqrt(1+mux*mux+muy*muy);
   //double rr=sqrt(mx*mx+my*my);
   //double ur=(mx*mux+my*muy)/rr;

   printf("Mean r=%f \t at t=%f\n",rr,taus[thisset]);
   printf("Mean drdt=%f ur=%f,ut=%f\n",drdt,ur,ut);
   //printf("Mean nr=%f \t at t=%f\n",nrr,taus[thisset+1]);
   /*
   
   for (double phip=0;phip<2*M_PI;phip+=0.1)
    {
      mx=gsl_spline_eval (xspline[thisset], phip, xacc[thisset]);
      my=gsl_spline_eval (yspline[thisset], phip, yacc[thisset]);
      double dpmx=gsl_spline_eval_deriv (xspline[thisset], phip, xacc[thisset]);
      double dpmy=gsl_spline_eval_deriv (yspline[thisset], phip, yacc[thisset]);
      firstintegrand(phip,parameters);
      printf("drdt %f with %f\n",rr*drdt,drdt);
      //printf("{%f,%f,%f},",phip,my,rr*sin(phip));
      //printf("{%f,%f,%f},",phip,dpmx,-rr*sin(phip));
      //printf("{%f,%f,%f}\n",phip,dpmx,-rr*sin(phip));
      //printf("{%f,%f,%f}\n",phip,dpmy,rr*cos(phip));
    }
   cout << endl;
   */

   double tt1,tt2;
   tt1=2*M_PI*drdt*rr*mt*gsl_sf_bessel_K1(mt/T*ut)*gsl_sf_bessel_I0(pt/T*ur);
   //tt1=2*M_PI*gsl_sf_bessel_K1(mt/T*ut)*gsl_sf_bessel_I0(pt/T*ur);
   tt2=2*M_PI*pt*rr*gsl_sf_bessel_K0(mt/T*ut)*gsl_sf_bessel_I1(pt/T*ur);
   //tt=2*M_PI*gsl_sf_bessel_I0(pt/T*ur);
    //tt2=gsl_sf_bessel_K0(mt/T*ut);

   
   /*
  for (double phip=0;phip<2*M_PI;phip+=0.1)
    {
     
      parameters[3]=phip;
      

      gsl_function F1,F2;
      F1.function = &firstintegrand;
      F1.params = &parameters;
  
      F2.function = &secondintegrand;
      F2.params = &parameters;

      gsl_integration_qag(&F1,0,2*M_PI,1e-10,1e-3,2000,3,w,&result1,&error1);
      gsl_integration_qag(&F2,0,2*M_PI,1e-3,1e-3,2000,3,w,&result2,&error2);
       //double temp1=integrate1(0.15,0.13957,0.1,phip,5);
      //printf("phip=%f t1= %f\t t2=%f \n",phip,temp1);
      printf("phip=%f \t res1 = %f\t res2=%f\n",phip,result1,result2);
    }
   */
  printf("comto\t res1=%f \t res2=%f\n",tt1,tt2);
   
  tt1/=2*M_PI;
  tt1/=2*M_PI;
  tt1/=2*M_PI;
  tt1*=2;
  tt1*=5.06842;
  tt1*=5.06842;
  tt1*=-5.06842;

  tt2/=2*M_PI;
  tt2/=2*M_PI;
  tt2/=2*M_PI;
  tt2*=2;
  tt2*=5.06842;
  tt2*=5.06842;
  tt2*=-5.06842;
  printf("normed\t res1=%f \t res2=%f\n",tt1,tt2);

   /*
    for (double phip=0;phip<2*M_PI;phip+=0.1)
    {
      mx=gsl_spline_eval (xspline[thisset], phip, xacc[thisset]);
      my=gsl_spline_eval (yspline[thisset], phip, yacc[thisset]);

      double r=sqrt(mx*mx+my*my);
      printf("x =%f \t vs. rcos %f\n",gsl_spline_eval (xspline[thisset], phip, xacc[thisset]),rr*cos(phip));
      //printf("r =%f \t vs. rr %f\n",r,rr);
    }
   */
}

void generatetab()
{
  double gsfact=1;
  double tempmass=0;
  double oldmass=0;
  char buffer[maxline];
  int i=0;
  double resbuff[PTASIZE][PHIPASIZE]; //length should match length of pt*phi array!
  double ptbuff[PTASIZE];
  double phipbuff[PHIPASIZE];

  switch (PHIPASIZE) {
  case 2:
    for(i=0;i<1;i++){
      phipbuff[1-i] = 0.25*M_PI*(gaus2x[i] + 1.0);
      phipbuff[i] = 0.25*M_PI*(1.0 - gaus2x[i]);
    }
    break;
  case 4:
    for(i=0;i<2;i++){
      phipbuff[3-i] = 0.25*M_PI*(gaus4x[i] + 1.0);
      phipbuff[i] = 0.25*M_PI*(1.0 - gaus4x[i]);
    }
    break;
  case 8:
    for(i=0;i<4;i++){
      phipbuff[7-i] = 0.25*M_PI*(gaus8x[i] + 1.0);
      phipbuff[i] = 0.25*M_PI*(1.0 - gaus8x[i]);
    }
    break;
  case 10:
    for(i=0;i<5;i++){
      phipbuff[9-i] = 0.25*M_PI*(gaus10x[i] + 1.0);
      phipbuff[i] = 0.25*M_PI*(1.0 - gaus10x[i]);
    }
    break;
  case 12:
    for(i=0;i<6;i++){
      phipbuff[11-i] = 0.25*M_PI*(gaus12x[i] + 1.0);
      phipbuff[i] = 0.25*M_PI*(1.0 - gaus12x[i]);
    }
    break;
  case 16:
    for(i=0;i<8;i++){
      phipbuff[15-i] = 0.25*M_PI*(gaus16x[i] + 1.0);
      phipbuff[i] = 0.25*M_PI*(1.0 - gaus16x[i]);
    }
    break;
  case 20:
    for(i=0;i<10;i++){
      phipbuff[19-i] = 0.25*M_PI*(gaus20x[i] + 1.0);
      phipbuff[i] = 0.25*M_PI*(1.0 - gaus20x[i]);
    }
    break;
  case 48:
    for(i=0;i<24;i++){
      phipbuff[47-i] = 0.25*M_PI*(gaus48x[i] + 1.0);
      phipbuff[i] = 0.25*M_PI*(1.0 - gaus48x[i]);
    }
    break;
  default:
    printf(" No abscissas for nPhi = %i !\n",PHIPASIZE);
    printf(" GOOD BYE AND HAVE A NICE DAY! \n");
    exit(0);
  }

  while (!massfile.eof())
    {
      massfile >> tempmass;
      gsfile >> gsfact;
      namesfile.getline(buffer,maxline,'\n');

      //printf("buffer0: %i\n",(int) buffer[0]);
     
      int carret=0;
      if (((int)buffer[0])!=carret)
	{
	  printf("Generating table for %s with mass %f and spin gs=%f\n",buffer,tempmass,gsfact);
	  if (tempmass!=oldmass)
	    {
	      int j=0;
	      //phip-table
	      for (double pt=0.01;pt<PTMAX;pt+=PTMAX/PTASIZE)
		{
		  for(int k=0;k<PHIPASIZE;k++)
		    {
		      //prepareint(TF,tempmass,pt,phipbuff[k]);
		      
		      ptbuff[j]=pt;
		      //printf("pt = %f phi =% f\n",pt,phipbuff[k]);
		      
		      
		      //get integral times spin degeneracy factor
		      resbuff[j][k]=gsfact*prepareint(TF,tempmass,pt,phipbuff[k]);
		 
		      /* Pade:
			 double pade=0;
			 double mocketa=ETAOS;
			 ETAOS=0.0001;
			 pade=gsfact*prepareint(TF,tempmass,pt,phipbuff[k]);
			 ETAOS=mocketa;
			 resbuff[j][k]=pade*pade/(2*pade-gsfact*prepareint(TF,tempmass,pt,phipbuff[k]));
		      */

		      
		      if (isnan(resbuff[j][k])!=0)
			  printf("Problem at %f %f\n",pt,phipbuff[k]);

		      //printf("result =%f\n",resbuff[j][k]);
		      //printf("pt = %f phi =% f res=%f\n",pt,phipbuff[k],resbuff[j][k]);
		    }
		  j++;
		}
	    }
	  else
	    {
	      //don't do anything, just repeat last result
	    }
	  
	  for (int k=0;k<PHIPASIZE;k++)
	    {
	      for (int j=0;j<PTASIZE;j++)
		{
		  //pttab << ptbuff[j] << "\t";
		  pttab << resbuff[j][k] << "\t";
		}
	      pttab << "\n";
	    }


	  oldmass=tempmass;
      
	  i++;
	}
    }

  printf("Finished main loop\n");

  //generate pt-table:
  
  ptfile << PTASIZE << endl;
  ptfile << PHIPASIZE << endl;
  for (int j=0;j<PTASIZE;j++)
    {
      ptfile << ptbuff[j] << "\n";
    }

  printf("Done!\n");
  
}



void singlept(double mass)
{
  
  double tempmass=mass;
  double oldmass=0;
  char buffer[maxline];
  int i=0;
  double resbuff[PTASIZE][PHIPASIZE]; //length should match length of pt*phi array!
  double ptbuff[PTASIZE];
  
  int j=0;
  //phip-table
  for (double pt=0.01;pt<PTMAX;pt+=PTMAX/PTASIZE)
    {
      int k=0;
      for (double phip=0;phip<M_PI;phip+=M_PI/PHIPASIZE)
	{
	  
	  prepareint(TF,tempmass,pt,PHIPASIZE);
	  
	  ptbuff[j]=pt;
	  //resbuff[j][k]=doint();
	  resbuff[j][k]=prepareint(TF,tempmass,pt,phip);
	  
	  printf("pt=%f phip=%f result =%f\n",ptbuff[j],phip,resbuff[j][k]);
	  k++;
	}
      j++;
    }

    
  for (int k=0;k<PHIPASIZE;k++)
    {
      for (int j=0;j<PTASIZE;j++)
	{
	  //pttab << ptbuff[j] << "\t";
	  pttab << resbuff[j][k] << "\t";
	}
      pttab << "\n";
    }
  
  
  //generate pt-table:
  
  ptfile << PTASIZE << endl;
  ptfile << PHIPASIZE << endl;
  for (int j=0;j<PTASIZE;j++)
    {
      ptfile << ptbuff[j] << "\n";
    }
  
  printf("Done!\n");
  
}




int main (void)
{

  extern void readParameters(const char*);

  readParameters("data/params.txt");

  //open data file

  freeze_out.open("data/freezeout.dat", ios::in);

  countsets();
  
  freeze_out.close();

  allocMem();

  cout << "total number of tau's " << totalnum << endl;

  dummy.open("data/freezeout.dat", ios::in);

  readsets();

  dummy.close();

  interpolate1();

  printf("First interpolation finished\n");
  //interpolate2();
  //printf("Second interpolation finished\n");


  massfile.open("pasim.dat", ios::in);
  namesfile.open("pasinames.dat", ios::in);
  gsfile.open("gslist.dat",ios::in);

  pttab.open("data/phipspectra.dat", ios::out);
  ptfile.open("data/ptarr.dat", ios::out);

  generatetab();

  //singlept(0.13957);

  massfile.close();
  namesfile.close();
  gsfile.close();

  //testing();

  freeMem();

  ptfile.close();
  pttab.close();

  return 0;
}
