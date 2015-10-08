#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <convert.h>
#include <math.h>

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

int Middle;

//controls value of tau_Pi
double COEFF=3.0;
double L1COEF=2.0;
int FREEZE=1;

int xx,yy,BF;
double px,py,Yme,m;

fstream freeze_out;

//takes x,y, \pm1, p_x,p_y,Y-\xi and m as arguments
//gives back f(tau,x,y,eta
int main(int argc, char **argv)
{
  
  extern void readParameters(const char*);
  
  readParameters("data/params.txt");
 
  Middle=(NUMT-1)/2+1;

  //get parameters
  
  if (argc<8)
    {
      printf("Wrong syntax!\n");
      printf("arguments: x,y (as integer lattice sites), +/- 1, px,py,Y-xi,m\n");
      exit(1);

    }

  xx=atoi(argv[1]);
  yy=atoi(argv[2]);
  BF=atoi(argv[3]);
  px=atof(argv[4]);
  py=atof(argv[5]);
  Yme=atof(argv[6]);
  m=atof(argv[7]);

  printf("got parameters: \n");
  printf("x=%i, y=%i, ",xx,yy);
  if (BF==1)
    printf("Bose statistics, ");
  if (BF==-1)
    printf("Fermi statistics, ");
  printf("px=%f, py=%f, Y-xi=%f, m=%f\n",px,py,Yme,m);
    

  //open file

  freeze_out.open("data/freezeout.dat", ios::in);

  double dummy[8];
  double temp;

  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	if ((sx==xx)&&(sy==yy))
	    {
	      freeze_out >> dummy[0]; //(sx-Middle)/5.06842*AT
	      freeze_out >> dummy[1]; //(sy-Middle)/5.06842*AT
	      freeze_out >> dummy[2]; //u[0][sx][sy]
	      freeze_out >> dummy[3]; //u[1][sx][sy]
	      freeze_out >> dummy[4]; //pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy]))
	      freeze_out >> dummy[5]; // pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) << "\t";
	      freeze_out >> dummy[6]; // piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) << "\t";
	      freeze_out >> dummy[7]; // T(sx,sy)/AT << "\n";
	      
	      printf("extracting info at x=%i (%f),y=%i (%f): \n",sx,dummy[0],sy,dummy[1]);
	      
	    }
	else
	  {
	    
	    freeze_out >> temp; //(sx-Middle)/5.06842*AT
	    freeze_out >> temp; //(sy-Middle)/5.06842*AT
	    freeze_out >> temp; //u[0][sx][sy]
	    freeze_out >> temp; //u[1][sx][sy]
	    freeze_out >> temp; //pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy]))
	    freeze_out >> temp; // pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) << "\t";
	    freeze_out >> temp; // piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) << "\t";
	    freeze_out >> temp; // T(sx,sy) << "\n";
	  }
	  
      }

  freeze_out >> argv[1];
  freeze_out >> temp;
  printf("Proper time tau=%f\n",temp);
  
  freeze_out.close();

  printf("ux=%f, uy=%f, pixx=%f, pixy=%f, piyy=%f, T=%f\n",dummy[2],dummy[3],dummy[4],dummy[5],dummy[6],dummy[7]);

  double f0,ut,pt;

  
  ut=sqrt(1+dummy[2]*dummy[2]+dummy[3]*dummy[3]); //alias u^tau
  temp=ut;
  pt=sqrt(m*m+px*px+py*py);//alias m_T
  temp*=pt; 
  temp*=cosh(Yme);
  temp+=px*dummy[2]+py*dummy[3];

  f0=1/(exp(temp)+BF);

  //printf("f0 %.12g\n",f0);

  double vx,vy;
  vx=dummy[2]/ut;
  vy=dummy[3]/ut;

  temp=(vx*vx*dummy[4]+2*vx*vy*dummy[5]+vy*vy*dummy[6])*pt*pt*cosh(Yme)*cosh(Yme);
  temp-=(vx*dummy[4]+vy*dummy[5])*2*pt*px*cosh(Yme);
  temp-=(vx*dummy[5]+vy*dummy[6])*2*pt*py*cosh(Yme);
  temp+=px*px*dummy[4];
  temp+=2*px*py*dummy[5];
  temp+=py*py*dummy[6];

  temp-=pt*pt*sinh(Yme)*sinh(Yme)*(dummy[4]*(1-vx*vx)-dummy[5]*2*vx*vy+dummy[6]*(1-vy*vy));
  
  temp/=2*dummy[7]*dummy[7];

  printf("f(x,p)= %.12g\n",f0*(1+(1-BF*f0)*temp));

}
  
