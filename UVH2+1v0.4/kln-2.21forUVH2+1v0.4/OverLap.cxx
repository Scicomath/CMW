#include <cstdlib>
#include <cmath>
#include <algorithm>
//#include "pythia/Pythia6.h"
#include "OverLap.h"
//#include "Const.h"

#include "Random.h"

using namespace std;

OverLap::OverLap(int a, double b, double signnin, int seed)
{
    //srand48(seed);

    A = (double)a;
    rad=1.12*pow(A,0.333333)-0.86/pow(A,0.333333);
    dr=0.54;
    rmaxCut = rad + 2.5;
    rwMax = 1.0/(1.0+exp(-rad/0.54));

    //rad=6.38;
    //dr=0.53;         // PHENIX

    density0 = 0.17;
    rmax = 10.0;
    nevent = 3000;
    TA = 0;
    sig = signnin/10.0;

    Nc=3;
    lambdaQCD2=0.2*0.2;

    zini = -10.0;
    zfin = 10.0;
    Gauss38(zini,zfin,z,zw);
    impactPar = b;
    TATAMax = getThickness(0.5*b,0.0)*getThickness(-0.5*b,0.0);
}

OverLap::~OverLap()
{
    if(TA) delete [] TA;
}

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
double OverLap::getNBC(double b)
{
    double xx[38],xxw[38],yy[38],yyw[38];
    double glauber = 0.0;

    Gauss38(zini,zfin,xx,xxw);
    Gauss38(zini,zfin,yy,yyw);
    //Gauss38(zini,zfin,z,zw);
       
    for(int ix=0;ix<38;ix++) {
	for(int iy=0;iy<38;iy++) {
	    double t1 = 0.0;
	    double t2 = 0.0;
	    for(int iz=0;iz<38;iz++) {
		double r1 = sqrt((xx[ix]-b/2.0)*(xx[ix]-b/2.0)
			+yy[iy]*yy[iy]+z[iz]*z[iz]);
		double r2 = sqrt((xx[ix]+b/2.0)*(xx[ix]+b/2.0)+yy[iy]*yy[iy]
			+z[iz]*z[iz]);
		double ws1 = exp(-(r1-rad)/dr)/(1.0+exp(-(r1-rad)/dr));
		double ws2 = exp(-(r2-rad)/dr)/(1.0+exp(-(r2-rad)/dr));
		t1 += density0*ws1*zw[iz];
		t2 += density0*ws2*zw[iz];
	    }
	    glauber += t1*t2*xxw[ix]*yyw[iy];
	}
    }

    return glauber;

}
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
double OverLap::getTAA(double x, double y, double b)
{
    double t1=0.0;
    double t2=0.0;
    double x1 = x - 0.5*b;
    double x2 = x + 0.5*b;
    for(int iz=0;iz<38;iz++) {
	double r1 = sqrt(x1*x1 + y*y + z[iz]*z[iz]);
	double r2 = sqrt(x2*x2 + y*y + z[iz]*z[iz]);
	double ws1 = exp(-(r1-rad)/dr)/(1.0+exp(-(r1-rad)/dr));
	double ws2 = exp(-(r2-rad)/dr)/(1.0+exp(-(r2-rad)/dr));
	t1 += density0*ws1*zw[iz];
	t2 += density0*ws2*zw[iz];
    }
    return t1*t2;

}

double OverLap::getThickness(double x, double y)
{
    //double zini = -10.0;
    //double zfin = 10.0;
 
    //Gauss38(zini,zfin,z,zw);
       
    double t1 = 0.0;
    for(int iz=0;iz<38;iz++) {
        double r1 = sqrt(x*x + y*y + z[iz]*z[iz]);
        //double ws1 = exp(-(r1-rad)/dr)/(1.0+exp(-(r1-rad)/dr));
        double ws1 = 1.0/(1.0+exp((r1-rad)/dr));
        t1 += density0*ws1*zw[iz];
    }
    return t1;
}

double OverLap::getNPart1(double x, double y,double b)
{
    double t1=getThickness(x+b/2.0,y);
    double t2=getThickness(x-b/2.0,y);
    return t1*(1.0 - pow((1.0-sig*t2/A),A));
}

double OverLap::getProb(double x, double y,double b)
{
    double t1=getThickness(x+b/2.0,y);
    return (1.0 - pow((1.0-sig*t1/A),A));
}

double OverLap::getNPart1(double b)
{
    double npart = 0.0;
    double npartsq = 0.0;
    for(int ix=0;ix<38;ix++) {
	for(int iy=0;iy<38;iy++) {
	    //npart += getNPart1(z[ix],z[iy],b)*zw[ix]*zw[iy];

	    double np = getNPart1(z[ix],z[iy],b);
	    npartsq += np*np*zw[ix]*zw[iy];
	    npart += np*zw[ix]*zw[iy];
	}
    }
    return npartsq/npart;
}

double OverLap::getNPartDensity(double x, double y,double b)
{
    double t1=getThickness(x-b/2.0,y);
    double t2=getThickness(x+b/2.0,y);
    return t1*(1.0 - pow((1.0-sig*t2/A),A))+t2*(1.0-pow(1.0-sig*t1/A,A));
}

double OverLap::getNPart(double b)
{
    double npart = 0.0;
    for(int ix=0;ix<38;ix++) {
	for(int iy=0;iy<38;iy++) {
	    npart += getNPartDensity(z[ix],z[iy],b)*zw[ix]*zw[iy];
	}
    }
    return npart;
}

double OverLap::getNPartDensity(double b)
{
    double npartsq = 0.0;
    double npart = 0.0;
    for(int ix=0;ix<38;ix++) {
	for(int iy=0;iy<38;iy++) {
	    double np = getNPartDensity(z[ix],z[iy],b);
	    npartsq += np*np*zw[ix]*zw[iy];
	    npart += np*zw[ix]*zw[iy];
	}
    }
    return npartsq/npart;
}

double OverLap::getNCollSq(double b)
{
    double ncl = 0.0;
    for(int ix=0;ix<38;ix++) {
	for(int iy=0;iy<38;iy++) {
	    double ta = getThickness(z[ix]+0.5*b,z[iy]);
	    ncl += ta*ta*zw[ix]*zw[iy];
	}
    }
    return ncl/A;
}



//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//....a binary collision point (x,y) is sampled randomly according to
//....number of binary collision profile with Wood-Saxon nuclear density.
//...Inputs:
void OverLap::getRandomBC(double b, double& x, double& y)
{
    //double xmax = rmax - 0.5*b;
    //double ymax = sqrt(rmax*rmax - 0.25*b*b);
    //double lmin = -10.0;
    //double lmax =  10.0;
    double xmax =  getXmax(b);
    double xmin =  - xmax;
    double ymax =  getYmax(b);
    double ymin =  - ymax;
    do {
	x = (xmax-xmin)*Random::getRand() + xmin;
	y = (ymax-ymin)*Random::getRand() + ymin;
    } while(Random::getRand() > 
	    getThickness(x+0.5*b,y)*getThickness(x-0.5*b,y)/TATAMax);

}

void OverLap::getRandomWS(double& x, double& y)
{
    double r = 0.0;
    do {
	r = rmaxCut*pow(Random::getRand(),1.0/3.0);
    } while(Random::getRand()*rwMax > 1.0/(1.0+exp((r-rad)/0.54)));

    //r *= radius/rmax;

    double cx = 1.0-2.0*Random::getRand();
    double sx = sqrt(1.0-cx*cx);
    double phi = 2*M_PI*Random::getRand();
    x = r*sx*cos(phi);
    y = r*sx*sin(phi);

}

//**********************************************************************
void OverLap::Gauss38(double xini,double xfin,double* xn,double* wn)
{
    double  x[38], w[38];

    x[37]=9.980499305357e-1;
    x[36]=9.897394542664e-1;
    x[35]=9.748463285902e-1;
    x[34]=9.534663309335e-1;
    x[33]=9.257413320486e-1;
    x[32]=8.918557390046e-1;
    x[31]=8.520350219324e-1;
    x[30]=8.065441676053e-1;
    x[29]=7.556859037540e-1;
    x[28]=6.997986803792e-1;
    x[27]=6.392544158297e-1;
    x[26]=5.744560210478e-1;
    x[25]=5.058347179279e-1;
    x[24]=4.338471694324e-1;
    x[23]=3.589724404794e-1;
    x[22]=2.817088097902e-1;
    x[21]=2.025704538921e-1;
    x[20]=1.220840253379e-1;
    x[19]=4.078514790458e-2;

//    .....   WEIGHT       ...........
    w[37]=5.002880749632e-3;
    w[36]=1.161344471647e-2;
    w[35]=1.815657770961e-2;
    w[34]=2.457973973823e-2;
    w[33]=3.083950054518e-2;
    w[32]=3.689408159400e-2;
    w[31]=4.270315850467e-2;
    w[30]=4.822806186076e-2;
    w[29]=5.343201991033e-2;
    w[28]=5.828039914700e-2;
    w[27]=6.274093339213e-2;
    w[26]=6.678393797914e-2;
    w[25]=7.038250706690e-2;
    w[24]=7.351269258474e-2;
    w[23]=7.615366354845e-2;
    w[22]=7.828784465821e-2;
    w[21]=7.990103324353e-2;
    w[20]=8.098249377060e-2;
    w[19]=8.152502928039e-2;

    for(int i=0;i<19;i++) {
	x[i] = -x[37-i];
	w[i] =  w[37-i];
    }
    for(int i=0;i<38;i++) {
         xn[i] =(xfin-xini)*x[i]/2.0+(xini+xfin)/2.0;
         wn[i] =(xfin-xini)*w[i]/2.0;
    }

}



//#define MAIN 1

#ifdef MAIN
#include <iomanip>
int main() {

    double sig=41.0; // 130GeV.
    //double sig=42.0; // 200GeV.
    OverLap* overlap = new OverLap(197,0.0,sig);

    double bmax=14.0;
    double db=1.0;
    int nn= (int)(14.0/db);
    for(int i=0;i<nn;i++) {
	double b = i*db;
	cout << setw(10) << b << setw(10) << overlap->getNBC(b)*sig/10
	     << setw(10) << overlap->getNPart(b)
	     <<	setw(10) << overlap->getNPartDensity(b)
	     << setw(10) << overlap->getSaturationScale(b)
	     << endl;
    }

    return 0;
}
#endif
