#include "MCSaturation.h"
#include "Regge96.h"

using namespace std;

MCSaturation::MCSaturation(double ecm,
	int maxx,int maxy,double dx0,double dy0)
{
    Maxx=maxx;
    Maxy=maxy;
    dx=dx0;
    dy=dy0;

    //....Estimate NN cross sections.
    sig = hadronxsec::totalXsection(ecm,0);
    sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
    sigin = sig - sigel;
    overlap = new OverLap(197,0.0,sigin);

    Qssq1 = new double* [Maxx];
    Qssq2 = new double* [Maxx];
    TA1 = new double* [Maxx];
    TA2 = new double* [Maxx];
    for(int ix=0;ix<Maxx;ix++) {
	Qssq1[ix] = new double[Maxy];
	Qssq2[ix] = new double[Maxy];
	TA1[ix] = new double[Maxy];
	TA2[ix] = new double[Maxy];
    }
}

MCSaturation::~MCSaturation()
{
    delete overlap;

    for(int ix=0;ix<Maxx;ix++) {
	delete [] Qssq1[ix];
	delete [] Qssq2[ix];
	delete [] TA1[ix];
	delete [] TA2[ix];
    }
    delete [] Qssq1;
    delete [] Qssq2;
    delete [] TA1;
    delete [] TA2;
}

void MCSaturation::generateNucleus(double b)
{
    for(int ia=0;ia<197;ia++) {
	double x,y;
	overlap->getRandomWS(x,y);
	x += b/2.0;
	//int ix= (int)floor(x+0.5);
	//int iy= (int)floor(y+0.5);

	nucl1.push_back(new Particle(x,y,0.0));

	overlap->getRandomWS(x,y);
	x -= b/2.0;
	nucl2.push_back(new Particle(x,y,0.0));
    }
}


void MCSaturation::getTA()
{
    double dsq = 0.1*sigin/M_PI;
    for(int ix=0;ix<Maxx;ix++) 
    for(int iy=0;iy<Maxy;iy++) {
	double x0=ix*dx;
	double y0=iy*dy;

	TA1[ix][iy]=0.0;
	for(int i=0;i<(int)nucl1.size();i++) {
	    double x = nucl1[i]->getX();
	    double y = nucl1[i]->getY();
	    double dc = (x-x0)*(x-x0) + (y-y0)*(y-y0);
	    if(dc<=dsq) TA1[ix][iy] += 10.0/sigin;
	}
	TA2[ix][iy]=0.0;
	for(int i=0;i<(int)nucl2.size();i++) {
	    double x = nucl2[i]->getX();
	    double y = nucl2[i]->getY();
	    double dc = (x-x0)*(x-x0) + (y-y0)*(y-y0);
	    if(dc<=dsq) TA2[ix][iy] +=  10.0/sigin;
	}

    }
}

void MCSaturation::getQs()
{
    for(int ix=0;ix<Maxx;ix++) 
    for(int iy=0;iy<Maxy;iy++) {
	if(TA1[ix][iy]*TA2[ix][iy]>0) {
	    Qssq1[ix][iy]=TA1[ix][iy];
	    Qssq2[ix][iy]=TA2[ix][iy];
	} else {
	    Qssq1[ix][iy]=0.0;
	    Qssq2[ix][iy]=0.0;
	}
    }
}

void MCSaturation::generate(double impactpar) 
{
    generateNucleus(impactpar);
    getTA();
    getQs();
    deleteNucleus();
}

void MCSaturation::deleteNucleus()
{
    for(int i=0;i<(int)nucl2.size();i++) {
	delete nucl1[i];
    }
    for(int i=0;i<(int)nucl2.size();i++) {
	delete nucl2[i];
    }
    nucl1.clear();
    nucl2.clear();
}

