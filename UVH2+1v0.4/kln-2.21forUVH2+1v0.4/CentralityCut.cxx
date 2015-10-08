
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include "CentralityCut.h"

using namespace std;

CentralityCut::CentralityCut(OverLap* ov)
{
    eps = 1e-3;
    overlap = ov;
    sig = overlap->getSigNNInelastic();

    crossAA = Gauss(0.0,20.0,eps);

}

double CentralityCut::getImpactParam(const double c,double bmin,double bmax)
{
    //double bmin = 0.0;
    //double bmax = 20.0;
    double b,cen;
    int itry=0;
    do {
	itry++;
	b = 0.5*(bmin+bmax);
	cen = Gauss(0.0,b,eps)/crossAA*100;
	if(cen < c) bmin = b;
	if(cen > c) bmax = b;
	//cout << "cen= " << cen << " bmin= " << bmin << " bmax= " << bmax
	//     << " c= " << c
	//    << endl;
	if(itry > 100) {
	    cerr << "CentralityCut::getImpactParam infinit loop? "
		 << itry
		 << endl;
	    exit(1);
	}
    } while(fabs(cen-c) > 1e-3);

    //cout << "# cen= " << cen << " itry= " << itry << endl;
    return b;
}


double CentralityCut::getNBC(const double bmin, const double bmax)
{
    double b[38],bw[38];
    OverLap::Gauss38(bmin,bmax,b,bw);
    double nbc = 0.0;
    for(int i=0;i<38;i++) {
	nbc += overlap->getNBC(b[i])*bw[i]*2*M_PI*b[i];
    }
    return nbc;
}

double CentralityCut::getTAA(const double bmin, const double bmax)
{
    double aa = M_PI*(bmax*bmax-bmin*bmin);
    return getNBC(bmin,bmax)/aa;
}

double CentralityCut::getNPart(const double bmin, const double bmax)
{
    double b[38],bw[38];
    OverLap::Gauss38(bmin,bmax,b,bw);
    double npart = 0.0;
    for(int i=0;i<38;i++) {
	npart += overlap->getNPart(b[i])*bw[i]*2*M_PI*b[i];
    }
    double aa = M_PI*(bmax*bmax-bmin*bmin);
    return npart/aa;
}

// Geometric cross section in which at least one inelastic NN collision
// occures
double CentralityCut::evalPar(const double b)
{
    return 2*M_PI*b*(1-exp(-overlap->getNBC(b)*sig));
}

double CentralityCut::print(const double cent1,const double cent2,double b1)
{
    if(b1<0.0) b1 = getImpactParam(cent1);
    double b2 = getImpactParam(cent2);
    double bave = getAverageImpactParameter(b1,b2);
    double taa = getTAA(b1,b2);
    double aa = M_PI*(b2*b2-b1*b1);

    cout << "# Centrality b1-b2  <b>  <TAA>  <Ncoll>  <Npart> Ncoll(b) Npart(b)"
	 << endl;
    cout << cent1 << "-" << cent2 << "%:"
       	 << setprecision(4) << setw(6) << b1 << " - "
	 << setprecision(4) << setw(6) << b2
         <<  setprecision(5) << setw(8) << bave
         <<  setprecision(5) << setw(8) << taa/10
         <<  setprecision(5) << setw(8) << taa*sig
         <<  setprecision(5) << setw(8) << getNPart(b1,b2)
	 <<  setprecision(5) << setw(8) << overlap->getNBC(bave)*sig
         <<  setprecision(5) << setw(8) << overlap->getNPart(bave)
         <<  setprecision(5) << setw(8) << taa*aa
	 << endl;

    return b2;
}


