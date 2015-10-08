#include "Saturation.h"
//#include "physicsbase/Const.h"
    const double hbarC = 0.197327053;
    const double hbarCsq = hbarC*hbarC;
#include <iostream>
#include <iomanip>
#include <algorithm>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef ALPHA
#define abs  fabs
#endif

using namespace std;

double Saturation::Nc = 3.0;
double Saturation::Nf = 3.0;
double Saturation::CF = (Nc*Nc-1.0)/(2*Nc);
double Saturation::FacSQ=2*M_PI*M_PI/CF;
double Saturation::Norm=CF/(2*M_PI*M_PI)/hbarCsq;
double Saturation::Beta0 = (33.0-2.0*Nf)/(12*M_PI);

Saturation::Saturation(OverLap* ov,double srt, int mode,UnintegPartonDist* f)
{
    overlap = ov;
    ecm = srt;
    model = mode;
    wavefunc = f;

    gammaS = 0.627;
    optFixAlpha = 0;
    alphaS = 0.5;

    dNdeta = 0;
    dEtdy = 0;
    impactPar = 2.0;
    QsCut = 0.25;
    //QsCut = 0.55;

    // Parameters for saturation model
    lambda=0.288;

    lambdaQCD=0.2;      // [GeV]
    lambdaQCD2=lambdaQCD*lambdaQCD;
    ccc = 1.0;

    KGlue=1.0;
    //mHadron = 0.14;
    mHadron = 0.4;

    //srt0 = 200.0;      //[GeV]
    srt0 = 130.0;      //[GeV]

    //histQ2=0;
    //histX=0;

}

Saturation::~Saturation()
{
    //if(histQ2) delete histQ2;
    //if(histX) delete histX;
}


void Saturation::loadHist()
{
    //histQ2 = new Book1("hist q2",50,0.0,8.0);
    //histX = new Book1("hist x",100,0.0,0.5);
}

void Saturation::printHist()
{
    //histQ2->print("q2.hist",0);
    //histX->print("x.hist",0);
}


double Saturation::getJacobian(double h)
{
    if(!dNdeta) return 1.0;

    double Qs = sqrt(Qs2);
    double meff = mHadron;
    double msq=2*Qs*meff;
    double pz=Qs*sinh(h);
    double e=sqrt(msq+Qs2+ pz*pz);

    rapidity=0.5*log((e+pz)/(e-pz));

    return Qs*cosh(h)/e;
}

double Saturation::xG(double q2)
{
    //static const double K=1.8;
    //static const double K=1.5;  // constant alpha in Qs

    //static const double K=1.2; //run 
    //static const double K=2.5; //run 

    //static const double K=2.6; //run 
    //static const double K=3.0; //run 

    // A.H.Mueller, Nucl. Phys. B643(2002)501,[hep-ph/0206216].
    // R. Baier, A. Kovner and U. A. Wiedemann, hep-ph/0305265.
    //static const double alpha = 0.5;
    //return K*alpha*(Nc*Nc-1)/(2*M_PI)*log(q2/lambdaQCD2+1.0/(9*lambdaQCD2));
    //return 3*getAlphaStrong(q2)*(Nc*Nc-1)/(Nc*2*M_PI)*log(q2/lambdaQCD2);

    //return KGlue*log(q2/lambdaQCD2 + 1.0/(9*lambdaQCD2));
    //static const double lam = 0.2*0.2;
    static const double lam = 0.0*0.0;
    return KGlue*log((q2+ lam )/lambdaQCD2);
    //return KGlue*log(q2/lambdaQCD2);

}

// Saturation scale from simple gluon density.
// before calling this, you need to set "NPartDensity"
double Saturation::getSaturationScale() // [GeV^2]
{
    double fac = FacSQ*NPartDensity*hbarCsq/2.0;
    double q0 = 2.0;
    double q2 = q0;
    int i=0;
    do {
	q0 = q2;
	if(q2 <= lambdaQCD2) return lambdaQCD2;

	//double alpha=0.5;
	//double alpha = min(0.5,getAlphaStrong(q2));
	double alpha = getAlphaStrong(q2);
	double xg;
	xg = xG(q2);
	q2 = fac*alpha*xg;
	if(i++ > 1000) {
	    cerr << "KLNModel::getSaturationScale: not converge "
		 << i
		 << " q2= " << q2  << endl;
	    return q2;
	}
    } while(abs(q0-q2) > 1e-3);

    //return q2*pow(ecm/srt0,lambda);
    return q2;
}

//...Local Q_s^2:saturation scale at a transverse coordinate (x,y)
double Saturation::getSaturationScale1()
{
    static const double beta0=(11.0-2.0*Nf/3)/4.0/M_PI;
    static const double fac = 2.0*M_PI*M_PI*Nc/(Nc*Nc-1.0)/beta0;

    double qs2 = fac*hbarC*hbarC*NPartDensity;
    return qs2*pow(ecm/srt0,lambda);

}


// calculate average local saturation scale in transverse plane.
double Saturation::getAverageQs(double b)
{

    /*
    double rmax=14.0;
    double xmax = rmax - 0.5*b;
    double ymax = sqrt(rmax*rmax - 0.25*b*b);
    double x[38],xw[38],y[38],yw[38];
    OverLap::Gauss38(-xmax,xmax,x,xw);
    OverLap::Gauss38(-ymax,ymax,y,yw);
    */

    double z[38],zw[38];
    double zini=-10.0;
    double zfin=10.0;
    OverLap::Gauss38(zini,zfin,z,zw);
       
    double q2 = 0.0;
    int ds=0;
    for(int ix=0;ix<38;ix++) {
	for(int iy=0;iy<38;iy++) {
	    NPartDensity = overlap->getNPartDensity(z[ix],z[iy],b);
	    //NPartDensity = overlap->getNPartDensity(x[ix],y[iy],b);
	    double q=getSaturationScale();
	    if(q>0.0) {
	    q2 += q*zw[ix]*zw[iy];
	    ds++;
	    }
	}
    }

    if(ds > 0) return q2/ds;
    return 0.0;

}

// test impact parameter dependence of saturation scale Qs
void Saturation::NpartNcollQs()
{
    double bmax=12.0;
    double db=1.0;
    int nn= (int)(bmax/db);
    cout << "# b (fm)   Ncoll    Npart      rho       Qs^2(GeV^2)" << endl;

    double sig = overlap->getSigNNInelastic();
    for(int i=0;i<nn;i++) {
	double b = i*db;
	double npart = overlap->getNPart(b);
	double npartdens =  overlap->getNPartDensity(b);
	setNParticipants(npart);
	setNPartDensity(npartdens);

	cout << setw(6) << b
	     << setw(10) << overlap->getNBC(b)*sig
	     << setw(10) << npart
	     <<	setw(10) << npartdens
	     << setw(10) << getSaturationScale()
	     << endl;
    }


}

double Saturation::SaturationScale0(double npart) // [GeV^2]
{
    //double fac = FacSQ*NPartDensity*hbarCsq/2.0;
    double fac = FacSQ*npart*hbarCsq;
    double q0 = 2.0;
    double q2 = q0;
    int i=0;
    double lam2 = 0.2*0.2;
    do {
	q0 = q2;
	double alpha = getAlphaStrong(q2);
	q2 = fac*alpha*KGlue*log(( q2 + lam2)/lambdaQCD2);
	if(i++ > 1000) {
	    cerr << "KLNModel::getSaturationScale: not converge "
		 << i
		 << " q2= " << q2  << endl;
	    return q2;
	}
    } while(abs(q0-q2) > 1e-3);
    return q2;
}
