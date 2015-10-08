#ifndef OVERLAP_h
#define OVERLAP_h

#include <iostream>
#include <cstdlib>
#include <cmath>

class OverLap
{
protected:
    double rad,rmaxCut,rwMax;
    double dr;
    double density0;
    double A;
    double impactPar;
    double Nc;
    double lambdaQCD2;
    double rmax;
    int    nevent;
    double* TA;
    double  sig;  // inelastic NN cross section [fm^2].
    double  TATAMax;

    // working area.
    double zini;
    double zfin;
    double z[38],zw[38];

public:
    // a=atomic number, b=impact paaramete [fm], sigin: [mb]
    OverLap(int a, double b, double signn, int seed=1231);
    virtual ~OverLap();
    double getSigNNInelastic() const {return sig;} // [fm^2]
    double getMassNumber() const {return A;}
    double getRadius() const {return rad;}
    double getDiffuseness() const {return dr;}
    void setRadius(double r) { rad=r;}
    void setDiffuseness(double d)  {dr=d;}

    double getNCollSq(double b);

    double getTAA(double x,double y, double b);
    double getThickness(double x, double y);
    double getNBC(double b); // [1/fm^2]
    double getNPart(double b);
    double getSaturationScale(double b); // [GeV^2]
    double getNPartDensity(double b);  //[1/fm^2]
    double getNPartDensity(double x, double y,double b);
    double getNPart1(double x, double y,double b);
    double getProb(double x, double y,double b);
    double getNPart1(double b);
    void   getRandomBC(double b, double& x, double& y);
    double getRmax(double b) {return rad+10.0*dr;}
    //double getRmax(double b) {return 10.0;}
    double getXmax(double b) {return getRmax(b)-0.5*b;}
    double getYmax(double b) {
	double rmax = getRmax(b);
	if(rmax*rmax - 0.25*b*b <=0.0) {
	    std::cerr << "OverLap:: funny rmax = " << rmax
	         << " b= " << b << std::endl;
	    exit(1);
	}
	return sqrt(rmax*rmax - 0.25*b*b);
    }
    void setTATAMax(double b) {
        TATAMax = getThickness(0.5*b,0.0)*getThickness(-0.5*b,0.0);
    }
    void getRandomWS(double& x, double& y);


    static void Gauss38(double xini,double xfin,double* xn,double* wn);


};



#endif

