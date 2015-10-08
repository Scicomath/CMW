#ifndef CentralityCut_h
#define CentralityCut_h

#include "Integral.h"
#include "OverLap.h"

class CentralityCut : public Integral
{
private:
    double eps;
    double crossAA;
    double sig;
    OverLap* overlap;

public:
    CentralityCut(OverLap* ov);
    double getImpactParam(const double c,double bmin=0.0, double bmax=20.0);
    double getTotalCross() const {return crossAA;}
    double getAverageImpactParameter(double b0, double b1) {
	return 2.0/3.0*(b1*b1*b1-b0*b0*b0)/(b1*b1-b0*b0);}
    double getTAA(const double bmin, const double bmax);
    double getNBC(const double bmin, const double bmax);
    double getNPart(const double bmin, const double bmax);
    double print(const double cent1,const double cent2,double b1);
private:
    double evalPar(const double b);
};
#endif

