#ifndef KLNGamfunc_h
#define KLNGamfunc_h
#include "Saturation.h"
#include "UnintegPartonDist.h"
#include <cmath>

class KLNGamfunc: public UnintegPartonDist
{
private:
    double fac;
    double lam;
    double gammaS;
    double lambda;
public:
    KLNGamfunc(double la) {
	fac = 1.0/Saturation::FacSQ/M_PI;
	lam=0.2*0.2;
	gammaS = 0.627;
	lambda = la;
    }
    double getFunc(double qs2, double x, double kt2, double alp) {
	if(kt2 <= qs2) return fac/alp;
	double gam= getGamma(x,kt2, qs2);
	return fac/alp*pow(qs2/kt2,gam);
    }
    double getGamma(double x, double q2, double qs2) {
	double y = log(1.0/x);
	double lgq = std::abs(log(q2/qs2));
	return gammaS + (1.0-gammaS)*lgq/(lambda*y + 1.3*sqrt(y) + lgq);
    }

};
#endif

