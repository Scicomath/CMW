#ifndef KLNfunc_h
#define KLNfunc_h
#include "Saturation.h"
#include "UnintegPartonDist.h"

class KLNfunc: public UnintegPartonDist
{
private:
    double fac;
    double lam;
public:
    KLNfunc() {
	fac = 1.0/Saturation::FacSQ/M_PI;
	//lam=0.2*0.2;
	lam=0.0;
    }
    double getFunc(double qs2, double x, double kt2, double alp) {
	if(kt2 <= qs2) return fac/alp*qs2/(qs2+lam);
	return fac*qs2/(kt2+lam)/alp;
    
	//return fac*qs2/(kt2+qs2+lam)/alp;
    }

};
#endif

