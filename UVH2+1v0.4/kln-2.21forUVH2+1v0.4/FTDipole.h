#ifndef FTDipole_h
#define FTDipole_h
#include <cmath>
#include "UnintegPartonDist.h"

//     does Arata's Mathematica script much faster
//     good old fortran !! 
//     I just translated directly 
//     hajo
//     I just translated to c++ directly from hajo's   Yasushi
  

class FTDipole : public UnintegPartonDist
{
private:
      int mxp;
      double *xg0, *wg0;
      double *xg1,*wg1;
      //double cxtot[61][3],cxqel[61][3],cxel[61][3],cxbel[61][3],cxinel[61][3];
      double lamb,x0,gammas,d;
      double cme;
      double Aeff;
      double yh;
      int    N;
      
protected:
      double ya(double xp_) {return log(1.0/xp_)+2.0*yh;}
      double Qs2(double xp_) {
	  return pow(Aeff,0.3333)*exp(2*lamb*yh)*pow(x0/xp_,lamb);
      }
      double Qs(double xp_) {return sqrt(Qs2(xp_));}
      double qt(double xp_) {return cme*exp(-yh)*xp_;}

      double cff(double xp_) {return 2.0/3.0*Qs(xp_)/qt(xp_);}
      double caa(double xp_) {return Qs(xp_)/qt(xp_);}

      double gamrrt(double xp_,double rt_) {
          double lgq = std::abs( log( qt(xp_)*qt(xp_)/(rt_*rt_*Qs2(xp_))  ) );
	  return gammas+(1.0-gammas)*lgq
	      /( lamb*ya(xp_) + d*sqrt(ya(xp_)) + lgq );
      }
      double gamma(double x,double rt,double qs2, double qt2) {
          //double lgq = std::abs( log( qt2/(rt*rt*qs2))   );
          double lgq = std::abs( log( (qt2+0.04)/(rt*rt*qs2))   );
	  double y = log(1.0/x);
	  return gammas+(1.0-gammas)*lgq/( lamb*y + d*sqrt(y) + lgq );
      }

public:
    FTDipole();
    ~FTDipole();
    void output();
    void gauleg(double x1,double x2,double *x,double *w,int n);
    void gausrange(int iop,double A,double B,double *xg1,double *wg1,int n);
    double BesselJ0(double x);
    double getFunc(double qs2, double x, double kt2, double alp);

};
#endif
