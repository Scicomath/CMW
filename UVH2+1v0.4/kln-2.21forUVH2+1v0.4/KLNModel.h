#ifndef KLNModel_h
#define KLNModel_h

#include "Saturation.h"
#include "Bases.h"
#include "Integral.h"
//#include "MVGlueDist.h"
//#include "KTGlueDist.h"

//...KLN kt factorization model for saturation.

class KLNModel : public Saturation , public Bases, public Integral
{
private:
    double Qsmin2,Qsmax2;
    double Np;
    double Ptmin, Ptmax;
    int  optLocal;
    int  optTA;
    double probA,probB;
    OverLap* overlap0;
    //MVGlueDist* MVglue;
    //KTGlueDist* KTglue;
    double expR1,expR2;
public:
    KLNModel(OverLap* ov,double srt,int mode,UnintegPartonDist* f);
    virtual ~KLNModel();


    void   setOptTA(int i) {optTA=i;}
    double getdNdy(double b,double x,double y,double h);
    double getdNdy(double b,double h);
    double getdNdy(double h);
    double getdNdyInteg(double b,double h);

    double getdNdpt(double x, double y,double b,double pt, double rap);
    double getdNdptInteg(double b,double pt,double h);

    double B,X,Y,Eta,Tau;
    void setB(double n)    {B=n;}
    void setX(double n)    {X=n;}
    void setY(double n)    {Y=n;}
    void setEta(double n)  {Eta=n;}

    double getB(){return B;}
    double getX(){return X;}
    double getY(){return Y;}
    double getEta(){return Eta;}

    double SaturationScaleX(double x,double npart);
    double SaturationScaleX();
protected:
    double dndyKL();
    double dndyKLN();
    double dndyKLNmodified();
    double dndyKLmodified();
    double integral();
    double waveFunction(double qssq,double x,double kt2,double prob);
    double func(double* x);
    double approxIntegral();
    double evalPar(const double pt);
};

#endif


