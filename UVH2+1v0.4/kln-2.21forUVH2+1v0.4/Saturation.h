#ifndef Saturation_h
#define Saturation_h

#include "UnintegPartonDist.h"
#include "OverLap.h"
//#include "IntegralMulti.h"
#include "Integral.h"
//#include "book/Book1.h"
#include <cmath>


//...kt factorization model for gluon saturation in AA collisions.

class Saturation
{
protected:
    OverLap* overlap;
    UnintegPartonDist* wavefunc;
    double ecm;
    double impactPar;
    int    model;
    int    dNdeta;       // =1: calculate dN/deta
    int    dEtdy;        // =1: calculate dE/dy
    int    optFixAlpha;
    double alphaS;      // alpha strong for fix case.
    double lambda;
    double lambdaQCD;
    double lambdaQCD2;
    double srt0;
    double ccc;
    double QsCut;
    double KGlue;
    double mHadron;


    double Qs2;
    double NParticipants;
    double NPartDensity;
    double participantsArea;
    double rapidity;
    double pT;
    double transEtaY;
    double NPart1, NPart2;

    //Book1 *histQ2, *histX;

    double tau0;
    double theta;
    double gammaS;

public:
    static double Nc;
    static double Nf;
    static double CF;
    static double FacSQ;
    static double Norm;
    static double Beta0;

public:
    Saturation(OverLap* ov,double srt,int mode, UnintegPartonDist* f);
    virtual ~Saturation();
    double getGamma(double x, double q2, double qs2) {
	double y = log(1.0/x);
	double lgq = std::abs(log(q2/qs2));
	return gammaS + (1.0-gammaS)*lgq/(lambda*y + 1.3*sqrt(y) + lgq);
    }

    void   loadHist();
    void   printHist();

    void setTau0(double t){tau0=t;}
    void setAlphaS(double s) {alphaS=s;}
    void setOptFixAlpha(int i) {optFixAlpha=i;}

    void setTheta(double t){theta=t;}
    double getTheta(){return theta;}

    virtual double getdNdy(double b,double x,double y,double h)=0;
    virtual double getdNdy(double b,double h)=0;
    virtual double getdNdy(double h)=0;
    virtual double getdNdyInteg(double b,double h)=0;

    virtual double getdNdpt(double x, double y,double b,double pt, double rap)
	=0;
    virtual double getdNdptInteg(double b,double pt,double h)=0;

    double getNpart1() {return NPart1;}
    double getNpart2() {return NPart2;}
    double getSaturationScale1();
    double getSaturationScale();
    double SaturationScale0(double npart);
    double xG(double q2);
    double getJacobian(double h);

    inline double getAlphaStrong(const double q2);


    double getQs2() const {return Qs2;}
    double getImpactPar() const {return impactPar;}

    void setNParticipants(double n)      {NParticipants=n;}
    void setNPartDensity(double n)       {NPartDensity=n;}
    void setSaturationScale(double q)    {Qs2=q;}
    void setImpactPar(double b)          {impactPar=b;}
    void setKGlue(double k)              {KGlue=k;}

    void setModel(int i)                 {model = i;}
    void setdNdeta(int i)                {dNdeta=i;}
    void setdEtdy(int i)                 {dEtdy=i;}
    void setLambdaQCD(double l)          {lambdaQCD=l; lambdaQCD2=l*l;}
    void setNoramlization(double c)      {ccc=c;}
    void setLambdaSaturation(double l)   {lambda = l;}
    void setMHadron(double m)            {mHadron=m;}

    int    getModel()            const {return model;}
    double getLambdaQCD()        const {return lambdaQCD;}
    double getNoramlization()    const {return ccc;}
    double getLambdaSaturation() const {return lambda;}
    double getKGlue()            const {return KGlue;}
    int    getdEtdy()            const {return dEtdy;}
    double getMHadron()          const {return mHadron;}


    // test functions.
    void NpartNcollQs();
    double getAverageQs(double b);

};

//inline double max(double a, double b) { return a>=b ? a:b;}
inline double Min(double a, double b) { return a<=b ? a:b;}

inline double Saturation::getAlphaStrong(const double q2)
{

    if(optFixAlpha) return alphaS;

    if(q2 <= lambdaQCD2) return 0.5;
    return Min(0.5,  1.0/( Beta0 * log( q2/lambdaQCD2 ) ) );

    //if(q2 <= lambdaQCD2) return 0.25;
    //return Min(0.25,  1.0/( Beta0 * log( q2/lambdaQCD2 ) ) );

    //return 1.0/( Beta0 * log( (q2 + 0.2*0.2) /lambdaQCD2 ) );

}

#endif


