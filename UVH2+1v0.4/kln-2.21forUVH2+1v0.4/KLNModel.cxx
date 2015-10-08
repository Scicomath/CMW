#include "KLNModel.h"
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
inline double max(double a, double b) { return a>=b ? a:b;}
inline double min(double a, double b) { return a<=b ? a:b;}
#endif

using namespace std;

KLNModel::KLNModel(OverLap* ov,double srt, int mode,UnintegPartonDist* f) 
      : Saturation(ov,srt,mode,f)
{
    optLocal = 1;
    optTA = 2;

    overlap0 = new OverLap(197,0.0,39.5277);

    //NPartDensity = overlap0->getNPartDensity(0.0,0.0,0.0);
    //Qs0 = getSaturationScale();

    //MVglue = new MVGlueDist();
    //KTglue = new KTGlueDist();
}


KLNModel::~KLNModel()
{
    //delete MVglue;
    //delete KTglue;
    delete overlap0;
}

// Local dN/dy use local Qs(x,y).
double KLNModel::getdNdy(double b,double x,double y,double h)
{
    //NPartDensity = overlap0->getNPartDensity(x,y,b);
    //NPartDensity = overlap0->getNPartDensity(0.0,0.0,0.0);
    //Qs0 = getSaturationScale();

    Qs2 = 0.0;
    NPartDensity = overlap->getNPartDensity(x,y,b);
    NParticipants = NPartDensity;
    optLocal = 1;
    double x0 = x;
    double y0 = y;
    probA = 1.0;
    probB = 1.0;

    // include free streaming, assuming massless particle (v=1)
    if(model==10){
	theta = getTheta();
	x0 = x-cos(theta)*tau0;
	y0 = y-sin(theta)*tau0;
    }

    if(optTA==1) {
	NPart1 = overlap->getThickness(x0+b/2.0,y0);
	NPart2 = overlap->getThickness(x0-b/2.0,y0);
    } else if(optTA==2) {
	probA = overlap->getProb(x0,y0,b);
	probB = overlap->getProb(x0,y0,-b);
//	if(probA < 1e-6 || probB < 1e-6) return 0.0;
	NPart1 = overlap->getThickness(x0+b/2.0,y0)/probA;
	NPart2 = overlap->getThickness(x0-b/2.0,y0)/probB;
	//cout << probA << " " <<  NPart1 << endl;
	//cin.get();
    } else {
	NPart1 = overlap->getNPart1(x0,y0,b);
	NPart2 = overlap->getNPart1(x0,y0,-b);
    }


	// pA
	//NPart1 = overlap->getThickness(x0-b/2.0,y0);
	//double qsp2= 0.6*0.6;
	//NPart2 = qsp2*1.53/KGlue;

    return getdNdy(h);
}

// dN/dy use average Qs(b) as original paper.
double KLNModel::getdNdy(double b, double h)
{
    NPartDensity = overlap->getNPartDensity(b);
    NParticipants = overlap->getNPart(b);
    optLocal = 0;
    NPart1 = overlap->getNPart1(b);
    NPart2 = overlap->getNPart1(-b);
    //Qs0 = SaturationScaleX(0.01,NParticipants);
    return getdNdy(h);
}

double KLNModel::getdNdy(double h)
{
    double p_density=0.0;
    rapidity = h;

    //Qs2 = getSaturationScale();

    transEtaY=1.0;
    double edep=1.0;
    Qsmin2=0.0;
    Qsmax2=0.0;

    if(model ==2 || model == 5) {
	Qs2 = getSaturationScale();
	transEtaY=Saturation::getJacobian(h);
	edep = pow(ecm/srt0,lambda);
    } else if(model ==3 || model == 4) {
	//Qs2 = SaturationScaleX();
	Qs2 = getSaturationScale();
	edep = pow(ecm/srt0,lambda);
	//edep = pow(sqrt(Qs2)/ecm,lambda);
	 //edep = pow(0.01*ecm/sqrt(Qs2),lambda);
	transEtaY=Saturation::getJacobian(h);
    } else if(model == 6) {
	Qs2 = SaturationScaleX();
	transEtaY=Saturation::getJacobian(h);
	double hh = rapidity;
	rapidity = abs(hh);
	Qsmin2 = SaturationScaleX();
	rapidity = -abs(hh);
	Qsmax2 = SaturationScaleX();
	rapidity = hh;
	if(Qsmin2 < 1e-5) return 0.0;
	if(Qsmax2 < 1e-5) return 0.0;
	//if(Qsmin2 <= lambdaQCD2) return 0.0;
	//if(Qsmax2 <= lambdaQCD2) return 0.0;
	if(Qsmin2 > Qsmax2) {
	    cout << "qsmin2= " << Qsmin2
		<< " qsmax2= " << Qsmax2
		<< endl;
	    return 0.0;
	}
    }

    Qs2 *= edep;

    //if(model <=5) Qs2 *= pow(ecm/srt0,lambda);

	//if(optLocal)  Qs2 *= pow(Qs0/Qs2,lambda/2);

    if(model <= 5) {

	double qs1 = SaturationScale0(NPart1)*exp(-lambda*abs(rapidity));
	double qs2 = SaturationScale0(NPart2)*exp(lambda*abs(rapidity));

	Qsmin2 = min(qs1,qs2);
	Qsmax2 = max(qs1,qs2);
	Qs2 = (Qsmin2+Qsmax2)/2;
	//cout << "qsmin= " << Qsmin2 << " qsmax= " << Qsmax2 << endl;

    //Qsmin2= Qs2*exp(-lambda*abs(rapidity));
    //Qsmax2= Qs2*exp(lambda*abs(rapidity));

    if(Qsmin2 <= lambdaQCD2) return 0.0;

    }

    participantsArea =NParticipants/NPartDensity;
    double norm = 1.0;
    if(optLocal) {
       	norm = participantsArea*Norm;
    }

    Np = Qsmin2;
    if(optLocal==0)
	Np = NParticipants*exp(-lambda*abs(rapidity))*edep;

    double et = 1.0;
    if(dEtdy && model <=6) et = sqrt(Qs2);

    switch (model) {
	case 2:
	    p_density = dndyKL();
	    break;
	case 3:
	    p_density= dndyKLN();
	    break;
	case 4:
	    p_density= dndyKLNmodified();
	    break;
	case 5:
	    p_density=dndyKLmodified();
	    break;
	case 6:
	    p_density= dndyKLN();
	    break;
        case 7: case 8: case 10:
	    p_density=integral();
	    if(optLocal==0) p_density *= participantsArea;
	    break;
 	case 9:
	    p_density=approxIntegral();
	    if(optLocal==0) p_density *= participantsArea;
	    break;
	default:
	    cerr << "(KLNModel::getdNdy:)";
	    cerr << "Invalid model= ";
	    cerr << model  << endl;
	    exit(1);
    }

    //    return ccc*norm*p_density*transEtaY*et;
    return ccc*norm*FacSQ*FacSQ*p_density*transEtaY*et;//from ver.2.332
}

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//D.~Kharzeev and M.~Nardi,
//``Hadron production in nuclear collisions at RHIC and high density QCD,''
// Phys.\ Lett.\ B {\bf 507}, 121 (2001) [arXiv:nucl-th/0012025].
//
//D.~Kharzeev and E.~Levin,
//``Manifestations of high density QCD in the first RHIC data,''
//Phys.\ Lett.\ B {\bf 523}, 79 (2001) [arXiv:nucl-th/0108006].
//
// D.~Kharzeev, E.~Levin and M.~Nardi,
//``The onset of classical QCD dynamics in relativistic heavy ion  collisions,''
//    arXiv:hep-ph/0111315.
//

double KLNModel::dndyKL()
{
    double xq = sqrt(Qsmax2)/ecm*exp(abs(rapidity));
    if ((xq > 1.0) ||  Qsmin2 <= lambdaQCD2) return 0.0;

    return Np*log(Qsmin2/lambdaQCD2)*(1+0.5*log(Qsmax2/Qsmin2)*pow(1-xq,4));

	/*
          return  pow(ecm/srt0,lambda)
	            * NParticipants*exp(-lambda*abs(rapidity))
                    * (log(Qs2/lambdaQCD2)-lambda*abs(rapidity))
                    * (1.0+lambda*abs(rapidity)*pow(1.0-xq,4));
	*/
}
    
// Ref. D.~Kharzeev, E.~Levin and M.~Nardi,
//    ``QCD saturation and deuteron nucleus collisions,''
//    arXiv:hep-ph/0212316.
double KLNModel::dndyKLN()
{
    double xqmax = sqrt(Qsmax2)/ecm*exp(abs(rapidity));
    double xqmin = sqrt(Qsmin2)/ecm*exp(abs(rapidity));
    //double xqmin = sqrt(Qsmin2)/ecm*exp(-abs(rapidity));

    if ((xqmax > 1.0) || (xqmin > 1.0)) return 0.0;

    /*
    double dndy= Np/getAlphaStrong(Qsmin2)
           *( pow(1-xqmin,4)*pow(1-xqmax,4)
		   + (log(Qsmax2/Qsmin2)+1)*pow(1-xqmax,4)*pow(1-xqmin,4) );
		   */
    //double dndy= Np/getAlphaStrong(Qsmin2)
     //   *( pow(1-xqmin,4) + 0.5*(log(Qsmax2/Qsmin2))*pow(1-xqmax,4) );


    double dndy= Np/getAlphaStrong(Qsmin2)
       *( pow(1-xqmin,4) +  0.5*log(Qsmax2/Qsmin2)*pow(1-xqmax,4) );

    // original.
    //double dndy= Np/getAlphaStrong(Qsmin2)
     //  *( pow(1-xqmin,4) + (log(Qsmax2/Qsmin2)+1)*pow(1-xqmax,4) );




    return dndy;

}

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
double KLNModel::dndyKLNmodified()
{
    double xq = sqrt(Qs2)/ecm*exp(abs(rapidity));
    if(xq > 1.0) return 0.0;

    return Np/getAlphaStrong(Qs2)
            *( pow(1-xq,4) + (log(Qsmax2/Qsmin2)+1)*pow(1-xq,4) );

}

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//...Jamal
double KLNModel::dndyKLmodified()
{
    double xq = sqrt(Qs2)/ecm*exp(abs(rapidity));
    if(xq > 1.0) return 0.0;
    return Np*log(Qs2/lambdaQCD2)*(1.0+lambda*abs(rapidity)*pow(1.0-xq,4));

}

double KLNModel::integral()
{
    Ptmin = 0.1;
    double npart = max(NPart1,NPart2);

    int opt_ptcut = 2;
    if(opt_ptcut == 1) {
    double xmin = Ptmin/ecm*exp(-lambda*abs(rapidity));
    double qq = SaturationScaleX(xmin,npart);
    if( qq <= lambdaQCD2) return 0.0;
    Ptmax = 3.0;

    } else {

    Ptmax = 30.0;
    //double npart = min(NPart1,NPart2);

    double qq=0.0;
    while(qq < 0.04) {
	Ptmax -= 0.1;
	double xmin = Ptmax/ecm*exp(abs(rapidity));
	qq = SaturationScaleX(xmin,npart);
	//if(qq >= lambdaQCD2) break;
	if(Ptmax <= Ptmin) return 0.0;
    };
    if(Ptmax <= Ptmin) return 0.0;
    if(Ptmax/sqrt(qq) > 10) Ptmax = 10*sqrt(qq);

    }


    Qs2 = SaturationScaleX(0.01,npart);
    transEtaY=Saturation::getJacobian(rapidity);

    expR1 = 1.0/ecm*exp(-rapidity);
    expR2 = 1.0/ecm*exp(rapidity);

    if(Ptmax <= Ptmin) return 0.0;

    bsinit();
    //const int ndim=2;
    int ndim=3;   // include phi dependence.
    int sample = 1000;

    if(model==10){
    ndim=3;   // include phi dependence with free streaming
    //    sample = 5000;
    //sample = 1000;
    sample = 1000;
    }
    double x_l[ndim],x_u[ndim];
    int jg[ndim];

    for(int i=0;i<ndim;i++) {
	x_l[i]=0.0;
	x_u[i]=1.0;
	jg[i]=1;
    }

    int nwild=ndim;
    //if(abs(rapidity)>4) sample = 5000;

    double tune = 1.5;
    int itr1=5;
    int itr2=5;
    double ac1=1.0;
    double ac2=1.0;
    setNoOfSample(sample);
    setTuneValue(tune);
    setIteration1(ac1,itr1);
    setIteration2(ac2,itr2);
    defineVariable(ndim,nwild,x_l,x_u,jg);

    // If you want to see the integration step status, set =2
    setPrint(0);
    double aa = bases();
    //double aa = getEstimate();
    //double err = getError();
    return aa;
}


//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double KLNModel::func(double* x)
{
    double pt= Ptmin + x[0]*(Ptmax-Ptmin);
    double ktmax = pt;
    //double ktmax = 2.0*pt;
    double kt= ktmax*x[1];

    // ignore phi dependence.
    //double ktsq2=(pt-kt)*(pt-kt);

    // include phi integral.
    double phi = 2*M_PI*x[2];
    double ktsq1 = 0.25*(pt*pt + kt*kt + 2*kt*pt*cos(phi));
    double ktsq2 = 0.25*(pt*pt + kt*kt - 2*pt*kt*cos(phi));

    double x1=pt*expR1;
    double x2=pt*expR2;
    if(x1 > 1.0 || x2 > 1.0) return 0.0;

    double qs2a= SaturationScaleX(x1,NPart1);
    double qs2b= SaturationScaleX(x2,NPart2);

    //if(ktsq1 > qs2a && ktsq2 > qs2b) return 0;
    //if(qs2a <= lambdaQCD2 && qs2b <= lambdaQCD2) return 0.0;
    /*
    if(qs2a <= lambdaQCD2 || qs2b <= lambdaQCD2) {
	cout << " pt= " << pt
	      << " kt= " << kt 
	      << endl;
	cout << " qs2a= " << qs2a
	     << " x1= " << x1
	     << " npart1= " << NPart1
	     << endl;
	cout << " qs2b= " << qs2b
	     << " x2= " << x2
	     << " npart2= " << NPart2
	     << endl;
	cin.get();
	return 0.0;
    }
    */

    double fnc1  = waveFunction(qs2a,x1,ktsq1,probA);
    double fnc2  = waveFunction(qs2b,x2,ktsq2,probB);
    double fnc = fnc1*fnc2;

    double scale = max(ktsq1,ktsq2);
    scale = max(scale,pt*pt);
    double alp = getAlphaStrong(scale);

    //    double result = 0.5*alp*kt*fnc*(Ptmax-Ptmin)*ktmax/pt;
    double result = alp*fnc; //alpha_s*phi_A*phi_B
    result *= 2.0*M_PI*(Ptmax-Ptmin)/pt; // pT integration
    result *= 2.0*M_PI*kt*ktmax;// kT integration
    result /= 4.0; //Jacobian due to symmetrization of integration
    if(dEtdy){
      return result*pt;
    }
    return result;

}

double KLNModel::SaturationScaleX(double x,double npart) // [GeV^2]
{
    //double fac = FacSQ*NPartDensity*hbarCsq/2.0*pow(x,-lambda);
    //double fac = FacSQ*npart*hbarCsq*pow(x,-lambda);

    double fac = FacSQ*npart*hbarCsq*pow(0.01/x,lambda);
    if(model==8 || model == 9 || model==10) fac *= pow(1-x,4);

    //double fac = FacSQ*NPartDensity*hbarCsq/2.0*pow(x,-lambda)*pow(1-x,4);

    return KGlue*npart/1.53*pow(0.01/x,lambda);
    //return KGlue*npart/1.53*pow(0.01/x,lambda)*pow(1-x,4)/0.9606;

    double q0 = 2.0;
    double q2 = q0;
    int i=0;
    static double lam=0.2*0.2;
    do {
	q0 = q2;
	double alpha = getAlphaStrong(q2);

	q2 = fac*alpha*KGlue*log(( q2 + lam )/lambdaQCD2);

	/*
	double ll = 0.5;
        if(q2 > lambdaQCD2) 
	    ll = max(0.5,  1.0/( Beta0 * log( q2/lambdaQCD2 ) ) );
	q2 = fac*alpha*KGlue*ll;
	*/

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

double KLNModel::SaturationScaleX() // [GeV^2]
{
    double fac = FacSQ*NPartDensity*hbarCsq/2.0;
    double q0 = 2.0;
    double q2 = q0;
    int i=0;
    double lam = 0.2;
    do {
	q0 = q2;
	if(q2 <= lambdaQCD2) return lambdaQCD2;
	double alpha = getAlphaStrong(q2);
	double x = sqrt(q2)/ecm*exp(rapidity);
	q2 = fac*alpha*KGlue*log(( q2 + lam*lam)/lambdaQCD2)*pow(x,-lambda);
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


//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//...simplest wave function in the saturation model
double KLNModel::waveFunction(double qssq,double x,double kt2,double prob)
{
    double rapdep=1.0;
    if(model==7) rapdep = pow(1.0-x,4);
    double alp = getAlphaStrong(qssq);
    //double alp = getAlphaStrong(kt2);

    return wavefunc->getFunc(qssq,x,kt2,alp)*rapdep*prob;

    // BFKL
    //double tmp = 0.5* (log(kt2) - log(qssq) );
    //return 0.03*exp( -tmp*tmp/(4*0.869) )/alp;
    //return 0.03*exp( -tmp*tmp/(4*0.169) )/alp;


    //    double fac = 1.0;
    //double fac = 1.0/FacSQ/M_PI;//from version 3.331
    //static double lam=0.2*0.2;

    // KL
    //if(kt2 <= qssq) return fac/alp*rapdep*qssq/(qssq+lam);
    //return fac*qssq/(kt2+lam)/alp*rapdep;

    //if(kt2 <= qssq) return fac/alp*rapdep;
    //double gam= getGamma(x,kt2, qssq);
    //return fac/alp*pow(qssq/kt2,gam)*rapdep;

    // GW
    //return kt2/qssq*exp(-kt2/qssq)*rapdep/alp;

    //...Itakura
    //return pow(qssq/kt2,0.64)/alp*rapdep;

    //return 1.0/alp*pow(qssq/(qssq+4*kt2),0.64)*rapdep;

    //return MVglue->getWS(sqrt(kt2),qssq);

    //return (1.0-exp(-kt2*qssq/4.0))/(kt2+0.001)*rapdep;

    //return KTglue->get(kt2+lam,qssq,Qs0,alp);

}

double KLNModel::getdNdyInteg(double b,double h)
{
    double z[38],zw[38];
    double zini=-10.0;
    double zfin=10.0;
    OverLap::Gauss38(zini,zfin,z,zw);
       
    double dNdy = 0.0;

    for(int ix=0;ix<38;ix++) {
	for(int iy=0;iy<38;iy++) {
	  double dndy = getdNdy(b,z[ix],z[iy],h);
	    dNdy += dndy*zw[ix]*zw[iy];
	}
    }

    /*
    for(int ix=0;ix<80;ix++) {
	double x = -12+ix*0.3;
	for(int iy=0;iy<80;iy++) {
	    double y = -12+iy*0.3;
	    double dndy = getdNdy(b,x,y,h);
	    dNdy += dndy*0.3*0.3;
	}
    }
    */

    return dNdy;

}

//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double KLNModel::evalPar(const double pt)
{
    double x1=pt/ecm*exp(-rapidity);
    double x2=pt/ecm*exp(rapidity);
    if(x1 > 1.0 || x2 > 1.0) return 0.0;
    double q2 = pt*pt;
    double alp = getAlphaStrong(q2);

    //double x0 = sqrt(Qs2)/130.0;
    //double qs2a= Qs2*pow(x0/x1,lambda);
    //double qs2b= Qs2*pow(x0/x2,lambda);
    //if(qs2a < lambdaQCD2 || qs2b < lambdaQCD2) return 0.0;

    double xg1,xg2;

    //double lam = 0.2*0.2;
    //xg1 = KGlue*log(( q2 + lam )/lambdaQCD2)*pow(x1,-lambda)*pow(1-x1,4);
    //xg2 = KGlue*log(( q2 + lam )/lambdaQCD2)*pow(x2,-lambda)*pow(1-x2,4);

    double qs2a= SaturationScaleX(x1,NPart1);
    double qs2b= SaturationScaleX(x2,NPart2);
    double alp1 = getAlphaStrong(qs2a);
    double alp2 = getAlphaStrong(qs2b);

    //if(kt2 <= qssq) return fac/alp*rapdep*qssq/(qssq+lam);
    //return fac*qssq/(kt2+lam)/alp*rapdep;

    if(q2 < qs2a) {
	//xg1 = q2/alp1*pow(1-x1,4);
	xg1 = q2/alp1;
    } else {
	//xg1 = qs2a/alp1*pow(1-x1,4);
	xg1 = qs2a/alp1;
    }
    if(q2 < qs2b) {
	//xg2 = q2/alp2*pow(1-x2,4);
	xg2 = q2/alp2;
    } else {
	//xg2 = qs2b/alp2*pow(1-x2,4);
	xg2 = qs2b/alp2;
    }

    //if(q2 > 2*qs2a && q2 > 2*qs2b) return 0.0;
    //if(x1 > 0.01 || x2 > 0.01) return 0.0;

    double result = alp*xg1*xg2/(pt*pt*pt);


    if(dEtdy) return result*pt;
    return result;

}

//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double KLNModel::approxIntegral()
{
    //double pmin=0.2;
    //double pmax=0.2*ecm;

    double pmin=0.0;
    double pmax=3;

    //double xg = 0.1;
    //double pmax= xg*ecm*exp(-abs(rapidity));
    //double pmax = 10;

    if(pmax <= pmin) return 0.0;
    static const double eps=1e-3;

    //NPartDensity = overlap->getNPartDensity(x,y,b);
    //double xg= pt/ecm*exp(rapidity);
    //Qs2 = getSaturationScale(xg);
    //if(Qs2 <= QsCut) return 0.0;

    double aa =  Gauss(pmin,pmax,eps);
    return aa;

}

double KLNModel::getdNdpt(double x, double y,double b,double pt, double rap)
{ return 0.0;}

double KLNModel::getdNdptInteg(double b,double pt,double h)
{return 0.0; }

