#include "FTDipole.h" 
#include <iostream>
#include <iomanip>

//Parameters:       
//cme ........ center of mass energy [GeV], 
//aeff ....... effective mass number; e.g. 18.5 gives the saturation scale \
//Q_s(x=0.01)~1 GeV at mid-rapidity,     
// lamb ..... 0.3 (GB-W fit to HERA data),     
//  x0 ...... 3.0*10^(-4) ( GB-W fit to HERA data),      
//  gammas .. 0.627 (Leading-order BFKL anomalous dimension with saturation
//                   boundary),     
//  d ....... a parameter in the DHJ anomalous dimension = 1.2
//              (DHJ fit d+Au RHIC data),     
// color factor .... color factor C_F/C_A = (2/3)^2 for N_F
//                                    and 1 for N_A at N_c=3     
//  Variables:     
//  yh ......... rapidity of the observed hadron,     
//  xp ......... momentum fraction carried by the projectile parton,     
//  ya[xp] ..... rapidity of the target-side gluon(s),     
//  qs[xp] ..... saturation scale for the target,     
//  qt[xp] ..... transverse momentum of the scattered parton,     
//  xfmin ...... Feynman-x (=p_t/cme*exp[yh]) at p_t=1 [GeV] (lowest-p_t),
//               the lower bound for the integral over xp,     
//  xpmax ...... the upper bound for the integral over xp,     
//  xpnum ...... total number of steps for xp,     
//  xpstep ..... step size in xp [= (xpmax-xfmin)/xnum]

using namespace std;

FTDipole::FTDipole()
{
    N=1000;                    // number of gauss points
    //mxp=10000;
    mxp=N;
    xg0 = new double [mxp];
    wg0 = new double [mxp];
    xg1 = new double [mxp];
    wg1 = new double [mxp];

    lamb=0.3;
    x0=3e-4;
    gammas=0.627;
    d=1.2;
    Aeff=18.5;
    yh=4.0;

//.....make gaus-values between standard -1,1 boundaries
//.....and transform to some other range
//
//     one can use the xg0/wg0  or the xg1/wg1 gauss points
//     see gausrange for what they mean

    gauleg(-1.0,1.0,xg0,wg0,N); // get gauss legendre points
    gausrange(1,0.0,100.0,xg0,wg0,N); // 1 is A to B

    gauleg(-1.0,1.0,xg1,wg1,N); // get gauss legendre points
    gausrange(2,1.0,0.0,xg1,wg1,N); // 2 is 0 to inf

}
FTDipole::~FTDipole()
{
    delete [] xg0;
    delete [] wg0;
    delete [] xg1;
    delete [] wg1;
}

double FTDipole::getFunc(double qs2, double xp, double kt2, double alp)
{
    double sum=0.0;
    for(int i=0;i<N;i++) {
        double rt=xg1[i];

        //double gam1=exp( -0.25*pow( qs2/kt2*rt*rt,gamma(xp,rt,qs2,kt2) ) );
	//double BJ0 = BesselJ0(rt);

        double gam1=exp( -0.25*pow( qs2*rt*rt,gamma(xp,rt,qs2,1.0) ) );
        //double gam1=exp( -0.25*qs2*rt*rt );
	double BJ0 = BesselJ0(sqrt(kt2)*rt);

        sum += wg1[i]*rt*BJ0*gam1;
    }
    //if(sum <0.0) return 0.0;
    if(sum <0.0) {
	//cout << " Negative Qs2= " << qs2 << " kt2= " << kt2
	 //   << " x= " << xp << endl;
	return 0.0;
    }
    return kt2/(2*M_PI*alp)*sum;
}

void FTDipole::output()
{
    cme=200.0;
    double xpmax=1; // 0.03 // for yh=0    
    int xpnum=25;
    double xpmin=1.0/cme*exp(yh);
    double xpstep=(xpmax-xpmin)/xpnum;

//.....Fourier transform with the Bessel function
      
      //do xp=xpmin,xpmax+xpstep/2,xpstep
    for(int ix=0;ix<=xpnum;ix++) {
	double xp = xpmin + ix*xpstep;
        double sum1=0.0;
        double sum2=0.0;
        double f=1.0/(2*M_PI)/(qt(xp)*qt(xp));
	for(int i=0;i<N;i++) {
          double rt=xg1[i];
          double gam1=f*exp( -0.25*  pow( cff(xp)*rt,2*gamrrt(xp,rt) ) );
          double gam2=f*exp( -0.25*  pow( caa(xp)*rt,2*gamrrt(xp,rt) ) );
	  double BJ0 = BesselJ0(rt);
          sum1 += wg1[i]*rt*BJ0*gam1;
          sum2 += wg1[i]*rt*BJ0*gam2;
	}
	/*
	cout.precision(9.8);
	cout << setw(15) << xp
	   << setw(15) << sum1
	   << setw(15) << sum2
	   << setw(15) << qt(xp)
	   << endl;
	   */

	cout << setprecision(8)  << xp
	   << setprecision(8) << setw(14) << sum1
	   << setprecision(8) << setw(14) << sum2
	   << setprecision(8) << setw(14) << qt(xp)
	   << endl;
    }

}

void FTDipole::gauleg(double x1,double x2,double *x,double *w,int n)
{
    static double EPS=3.e-14;

    int m=(n+1)/2;
    double  xm=0.5*(x2+x1);
    double xl=0.5*(x2-x1);

    for(int i=0;i<m;i++) {
    //do i=1,m
        double z=cos(M_PI*(i+1-0.25)/(n+0.5));
	double z1=z;
	double pp;
	do {
	    double p1=1.0;
	    double p2=0.0;
	    for(int j=1;j<=n;j++) {
		double p3=p2;
		p2=p1;
		p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
	    }
	    pp=n*(z*p1-p2)/(z*z-1.0);
	    z1=z;
	    z=z1-p1/pp;
	}while(abs(z-z1) > EPS);

        x[i]=xm-xl*z;
        //x[n+1-i]=xm+xl*z;
        x[n-1-i]=xm+xl*z;
        w[i]=2.0*xl/((1.0-z*z)*pp*pp);
        //w[n+1-i]=w[i];
        w[n-1-i]=w[i];

    }

}

void FTDipole::gausrange(int iop,double A,double B,double *xg1,double *wg1,int n)
{
//     transform gausspoints to other range than [-1;1]
//     iop = 1  [A,B]       uniform
//     iop = 2  [0,inf]     A is midpoint
//     opt = 3  [-inf,inf]  scale is A 
//     opt = 4  [B,inf]     A+2B is midoint
//     opt = 5  [0,B]     AB/(A+B)+ is midoint
       
    double xp,wp;
    for(int i=0;i<N;i++) {
//      do i=1,N
        if(iop==1) {
//...... A to B 
          xp=(B+A)/2+(B-A)/2*xg1[i];
          wp=(B-A)/2*wg1[i];
	} else if(iop==2) {
//...... zero to infinity
          xp=A*(1+xg1[i])/(1-xg1[i]);
          wp=2.*A/(1-xg1[i])/(1-xg1[i])*wg1[i];
	} else if(iop==3) {
//...... -inf to inf
          xp=A*(xg1[i])/(1-xg1[i]*xg1[i]);
          double tmp=(1-xg1[i]*xg1[i]);
          wp=A*(1+xg1[i]*xg1[i])/(tmp*tmp)*wg1[i];
	} else if(iop==4) {
//......  B to inf,  A+2B is midoint
          xp=(A+2*B+A*xg1[i])/(1-xg1[i]);
          wp=2.*(B+A)/(1-xg1[i])/(1-xg1[i])*wg1[i];
	} else if(iop==5) {
//...... 0 to B , AB/(A+B) is midpoint
          xp=A*B*(1+xg1[i])/(B+A-(B-A)*xg1[i]);
          double tmp = B+A-(B-A)*xg1[i];
          wp=2*A*B*B/(tmp*tmp)*wg1[i];
	} else if(iop==3) {
	    cerr << " invalid option iop= " << iop << endl;
	    exit(1);
	}
        xg1[i]=xp;
        wg1[i]=wp;
    }

}

/*
        subroutine MJY01A
C
C       =========================================================
C       Purpose: This program computes the Bessel functions  
C                Jn(x) and Yn(x) ( n=0,1 ) and their derivatives 
C                using subroutine JY01A
C       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ò 0 )
C       Output:  BJ0 --- J0(x)
C                DJ0 --- J0'(x)
C                BJ1 --- J1(x)
C                DJ1 --- J1'(x)
C                BY0 --- Y0(x)
C                DY0 --- Y0'(x)
C                BY1 --- Y1(x)
C                DY1 --- Y1'(x)
C       Example:
C
C        x       J0(x)        J0'(x)       J1(x)        J1'(x)
C       ---------------------------------------------------------
C        1     .76519769   -.44005059    .44005059    .32514710
C        5    -.17759677    .32757914   -.32757914   -.11208094
C       10    -.24593576   -.04347275    .04347275   -.25028304
C       20     .16702466   -.06683312    .06683312    .16368301
C       30    -.08636798    .11875106   -.11875106   -.08240961
C       40     .00736689   -.12603832    .12603832    .00421593
C       50     .05581233    .09751183   -.09751183    .05776256
C
C        x       Y0(x)        Y0'(x)       Y1(x)        Y1'(x)
C      ---------------------------------------------------------
C        1     .08825696    .78121282   -.78121282    .86946979
C        5    -.30851763   -.14786314    .14786314   -.33809025
C       10     .05567117   -.24901542    .24901542    .03076962
C       20     .06264060    .16551161   -.16551161    .07091618
C       30    -.11729573   -.08442557    .08442557   -.12010992
C       40     .12593642    .00579351   -.00579351    .12608125
C       50    -.09806500    .05679567   -.05679567   -.09692908
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x'
        READ(*,*)X
        WRITE(*,20)X
        WRITE(*,*)'  x          J0(x)          J0''(x)         J1(x)',
     &            '          J1''(x)'
        WRITE(*,*)'------------------------------------------',
     &            '----------------------------'
        CALL JY01A(X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
        WRITE(*,10)X,BJ0,DJ0,BJ1,DJ1
        WRITE(*,*)
        WRITE(*,*)'  x          Y0(x)          Y0''(x)         Y1(x)',
     &            '          Y1''(x)'
        WRITE(*,*)'------------------------------------------',
     &            '----------------------------'
        WRITE(*,10)X,BY0,DY0,BY1,DY1
10      FORMAT(1X,F5.1,4E16.8)
20      FORMAT(3X,'x =',F5.1)
        END
*/


double FTDipole::BesselJ0(double x)
{
//       =======================================================
//       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
//                Y1(x), and their derivatives
//       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ò 0 )
//       Output:  BJ0 --- J0(x)
//                DJ0 --- J0'(x)
//                BJ1 --- J1(x)
//                DJ1 --- J1'(x)
//                BY0 --- Y0(x)
//                DY0 --- Y0'(x)
//                BY1 --- Y1(x)
//                DY1 --- Y1'(x)
//       =======================================================

    double PI=3.141592653589793;
    double RP2=0.63661977236758;

    double x2=x*x;
    if (x==0.0) return 1.0;
    if (x <= 12.0) {
	double BJ0=1.0;
        double r=1.0;
	for(int k=1;k<=30;k++) {
              r = -0.25*r*x2/(k*k);
              BJ0 += r;
              if (abs(r) < abs(BJ0)*1.0e-15) return BJ0;
	}
	return BJ0;

    } else {
	static double A[12]={
                  -.7031250000000000e-01,.1121520996093750e+00,
                  -.5725014209747314e+00,.6074042001273483e+01,
                  -.1100171402692467e+03,.3038090510922384e+04,
                  -.1188384262567832e+06,.6252951493434797e+07,
                  -.4259392165047669e+09,.3646840080706556e+11,
                  -.3833534661393944e+13,.4854014686852901e+15};
        static double B[12]={
	           .7324218750000000e-01,-.2271080017089844e+00,
                   .1727727502584457e+01,-.2438052969955606e+02,
                   .5513358961220206e+03,-.1825775547429318e+05,
                   .8328593040162893e+06,-.5006958953198893e+08,
                   .3836255180230433e+10,-.3649010818849833e+12,
                   .4218971570284096e+14,-.5827244631566907e+16};

	int k0=12;
        if (x >= 35.0) k0=10;
        if (x >= 50.0) k0=8;
        double t1=x-0.25*PI;
        double p0=1.0;
        double q0=-0.125/x;
	for(int k=1;k<=k0;k++) {
	    p0 += A[k-1]*pow(x,(-2*k));
            q0 += B[k-1]*pow(x,(-2*k-1));
	}
        double cu=sqrt(RP2/x);
        return cu*(p0*cos(t1)-q0*sin(t1));
    }

}

//#define MAIN

#ifdef MAIN
int main() {
    FTDipole* dipole = new FTDipole();
    dipole->output();
    /*
    int pmax=1000;
    double qs2=1.0;
    for(int i=0;i<pmax;i++) {
	double p = 0.001 + i*0.1;
        cout << p/sqrt(qs2)
	    << " " << dipole->getFunc(qs2, x, p*p, 0.5 ) << endl;
    }
    */
    return 0;
}
#endif



