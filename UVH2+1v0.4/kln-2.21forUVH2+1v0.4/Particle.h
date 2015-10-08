#ifndef PARTICLE_h
#define PARTICLE_h

#include <iostream>
#include <cstdlib>
#include <cmath>

class Particle
{
protected:
    double x,y,z;
    double vx,vy,vz;
    double absorptionProb;
    int freezeOut;
    double phi;

public:
    Particle(double x0,double y0, double z0) {
	x = x0; y = y0; z = z0;
	vx = 0.0; vy = 0.0; vz = 0.0;
	absorptionProb = 0.0;
	freezeOut=0;
	phi=0.0;

    }
    Particle(double x0,double y0, double vx0, double vy0,double ph) {
	 x = x0;   y = y0;   z = 0.0;
	vx = vx0; vy = vy0; vz = 0.0;
	absorptionProb = 0.0;
	freezeOut=0;
	phi=ph;

    }
    ~Particle() {};

    double getX() {return x;}
    double getY() {return y;}
    double getZ() {return z;}
    double getVx() {return vx;}
    double getVy() {return vy;}
    double getVz() {return vz;}
    double updateX(double dt) {return x += vx*dt;}
    double updateY(double dt) {return y += vy*dt;}
    double updateZ(double dt) {return z += vz*dt;}
    void   updateAbsProb(double a) {absorptionProb += a;}
    double getAbsProb() {return absorptionProb;}
    int    isFreezeOut() {return freezeOut;}
    void   setFreezeOut(int i) {freezeOut=i;}
    double getPhi() {return phi;}
};

#endif

