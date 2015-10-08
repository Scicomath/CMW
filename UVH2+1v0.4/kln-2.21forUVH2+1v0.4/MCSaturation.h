#ifndef MCSATURATION_h
#define MCSATURATION_h
#include <vector>
#include "Particle.h"
#include "OverLap.h"

class MCSaturation
{
protected:
    OverLap* overlap;
    std::vector<Particle*> nucl1,nucl2;
    double **TA1,**TA2;
    double **Qssq1,**Qssq2;
    int Maxx,Maxy;
    double dx,dy;
    double sig,sigel,sigin;

public:
    MCSaturation(double ecm,int maxx,int maxy,double dx0,double dy0);
    ~MCSaturation();

    double getQs1(double x,double y) {return Qssq1[(int)(x/dx)][(int)(y/dy)];}
    double getQs2(double x,double y) {return Qssq2[(int)(x/dx)][(int)(y/dy)];}
    double getTA1(double x,double y) {return TA1[(int)(x/dx)][(int)(y/dy)];}
    double getTA2(double x,double y) {return TA2[(int)(x/dx)][(int)(y/dy)];}
    double getNpart1(double x,double y,double b){
	return overlap->getNPart1(x,y,b);}
    double getNpart2(double x,double y,double b){
	return overlap->getNPart1(x,y,-b);}
    double getTA1(double x,double y,double b) {
	return overlap->getThickness(x+b/2.0,y);}
    double getTA2(double x,double y,double b) {
	return overlap->getThickness(x-b/2.0,y);}

    void generateNucleus(double b);
    void deleteNucleus();
    void getTA();
    void getQs();
    void generate(double impactpar);
};

#endif

