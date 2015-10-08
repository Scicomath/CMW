#ifndef Random_h
#define Random_h

#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

class Random
{

protected:
    int seed;
    static Random*  srand;
    //static bool     first;

public:
    Random() {seed = ::time(0); srand=0;}
    Random(int s) { seed = s;srand=0;}
    virtual ~Random() { }
    virtual void   setSeed()      {seed = ::time(0); srand48(seed);}
    virtual void   setSeed(int s) {srand48(s); seed=s;}
    virtual double rand()         { return drand48();}
    double Gauss(double mean, double sigma);

    static void setRandom(Random* r) {
	if(srand) {
	    cerr << "warning: setRandom should be called once!!!"
	       	<< endl;
	    return;
	}
	srand = r;
    }

    static Random* getRandom()   {return srand;}
    static double getRand()   {return srand->rand();}
};

//extern Random *gRandomX;

#endif
