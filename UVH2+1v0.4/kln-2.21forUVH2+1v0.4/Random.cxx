#include "Random.h"

#include <cmath>

Random* Random::srand = 0;
//bool EventBase::first = true;

double Random::Gauss(double mean, double sigma)
{
    double  x, y, z, result;
    do {
	y = rand();
    } while (!y);
    z = rand();
    x = z * 6.283185;
    result = mean + sigma*sin(x)*sqrt(-2*log(y));
   return result;

}

