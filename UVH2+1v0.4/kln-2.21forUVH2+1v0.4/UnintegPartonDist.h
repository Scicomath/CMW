#ifndef UnintegPartonDist_h
#define UnintegPartonDist_h

class UnintegPartonDist
{
public:
    UnintegPartonDist() {}
    virtual ~UnintegPartonDist() {}
    virtual double getFunc(double qs2, double x, double pt2,double alp)=0;
};
#endif
