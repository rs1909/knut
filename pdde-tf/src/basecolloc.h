#ifndef BASECOLLOC_H
#define BASECOLLOC_H

class Vector;

class BaseColloc
{
  public:
    virtual void   Init(const Vector& sol, const Vector& par) = 0;
    virtual void   Star(Vector& out, const Vector& sol) = 0;
    virtual double Integrate(const Vector& v1, const Vector& v2) = 0;
    virtual void   meshAdapt(Vector& newprofile, const Vector& profile, Vector& newtangent, const Vector& tangent) = 0;
};

#endif /*BASECOLLOC_H*/
