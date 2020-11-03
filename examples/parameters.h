#ifndef PARAMETERS_H
#define PARAMETERS_H


#include<string>
#include<vector>

struct Parameters
{
    Parameters(const char filename[]);

    std::string model="KHCN";
    int length=5;
    double t1=0.5;
    double t2=0.5;
    double U=4;
    double V=2;
    double J=0;
    int nPart=10;
    double phi=0;

    bool periodic=false;
    double t=0.0;
    double tp=0.0;
    double tpp=0.0;
    double tOO=0.0;
    double U3=0.0;
    double mu=0.0;
    double muCu=0.0;
    int Sz=0;
};


#endif // PARAMETERS_H
