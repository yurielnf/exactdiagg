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
};


#endif // PARAMETERS_H
