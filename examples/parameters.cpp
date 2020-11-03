#include "parameters.h"

#include"parameters.h"
#include<fstream>
#include<iostream>
#include<string>
#include<stdexcept>

using namespace std;

Parameters::Parameters(const char filename[])
{
    ifstream in(filename);
    if (!in.is_open())
        throw std::invalid_argument("I couldn't open parameter file");
    string param;
    while (!in.eof())
    {
        in>>param;
        if (param=="model")   in>>model;
        else if(param=="length") in>>length;
        else if(param=="t1") in>>t1;
        else if(param=="t2") in>>t2;
        else if(param=="U") in>>U;
        else if(param=="V") in>>V;
        else if(param=="J") in>>J;
        else if(param=="nPart") in>>nPart;
        else if(param=="phi") in>>phi;

        else if(param=="periodic") in>>periodic;
        else if(param=="t") in>>t;
        else if(param=="tp") in>>tp;
        else if(param=="tpp") in>>tpp;
        else if(param=="tOO") in>>tOO;
        else if(param=="U3") in>>U3;
        else if(param=="mu") in>>mu;
        else if(param=="muCu") in>>muCu;
        else if(param=="Sz") in>>Sz;

        in.ignore(1000,'\n');
    }
}
