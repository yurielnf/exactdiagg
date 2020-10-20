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

        in.ignore(1000,'\n');
    }
}
