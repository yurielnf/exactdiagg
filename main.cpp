#include<iostream>
#include<iomanip>
#include"examples/example_tb.h"
#include"examples/torou.h"

using namespace std;

int main(int argc,const char *argv[])
{
    cout << "Hello World!" << endl << setprecision(15);
//    TestHamiltonianTB();

    if (argc==3 && string(argv[1])=="basic")  // ./a.out basic U
    {
        double U=atof(argv[2]);
        TestHubbardBasic(U);
    }
    else if (argc==6 && string(argv[1])=="quench") // ./a.out quench U1 U2 dt nt
    {
        double U1=atof(argv[2]);
        double U2=atof(argv[3]);
        double dt=atof(argv[4]);
        int nt=atoi(argv[5]);
        TestHubbardQuench(U1,U2,dt,nt);
    }
    return 0;
}
