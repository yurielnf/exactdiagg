#include <iostream>
#include"exactdiagg/all.h"

#define NQUBIT 8

using namespace std;

QOperator HamiltonianTB(int L,bool periodic)
{
    QOperator h;
    for(int i=0;i<L-1+periodic;i++)
    {
        h.Add(Create(i)*Destroy((i+1)%L),-1);
        h.Add(Create((i+1)%L)*Destroy(i),-1);
    }
    return h;
}

double ExactEnergyTB(int L, int nPart)
{
    std::vector<double> evals(L);
    for(int k=0;k<L;k++)
        evals[k]=2*cos(2*M_PI*k/L);
    std::sort(evals.begin(),evals.end());
    double sum=0;
    for(int k=0;k<nPart;k++)
        sum+=evals[k];
    return sum;
}

void TestTBHamiltonian()
{
    const int L=NQUBIT;
    int nPart=L/2;

    auto Gr=Z2_Group<L>( ReflectionOp<L> );
    auto Geh=Z2_Group<L> ( ParticleHoleOp<L> );
    auto T1=TranslationOp<L>(1);
    auto Gt=CyclicGroupPow<L>(T1, L);
    auto G=Gr;

    auto ham=HamiltonianTB(L,true);

    EigenStateG<double> gs;
    auto b=FockBasisFixedChargeG<L>(nPart,G,0);
    for(int nu=0;nu<G.nSym();nu++,b.SetSym(G,nu))
    {
        auto gsn=FindGS<L>(ham,b,G,nPart);
        gs=std::min( gs, gsn );
    }
    cout<<"GS: nPart="<<gs.nPart<<" sym="<<gs.sym<<" ener="<<gs.ener<<endl;
    cout<<"exac="<<ExactEnergyTB(L,nPart)<<endl;
}


int main()
{
    cout << "Hello World!" << endl;
    TestTBHamiltonian();
    return 0;
}
