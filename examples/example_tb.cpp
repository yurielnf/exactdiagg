#include "example_tb.h"
#include"exactdiagg/all.h"

#define NQUBIT_TB 8

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

void TestHamiltonianTB()
{
    const int L=NQUBIT_TB;
    int nPart=L/2;

    auto Gr=Z2_Group<L>( ReflectionOp<L> );
    auto Geh=Z2_Group<L> ( ParticleHoleOp<L> );
    auto T1=TranslationOp<L>(1);
    auto Gt=CyclicGroupPow<L>(T1, L);

    auto G=Gr.DirectProd(Geh);
//    auto G=Gt;

    auto ham=HamiltonianTB(L,true);

    auto b=FockBasisFixedChargeG<L>(nPart,G,0);
    auto gs=FindGS<L>(ham,b,G,nPart);
    for(int nu=1;nu<G.nSym();nu++,b.SetSym(G,nu))
    {
        auto gsn=FindGS<L>(ham,b,G,nPart);
        gs=std::min( gs, gsn );
    }
    cout<<"GS: nPart="<<gs.nPart<<" sym="<<gs.sym<<" ener="<<gs.ener<<endl;
    cout<<"exac="<<ExactEnergyTB(L,nPart)<<endl;
}
