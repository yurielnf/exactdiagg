#include"exactdiagg/all.h"
#include<string>
#include<valarray>

#define NQUBIT_2QD 16

using namespace std;

template<int L>
struct DoubleQD
{
    double U, t=25, gamma=5, tp=0.5;

    DoubleQD(double U) : U(U) {}
    int toId(int i,int s) const  { return i+s*L ; }
    FermiOp Create(int i,int s) const {return FermiOp(toId(i,s),true);}
    FermiOp Destroy(int i,int s) const {return FermiOp(toId(i,s),false);}

    QOperator Ham() const
    {
        QOperator h;
        // kinetic energy bath
        for(int i=1;i<L-1; i++)
            for(int s=0;s<4;s++)
            {

                h .Add( Create(i,s)*Destroy(i+1,s),t );
                h .Add( Create(i+1,s)*Destroy(i,s),t );
            }
        {// hoping to impurity
            int i=0;
            for(int s=0;s<4;s++)
            {

                h .Add( Create(i,s)*Destroy(i+1,s),gamma );
                h .Add( Create(i+1,s)*Destroy(i,s),gamma );
            }
        }
        {// hoping inter-dot
            int i=0;
            h .Add( Create(i,0)*Destroy(i,2),tp );
            h .Add( Create(i,1)*Destroy(i,3),tp );
        }

        {//interaction
            int i=0;
            h .Add( Create(i,0)*Destroy(i,0)*
                    Create(i,1)*Destroy(i,1),U );
            h .Add( Create(i,2)*Destroy(i,2)*
                    Create(i,3)*Destroy(i,3),U );
        }
        return h;
    }

    SymmetryGroup<4*L,cmpx> Sym4Channel()
    {
        const int Lt=4*L;
        auto T=TranslationOp<Lt>(L);
        return CyclicGroupPow<Lt>(T, 4);
    }
};


void Test2Qd(double U)
{
    const int Lt=NQUBIT_2QD;
    const int L=Lt/4;
    int nPart=2*L;

    auto sys=DoubleQD<L>(U);
    auto H=sys.Ham();
    auto G=sys.Sym4Channel();

    auto b=FockBasisFixedChargeG<Lt>(nPart,G,0);
    auto gs=FindGS<Lt>(H,b,G,nPart);
    for(int nu=1;nu<G.nSym();nu++,b.SetSym(G,nu))
    {
        auto gsn=FindGS<Lt>(H,b,G,nPart);
        gs=std::min( gs, gsn );
    }
    cout<<"GS: nPart="<<gs.nPart<<" sym="<<gs.sym<<" ener="<<gs.ener<<"\n";
    b.SetSym(G,gs.sym);
    //Measurements
    ofstream out("QvsU.txt",std::ios_base::app);
    out<<U;
    for(int s=0;s<1;s++)
    {
        QOperator n0=sys.Create(0,s)*sys.Destroy(0,s);
        auto n0v=cdot(gs.state,  n0.toMatrix(b,G,b) * gs.state);
        out<<" "<<real(n0v);
    }
    out<<endl;
}
