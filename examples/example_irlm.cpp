#include"exactdiagg/all.h"
#include<string>
#include<valarray>

#define NQUBIT_2QD 10

using namespace std;

template<int L>
struct IRLM
{
    arma::mat tmat, Pmat;
    double U;

    IRLM(string tFile,string PFile, double U_)
        :U(U_)
    {
        tmat.load(tFile);
        Pmat.load(PFile);
    }
    FermiOp Create(int i) const {return FermiOp(i,true);}
    FermiOp Destroy(int i) const {return FermiOp(i,false);}

    QOperator Ham() const
    {
        QOperator h;
        // kinetic energy bath
        for(int i=0;i<L; i++)
            for(int j=0;j<L; j++)
        {
            if (fabs(tmat(i,j))<1e-15) continue;
            h .Add( Create(i)*Destroy(j),tmat(i,j) );
        }

        {//interaction
            for(int a=0;a<L; a++)
                for(int b=0;b<L; b++)
                    for(int c=0;c<L; c++)
                        for(int d=0;d<L; d++)
                            h .Add( Create(a)*Create(b)*
                                    Create(c)*Destroy(d),U*Pmat(0,a)*Pmat(1,b)*Pmat(1,c)*Pmat(0,d) );
        }
        return h;
    }

    SymmetryGroup<L,cmpx> Translation()
    {
        auto T=TranslationOp<L>(1);
        return CyclicGroupPow<L>(T, L);
    }
};


void TestIRLM(double U)
{
    const int Lt=NQUBIT_2QD;
    const int L=Lt;
    int nPart=L/2;

    auto sys=IRLM<L>("","",U);
    auto H=sys.Ham();
    auto G=sys.Translation();

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
    {
        QOperator n0=sys.Create(0)*sys.Destroy(0);
        auto n0v=cdot(gs.state,  n0.toMatrix(b,G,b) * gs.state);
        out<<" "<<real(n0v);
    }
    out<<endl;
}
