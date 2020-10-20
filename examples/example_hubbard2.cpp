#include "all.h"
#include"exactdiagg/all.h"
#include<string>

#define NQUBIT_HUBBARD2 8
#define NQUBIT_HUBBARD2EXTSQ 8

using namespace std;

template<int L>
struct Hubbard2
{
    bool periodic=false;
    double U,t=1, tp;
    Hubbard2(double U, double tp,bool periodic) : periodic(periodic), U(U), tp(tp) {}
    int toId(int i,int s) const  { i=i%L; return s+i*2 ; }
    FermiOp Create(int i,int s) const {return FermiOp(toId(i,s),true);}
    FermiOp Destroy(int i,int s) const {return FermiOp(toId(i,s),false);}

    QOperator Ham() const
    {
        QOperator h;
        for(int i=0;i<L-1+periodic; i++)
        {
            for(int s=0;s<2;s++)
            {
                h .Add( Create(i,s)*Destroy(i+1,s),t );
                h .Add( Create(i+1,s)*Destroy(i,s),t );

                if (!periodic && i+2>=L) continue;
                h .Add( Create(i,s)*Destroy(i+2,s),tp );
                h .Add( Create(i+2,s)*Destroy(i,s),tp );
            }
            h .Add( Create(i,0)*Destroy(i,0)*
                    Create(i,1)*Destroy(i,1),U );
        }
        return h;
    }

    SymmetryGroup<2*L,cmpx> SymTraslation()
    {
        const int Lt=2*L;
        auto T2=TranslationOp<Lt>(2);
        return CyclicGroupPow<Lt>(T2, L);
    }
    SymmetryGroup<2*L,double> SymRefl()
    {
        const int Lt=2*L;
        auto Refl=TensorPow<L,2,ElementaryOp<L>> ( ReflectionOp<L> );
        return Z2_Group<Lt>(Refl);
    }

    QOperator SiSj(int i,int j) const
    {
        QOperator s;
        for(int a=0;a<2;a++)
            for(int b=0;b<2;b++)
            {
                s .Add( Create(i,a)*Destroy(i,b)*Create(j,b)*Destroy(j,a),0.5 );
                s. Add( Create(i,a)*Destroy(i,a)*Create(j,b)*Destroy(j,b),-0.25);
            }
        return s;
    }
};


void TestHubbard2Basic(double U, double tp)
{
    const int Lt=NQUBIT_HUBBARD2;
    const int L=Lt/2;
    int nPart=L;

    auto hub=Hubbard2<L>{U,tp,true};
    auto H=hub.Ham();
    auto G=hub.SymTraslation(); //or hub.SymRefl();
    auto b=FockBasisFixedChargeG<Lt>(nPart,G,0);

    auto gs=FindGS<Lt>(H,b,G,nPart);
    b.SetSym(G,1);
    for(int nu=1;nu<G.nSym();nu++,b.SetSym(G,nu))
    {
        auto gsn=FindGS<Lt>(H,b,G,nPart);
        gs=std::min( gs, gsn );
    }
    cout<<"GS: nPart="<<gs.nPart<<" sym="<<gs.sym<<" ener="<<gs.ener<<"\n";
    b.SetSym(G,gs.sym);
    //Measurements
    {
        QOperator nd=hub.Create(0,0)*hub.Destroy(0,0)*
                hub.Create(0,1)*hub.Destroy(0,1);
        auto ndv=cdot(gs.state,  nd.toMatrix(b,G,b) * gs.state);
        cout<<"double occupancy="<<real(ndv)<<endl;
    }
    {
        QOperator nd=hub.Create(0,0)*hub.Destroy(1,0);
        auto ndv=cdot(gs.state,  nd.toMatrix(b,G,b) * gs.state);
        cout<<"c0c1="<<real(ndv)<<endl;
    }
    {
        QOperator nd=hub.Create(0,0)*hub.Destroy(2,0);
        auto ndv=cdot(gs.state,  nd.toMatrix(b,G,b) * gs.state);
        cout<<"c0c2="<<real(ndv)<<endl;
    }
    for(int i=1;i<L;i++)
    {
        auto s0si=cdot(gs.state,hub.SiSj(0,i).toMatrix(b,G,b)*gs.state);
        cout<<"S0.S"<<i<<"="<<real(s0si)<<endl;
    }
    gs.Print(b);
}

//------------------------------- Hubbard con V a primeros vecinos --------------------------

template<int L>
struct Hubbard2ExtSquare
{ //no esta puesto el hopping a segundos vecinos
    bool periodic=true;
    double U, V, tp;
    double t=0.25, epsilon=-(0.5*U+2*V);
    Hubbard2ExtSquare(double U, double V, double tp, bool periodic) : periodic(periodic), U(U), V(V), tp(tp) {}
    int toId(int i,int s) const  { i=i%L; return s+i*2 ; }      //1updw_2updw_3updw_4updw
    FermiOp Create(int i,int s) const {return FermiOp(toId(i,s),true);}
    FermiOp Destroy(int i,int s) const {return FermiOp(toId(i,s),false);}

    QOperator Ham() const
    {
        QOperator h;
    //---------------------1__2---------------------------------------
    //---------------------4__3---------------------------------------
        for(int i=0;i<L;i++)
        {
            h .Add( Create(i,0)*Destroy(i,0)*
                    Create(i,1)*Destroy(i,1),U );
            for(int s=0;s<2;s++)
                h .Add( Create(i,s)*Destroy(i,s),epsilon );
        }
        for(int i=0;i<L-1+periodic; i++)
        {
            for(int s=0;s<2;s++)
            {
                h .Add( Create(i,s)*Destroy(i+1,s),t );
                h .Add( Create(i+1,s)*Destroy(i,s),t );
                for(int sp=0;sp<2;sp++)
                    h .Add( Create(i,s)*Destroy(i,s)*
                            Create(i+1,sp)*Destroy(i+1,sp),V );
                if (!periodic && i+2>=L) continue;
                h .Add( Create(i,s)*Destroy(i+2,s),tp );
                h .Add( Create(i+2,s)*Destroy(i,s),tp );
            }
        }

        return h;
    }

    SymmetryGroup<2*L,cmpx> SymTraslation()
    {
        const int Lt=2*L;
        auto T2=TranslationOp<Lt>(2);             //traslada en 2, con dimension total de qubits
        return CyclicGroupPow<Lt>(T2, L);
    }
    SymmetryGroup<2*L,double> SymRefl()
    {
        const int Lt=2*L;
        auto Refl=TensorPow<L,2,ElementaryOp<L>> ( ReflectionOp<L> );
        return Z2_Group<Lt>(Refl);
    }

    QOperator SiSj(int i,int j) const
    {
        QOperator s;
        for(int a=0;a<2;a++)
            for(int b=0;b<2;b++)
            {
                s .Add( Create(i,a)*Destroy(i,b)*Create(j,b)*Destroy(j,a),0.5 );
                s. Add( Create(i,a)*Destroy(i,a)*Create(j,b)*Destroy(j,b),-0.25);
            }
        return s;
    }
};

void TestHubbard2V(double U, double V, double tp)
{
    const int Lt=NQUBIT_HUBBARD2EXTSQ;
    const int L=Lt/2;
    int nPart=L;

    auto hub=Hubbard2ExtSquare<L>(U,V,tp,true);
    auto H=hub.Ham();
    auto G=hub.SymTraslation();

//    auto b=FockBasisFixedCharge<Lt>(nPart);
//    auto gs=FindGS<Lt>(H,b,nPart);
    EigenStateG<cmpx> gs;
    for(int nu=0;nu<G.nSym();nu++)
    {
        auto b=FockBasisFixedChargeG<Lt>(nPart,G,nu);
        auto gsn=FindGS<Lt>(H,b,G,nPart);
        gs=std::min( gs, gsn );
    }
    cout<<"GS: nPart="<<gs.nPart<<" sym="<<gs.sym<<" ener="<<gs.ener<<"\n";
    auto b=FockBasisFixedChargeG<Lt>(nPart,G,gs.sym);
    gs.Print(b);

    //Measurements
    // QOperator nd=hub.Create(0,0)*hub.Destroy(0,0)*
    //              hub.Create(0,1)*hub.Destroy(0,1);
    // auto ndv=cdot(gs.state,  nd.toMatrix(b,G,b) * gs.state);    //sin simetrías qué pongo acá en vez de G?
    // cout<<"double occupancy= "<<real(ndv)<<endl;
    // for(int i=1;i<L;i++)
    // {
    //     auto s0si=cdot(gs.state,hub.SiSj(0,i).toMatrix(b,G,b)*gs.state);
    //     cout<<"S0.S"<<i<<"= "<<real(s0si)<<endl;
    // }
}


