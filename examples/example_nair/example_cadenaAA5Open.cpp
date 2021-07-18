#include "examples/all.h"
#include"exactdiagg/all.h"
#include<string>
#include<complex>
#include<math.h>
#include"cadenitaaa5open.h"
#include"examples/parameters.h"

using namespace std;

//template<int Lt>
//struct HasSz
//{
//    int sz;
//    function<int(int,int,int)> toInt;
//    HasSz(int sz, function<int(int,int)> toInt):sz(sz),toInt(toInt){}
//    bool operator()(const FockState<Lt>& f)
//    {
//        int sum=0;
//        const int nOrb=Lt/2;
//        for(int i=0;i<nOrb;i++)
//            sum-=f.test(toInt(i,0));
//        for(int i=0;i<nOrb;i++)
//            sum+=f.test(toInt(i,1));
//        return sum==sz;
//    }
//};


void TestGS_CadenitaAA5(const Parameters& par)//,int nTwist,int id, int id_last)
{
    ofstream out(string("energ.txt"));
    const int Lt=4*10;  //  Lt = nq*L
    const int nq=10; //number of qubits per unit cell


//    auto T=TranslationOp<Lt>(nq);  //este traslada en la celda completa. no cambiar.
//    auto Gt=CyclicGroupPow<Lt>(T, Lt/nq);
//    auto G=Gt;

//    const int L=Lt/nq;
//    auto Refl=TensorPow<nq,L,ElementaryOp<nq>> ( ReflectionOp<nq> );
//    auto G = Z2_Group<Lt,cmpx>(Refl);



//    auto T1=TranslationOp<Lt/nq>(1);
//    auto T=TensorPow<Lt/nq,nq> ( T1 );
//    auto Gt=CyclicGroupPow<Lt>(T, Lt/nq);
//    auto G = Gt;

    bool Phi=par.phi;
    CadenitaAA5Open hnn(Lt/nq);
    {
        out << "Parameters\n";
        hnn.periodic=par.periodic;
        out << "periodic="<<hnn.periodic;
        hnn.mu=par.mu;
        out << " mu="<<hnn.mu;
        hnn.J=par.J;
        out << " J="<<hnn.J;
        hnn.muCu=par.muCu;
        out << " muCu="<<hnn.muCu;
        hnn.tp=par.tp;
        out << "\ntp="<<hnn.tp;
        hnn.tpp=par.tpp;
        out << " tpp="<<hnn.tpp;
        hnn.tOO=par.tOO;
        out << " tOO="<<hnn.tOO;
        hnn.U3=par.U3;
        out << " U3="<<hnn.U3;
        if ( Phi )
            hnn.phi=M_PI*nq/Lt; // <------ chequear con cuidado!!!
        else hnn.phi=0;
        out << "\nphi="<<hnn.phi*Lt/nq;
    }

    hnn.Initialize();
    int sz=par.Sz;
    int Sz=2*sz; // <--------- Ojo, si el numero de particulas es impar, esto no funciona
      out << " Sz="<<Sz;


    cout<<setprecision(15);

//    auto G=hnn.SymTraslation<Lt>();
    auto G=hnn.SymReflectionOnSite<Lt>();
    EigenStateG<double> gs;

//    out<<setprecision(15)<<" ang_spin*Lt/nq "<<ang_spin*Lt/nq<<"\n";
    out<<"\n";
    out<<"Energy\n";

    auto hasSz=hnn.HasSz<Lt>(Sz);
    for(int nu=0;nu<G.nSym();nu++)
    {
        auto b=FockBasisFixedChargeG<Lt>(par.nPart,G,nu,hasSz); //FockBasisFixedChargeG<Lt>(par.nPart,G,nu);
        cout<<"basis size="<<b.Size()<<endl; cout.flush();
        auto gsn=FindGS<Lt>(hnn.Ham(),b,G,par.nPart);
        out<< "nu "<< nu << " ";
        out<<setprecision(9)<<gsn.ener << "\n";
        gs=std::min( gs, gsn );
    }
    out<<"\n";
    //if ((id+1)%nTwist==0) out<<"\n\n";

//    vector<double> n_part(Lt);
    auto b=FockBasisFixedChargeG<Lt>(par.nPart,G,gs.sym,hasSz) ; //FockBasisFixedChargeG<Lt>(par.nPart,G,gs.sym);


//    for(int i=0;i<Lt;i++)
//    {
//        QOperator ni=Create(i)*Destroy(i);
//        auto niv=Dot(gs.state,  ni.toMatrix(b,G,b) * gs.state);
//        n_part[i]=real(niv);
//    }
    cout<<"GS: nPart="<<gs.nPart<<" sym="<<gs.sym<<" ener="<<gs.ener<<"\n";
//    for(int i=0;i<Lt;i++)
//        cout<<i<<" "<<n_part[i]<<endl;
//    gs.Print(b,1e-1);
}
