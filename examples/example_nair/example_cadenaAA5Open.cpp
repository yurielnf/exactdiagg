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
    auto Gr=hnn.SymReflectionY<Lt>();
    auto Gsp=hnn.SymSpinFlip<Lt>();
    auto G= Sz==0 ? Gr.DirectProd(Gsp) : Gr;
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

/* Este es el output de referencia cuando habia simetria de traslacion

Hello World!
basis size=9169
matrix nonzero=172319
-23.1512339094863
-23.1457921220032
-23.1438977728969
-19.4371776023077
-19.3540652928477
-19.3540652928476
-19.2658049992995
-19.1082816498055
-19.1082816498055
-19.0935888760008
-19.0905666419482
-18.9105792369143
dim=9169 nPart=4 sym=0 ener=-23.1512339094863
basis size=9169
matrix nonzero=168715
-23.1487133007212
-23.1460670405201
-23.1460670405201
-19.4149428044646
-19.3501348401211
-19.350134840121
-19.2400530258198
-19.1226318355738
-19.1008734647146
-19.1008734647146
-19.0901272171505
-18.9931716221981
-18.8897402095018
-18.8897402095017
dim=9169 nPart=4 sym=1 ener=-23.1487133007212
basis size=9169
matrix nonzero=170652
-18.9017767040417
-18.8864435024124
-18.8836752375494
-18.8836752375494
-18.8804642715322
-18.8683993993848
-18.8683993993848
-18.8322643993158
-18.7865292720648
-18.7865292720647
-18.7798067315794
-18.7709212584709
-18.7260668052498
dim=9169 nPart=4 sym=2 ener=-18.9017767040417
basis size=9169
matrix nonzero=170416
-18.9110592405034
-18.9084445058908
-18.9084445058908
-18.8947520821377
-18.8426193685099
-18.8243215668493
-18.8124853835775
-18.8124853835775
-18.8066575800347
-18.8066575800347
-18.7949669926111
-18.7827294370608
-18.7572084879024
dim=9169 nPart=4 sym=3 ener=-18.9110592405034
GS: nPart=4 sym=0 ener=-23.1512339094863


---------------------------- Y este con condiciones abiertas -----------------------------

Hello World!
basis size=9169
matrix nonzero=130747
-19.125636125758
-19.029988546021
-19.0180719924479
-18.6553353993989
-18.3219943493753
-18.2274363173825
-18.0347411957279
-17.9730785986892
-17.6272715742018
-17.3955830532868
dim=9169 nPart=4 sym=0 ener=-19.125636125758
basis size=9169
matrix nonzero=127964
-19.1050100369767
-19.0542471025218
-19.0090236963502
-18.6441206833301
-18.3785076290227
-18.1971787001959
-18.0182173979389
-17.9949615205755
-17.6112529818011
-17.3759006339667
dim=9169 nPart=4 sym=1 ener=-19.1050100369767
basis size=9169
matrix nonzero=129488
-18.661755438898
-18.6308740192071
-18.5948234071507
-17.6587159369247
-17.6220057651487
-17.6090645966068
-17.057335513471
-17.049745948653
-17.0353653517587
-17.033530373207
dim=9169 nPart=4 sym=2 ener=-18.661755438898
basis size=9169
matrix nonzero=129260
-18.6526350212339
-18.6433751041951
-18.5914927858142
-17.6501935038119
-17.6409422200759
-17.5985546672448
-17.0593000154596
-17.0487877468219
-17.0347906556225
-17.0315465378564
dim=9169 nPart=4 sym=3 ener=-18.6526350212339
GS: nPart=4 sym=0 ener=-19.125636125758


*/
