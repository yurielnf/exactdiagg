#include "examples/all.h"
#include"exactdiagg/all.h"

#include <iostream>
#include <complex.h>
#include "parameters.h"

using namespace std;

double tb_exact(int L){
    double e=0, k, ek;
    for(int a=0;a<L;a++){
        k=2*M_PI*a/L;
        ek=-2*cos(k);
        if(ek<0){
            e+=ek;
        }
    }
    return(e);
}



std::string order="OSP";


const int lonx=5;
const int lenx=4*lonx;

int toIdx(int pos, int orb, int spin)
{
    int idx=0;

    if(order=="POS"    ) idx= 4*pos + 2*orb + spin      ;
    else if(order=="OSP"    ) idx=   pos +(2*orb + spin)*lonx;
    else {}; // something wrong!!

    return idx;
}

inline FermiOp Fock(int pos, int orb, int spin, bool occ     ) { return FermiOp(toIdx(pos, orb, spin), occ);  }




QOperatorG<cmpx> Ham_KH_fluxed(int L, double t1, double t2, double U, double V, double J, double dmu, double phi)
{
    cout<<"building KH model with magnetic flux"<<endl;

    QOperatorG<cmpx> H;

    complex<double> imag(0.0, 1.0);


    double e0 = (-0.5*U -V +0.5*J), mu = dmu - e0;

    for(int pos=0; pos<L; pos++)
    {
        for(int orb=0; orb<2; orb++) H.Add( Fock(pos, orb, 0, true)*Fock(pos, orb, 0, false)
                                           *Fock(pos, orb, 1, true)*Fock(pos, orb, 1, false), U);

        for(int s1 =0; s1 <2;  s1++)
        for(int s2 =0; s2 <2;  s2++)
        {
            H.Add( Fock(pos, 0, s1, true)*Fock(pos, 0, s1, false)
                  *Fock(pos, 1, s2, true)*Fock(pos, 1, s2, false), V - ((s1==s2)? J:0) );

            if(s1!=s2 && J!=0)
            H.Add( Fock(pos, 0, s1, true)*Fock(pos, 0, 1-s1, false)
                  *Fock(pos, 1, s2, true)*Fock(pos, 1, 1-s2, false), -J);
        }

        if(J!=0)
        for(int orb=0; orb<2; orb++) H.Add( Fock(pos,   orb, 0, true )*Fock(pos,   orb, 1, true )
                                           *Fock(pos, 1-orb, 0, false)*Fock(pos, 1-orb, 1, false), -J);


        for(int  orb=0;  orb<2;  orb++) for(int spin=0; spin<2; spin++)
        {
            H.Add( Fock( pos     , orb, spin, true)*Fock((pos+1)%L, orb, spin, false), -((orb==0)? t1:t2) * exp( imag*(phi/L)) );
            H.Add( Fock((pos+1)%L, orb, spin, true)*Fock( pos     , orb, spin, false), -((orb==0)? t1:t2) * exp(-imag*(phi/L)) );

            if(mu!=0)
            H.Add( Fock( pos     , orb, spin, true)*Fock( pos     , orb, spin, false), -mu );
        }
    }

    return H;
}


template<int L>
SymmetryGroup<4*L,cmpx> KH_Group()
{
    const int Lt=4*L;

    auto    chain_T1 = TranslationOp<L> (1);
    auto   system_T1 = TensorPow<L,4> (chain_T1);
    auto trasl_group = CyclicGroupPow<Lt> (system_T1, L);

    auto   system_OF = TranslationOp<Lt> (2*L);
    auto    OF_group = Z2_Group<Lt, cmpx> (system_OF);

    auto      orb_SF = TranslationOp<2*L> (L);
    auto      sys_SF = TensorPow<2*L,2> (orb_SF);
    auto    SF_group = Z2_Group<Lt, cmpx> (sys_SF);

    auto  chain_rflx = ReflectionOp<L>;
    auto system_rflx = TensorPow<L, 4, ElementaryOp<L> > (chain_rflx);
    auto  rflx_group = Z2_Group<Lt, cmpx> (system_rflx);

    auto   system_PH = ParticleHoleOp<Lt>;
    auto    PH_group = Z2_Group<Lt, cmpx> (system_PH);



    auto    G = trasl_group;

            G = G.DirectProd(SF_group);
            G = G.DirectProd(OF_group);

            G = G.DirectProd(rflx_group); //only for \Phi=0
            G = G.DirectProd(PH_group);   //only for half-filled

    return G;
}



void KHC_N(Parameters param)
{
    const int L=lonx;
    const int Lt=lenx;

    int Le=L;

    double t1=param.t1, t2=param.t2, U=param.U, V=param.V, J=param.J, nPart=param.nPart, phi=param.phi;
    cout<<"L="<<L<<"\tt1="<<t1<<"\tt2="<<t2<<"\tU="<<U<<"\tV="<<V<<"\tJ="<<J<<"\tphi="<<phi<<"\tnPart="<<nPart<<endl;

    auto H = Ham_KH_fluxed(L, t1, t2, U, V, J, 0, phi);
    auto G = KH_Group<L>();

    cout<<"nSym= "<<G.nSym()<<endl<<endl;


        auto b=FockBasisFixedChargeG<Lt>(nPart, G);
        auto gs=FindGS<Lt>(H, b, G, nPart);

        for(int nu=1; nu<G.nSym(); nu++)
        {
            b.SetSym(G, nu);
            auto gsn = FindGS<Lt>(H, b, G, nPart);
            gs = std::min( gs, gsn );
        }
        cout<<"GS(nPart="<<gs.nPart<<") --> sym="<<gs.sym<<" ener="<<gs.ener<<endl<<endl;


    ofstream E0stream("E0.dat");
    E0stream<<gs.ener;
}

/*
void KHC_mu(Parameters param)
{
    const int L=lonx;
    const int Lt=lenx;

    int Le=L;

    double t1=param.t1, t2=param.t2, U=param.U, V=param.V, J=param.J, dmu=param.dmu, phi=param.phi;
    cout<<"L="<<L<<"\tt1="<<t1<<"\tt2="<<t2<<"\tU="<<U<<"\tV="<<V<<"\tJ="<<J<<"\tphi="<<phi<<"\tdmu="<<dmu<<endl;

    auto H = Ham_KH_fluxed(L, t1, t2, U, V, J, dmu, phi);
    auto G = KH_Group<L>();

    auto b0 = FockBasisFixedChargeG<Lt>(Lt/2, G, 0);
    auto EF = FindGS<Lt>(H, b0, G, Lt/2);
    cout<<"starting with half-filled E0="<<EF.ener<<endl;

    for(int nPart=0; nPart<=Lt; nPart++)
    {
        auto b=FockBasisFixedChargeG<Lt>(nPart, G, 0);
        auto gs=FindGS<Lt>(H,b,G,nPart);

        for(int nu=1; nu<G.nSym(); nu++)
        {
            b.SetSym(G, nu);
            auto gsn = FindGS<Lt>(H, b, G, nPart);
            gs = std::min( gs, gsn );
        }
        cout<<"nPart="<<gs.nPart<<" --> sym="<<gs.sym<<" ener="<<gs.ener<<endl<<endl;

        EF = std::min( EF, gs );
    }


    cout<<"GS: nPart="<<EF.nPart<<" sym="<<EF.sym<<" ener="<<EF.ener<<endl;

    ofstream E0stream("E0.dat");
    E0stream<<EF.ener;
}
*/


void Hernan_ED_gs(const char filename[])
{
    auto param=Parameters(filename);
//    if(param.model=="KHC" ) KHC_mu(param);
    if(param.model=="KHCN") KHC_N (param);

}
