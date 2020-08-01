#ifndef HAMSOLVER_H
#define HAMSOLVER_H

#include<iostream>
#include<iomanip>
#include<vector>
#include<complex>
#include<fstream>
#include<string>
#include<armadillo>

#include"fockbasis.h"
#include"qoperator.h"
#include"symmetrygroup.h"


template<class scalar=double>
struct EigenStateG
{
    double ener=std::numeric_limits<double>::max();
    Col<scalar> state;
    int nPart;
    int sym;

    EigenStateG(){}
    EigenStateG(double ener,const Col<scalar>& state,int nPart,int sym)
        :ener(ener),state(state), nPart(nPart), sym(sym) {}

    void Save(std::string filename) const
    {
        std::ofstream out(filename);
        out<<std::setprecision(14);
        out<<ener<<" ";
        out<<nPart<<" ";
        out<<sym<<endl;
        state.save(out);
    }
    bool Load(std::string filename)
    {
        std::ifstream in(filename);
        if (!in) return false;
        in>>ener;
        in>>nPart;
        in>>sym;
        state.load(in, arma::arma_binary);
        return true;
    }
    template<int L>
    void Print(const FockBasis<L> &b,double tol=1e-3) const
    {
        uvec id=sort_index(abs(state),"descend");
        for(size_t i=0;i<state.size();i++)
            if(fabs(state(id(i)))>tol)
            {
                auto p=id(i);
                cout<<b.vec[p]<<" "<< abs(state(p))<<" "<<state(p)<<endl;
            }
    }
};

template<class scalar>
bool operator<(const EigenStateG<scalar>& e1,const EigenStateG<scalar>& e2)
{
    return e1.ener<e2.ener;
}

using EigenState=EigenStateG<>;

template<int L,class T>
EigenState FindGS(QOperatorG<T> ham,FockBasisFixedChargeG<L>& b,
                  const SymmetryGroup<L,double>& G, int nPart)
{
    auto H=ham.template toMatrix<L>(b,G);

    vec eval;
    mat evec;
    eigs_sym(eval,evec,H ,std::min(b.Size()-2,L*2),"sa");
    uvec ind = sort_index(eval);
    for(uint i=0;i<eval.size();i++)
        if (eval(i)-eval(ind(0))<1e-10)
            cout<<eval[i]<<"\n";
    cout<<"dim="<<b.Size()<<" nPart="<<nPart<<" sym="<<b.sym<<" ener="<<eval(ind(0))<<endl;
    return {eval(ind(0)),evec.col(ind(0)),nPart,b.sym};
}


template<int L>
EigenStateG<cmpx> FindGS(QOperator ham,FockBasisFixedChargeG<L>& b,
                  const SymmetryGroup<L,cmpx>& G, int nPart)
{
//        b.SetSym(G,sym);
    auto H=ham.toMatrix<L,cmpx>(b,G);

    cx_vec evalz;
    cx_mat evec;
    eigs_gen(evalz,evec,H ,std::min(b.Size()-2,L*2),"sr");
    vec eval=arma::real(evalz);
    uvec ind = sort_index(eval);
    for(uint i=0;i<eval.size();i++)
        if (eval(i)-eval(ind(0))<1e-10)
            cout<<eval[i]<<"\n";
    cout<<"dim="<<b.Size()<<" nPart="<<nPart<<" sym="<<b.sym<<" ener="<<eval(ind(0))<<endl;
    return {eval(ind(0)),evec.col(ind(0)),nPart,b.sym};
}


template<int L>
EigenState FindGS(QOperator ham,FockBasisFixedCharge<L>& b,int nPart)
{
    auto H=ham.toMatrix<L>(b);

    vec eval;
    mat evec;
    eigs_sym(eval,evec,H ,std::min(b.Size()-2,100),"sa");
    uvec ind = sort_index(eval);
    for(uint i=0;i<eval.size();i++)
        if (eval(i)-eval(ind(0))<1e-10)
            cout<<eval[i]<<"\n";
    cout<<"dim="<<b.Size()<<" nPart="<<nPart<<" sym="<<b.sym<<" ener="<<eval(ind(0))<<endl;
    return {eval(ind(0)),evec.col(ind(0)),nPart,b.sym};
}




//----------------------------------------- Time evolution---------------------


struct TimeEvolution
{
    cx_mat evec;
    vec eval;
    cx_vec psi0_n; //to calculate <n|psi0>
    TimeEvolution(const sp_cx_mat &H, const cx_vec &psi0)
        :psi0_n(psi0.size())
    {
        eig_sym(eval,evec,cx_mat(H));
        for(uint n=0;n<eval.size();n++)
            psi0_n[n]=cdot(evec.col(n),psi0);
    }
    cx_vec EvolveTo(double t) const
    {
        cx_vec gt(eval.size());
        for(uint n=0;n<eval.size();n++)
        {
            cmpx q={0,-eval[n]*t};
            gt(n)=psi0_n[n]*exp(q);
        }
        return evec*gt;
    }
};



#endif // HAMSOLVER_H
