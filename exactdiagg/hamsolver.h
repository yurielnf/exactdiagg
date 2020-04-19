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
        for(size_t i=0;i<state.size();i++)
            if(fabs(state(i))>tol)
                cout<<b.vec[i]<<" "<<state(i)<<endl;
    }
};

template<class scalar>
bool operator<(const EigenStateG<scalar>& e1,const EigenStateG<scalar>& e2)
{
    return e1.ener<e2.ener;
}

using EigenState=EigenStateG<>;

template<int L>
EigenState FindGS(QOperator ham,FockBasisFixedChargeG<L>& b,
                  const SymmetryGroup<L,double>& G, int nPart)
{
    auto H=ham.toMatrix<L>(b,G);

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
EigenStateG<cmpx> FindGS(QOperatorG<cmpx> ham,FockBasisFixedChargeG<L>& b,
                  const SymmetryGroup<L,cmpx>& G, int nPart)
{
//        b.SetSym(G,sym);
    auto H=ham.toMatrix<L>(b,G);

    cx_vec evalz;
    cx_mat evec;
    eigs_gen(evalz,evec,H ,std::min(b.Size()-2,L*2),"sr");
    vec eval=arma::real(evalz);
    uvec ind = sort_index(eval);
    for(uint i=0;i<eval.size();i++)
        if (eval(i)-eval(ind(0))<1e-10)
            cout<<eval[i]<<"\n";
    cout<<"dim="<<b.Size()<<" nPart="<<nPart<<" sym="<<b.sym<<" ener="<<eval(ind(0))<<endl;
    return EigenStateG<cmpx>(eval(ind(0)),evec.col(ind(0)),nPart,b.sym);
}

#endif // HAMSOLVER_H
