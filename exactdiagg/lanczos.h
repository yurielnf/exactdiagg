#ifndef LANCZOS_H
#define LANCZOS_H

#include<iostream>
#include<vector>
#include<stdexcept>
#include<cmath>
#include<complex>
#include<armadillo>

using cmpx=std::complex<double>;

inline double Norm(const arma::vec& x) {return arma::norm(x);}
inline double Norm(const arma::cx_vec& x) {return arma::norm(x);}
inline double Dot(const arma::vec& x,const arma::vec& y) {return arma::cdot(x,y);}
inline cmpx Dot(const arma::cx_vec& x,const arma::cx_vec& y) {return arma::cdot(x,y);}

//---------------------- binder -------------


template<class scalar>
struct EigenPair
{
    std::vector<double> eval;
    std::vector<scalar> evec;
};

inline EigenPair<double> DiagonalizeTridiagonal(double *an, double *bn, int size,int nEvals=1)
{
    arma::sp_mat h(size,size);
    h.diag()=arma::vec(an,size);
    if (size>=2) {
        h.diag(-1)=arma::vec(bn,size-1);
        h.diag(1)=h.diag(-1).t();
    }
    arma::vec eval;
    arma::mat evec;
    if (size>10)
        arma::eigs_sym(eval,evec,h ,std::min(size-2,101),"sa");
    else
        eig_sym(eval,evec,arma::mat(h));
    //uvec ind = sort_index(eval);                          <--- TODO: sort eval
    return {{eval.begin(),eval.begin()+nEvals},
            {evec.begin(),evec.begin()+nEvals*size}};
}

inline EigenPair<cmpx> DiagonalizeTridiagonal(cmpx *an,cmpx *bn, int size,int nEvals=1)
{
    arma::sp_cx_mat h(size,size);

    h.diag()=arma::cx_vec(an,size);
    if (size>=2) {
        h.diag(-1)=arma::cx_vec(bn,size-1);
        h.diag(1)=h.diag(-1).t();
    }
    if (!h.is_hermitian(1e-14) || arma::norm(arma::imag(h))>1e-10)
    {
        h.print("h=");
        throw std::invalid_argument("lanczos non Hermitian");
    }
    arma::vec eval;
    arma::cx_mat evec;
//    if (size>10)
//    {
//        arma::cx_vec evalz;
//        arma::eigs_gen(evalz,evec,h,std::min(size-2,101),"sr");
//        eval=arma::real(evalz);
//    }
//    else
        arma::eig_sym(eval,evec,arma::cx_mat(h));
    //arma::uvec ind = arma::sort_index(eval);                        //  <--- TODO: sort eval
    return {{eval.begin(),eval.begin()+nEvals},
            {evec.begin(),evec.begin()+nEvals*size}};
}




//----------------------Lanczos ------------
template<class Ket>
void Orthogonalize(Ket &t,const std::vector<Ket>& basis,int size)
{
    const double k=0.25;
    double tauin=Norm(t);
    for (int i = 0; i < size; ++i)
        t-=basis[i]*Dot(basis[i],t);
    if ( Norm(t)/tauin > k ) return;
    for (int i = 0; i < size; ++i)
        t-=basis[i]*Dot(basis[i],t);
}




template<class Operator, class Ket,class Num>
struct Lanczos
        //Find the lowest eigenvalue lambda0, and eigenvector x0
        //for the eigen-problem A x = lambda x, A Hermitian
{
    const Operator& A;
    std::vector<Ket> basis;
    Ket x0;                             //ground state
    double ener0, error;
    std::vector<double>  an,bn,evec;
    int nInnerIter,nMaxIter,cIter;

    double tol;

    const double orthogonality=1e-15;

    Lanczos(const Operator& A,const Ket& v0,int nInnerIter,int nMaxIter,double tol):
        A(A),x0(v0),nInnerIter(nInnerIter),nMaxIter(nMaxIter),cIter(0),tol(tol)
    {}
    int Iterate(int nIter)              //It returns the number of iterations
    {
        an.resize(nIter);
        bn.resize(nIter);
        basis.resize(nIter);

        Ket r=x0;
        bn[0]=Norm(r);
        for(int j=0;j<nIter-1;j++)
        {
            basis[j]=r*(1.0/bn[j]);
            r=A*basis[j]; cIter++;
            if (j>0) r-=basis[j-1]*bn[j];
            auto dt=Dot(basis[j],r);
            if (std::imag(dt)>tol) throw std::invalid_argument("Lanczos: non Hermitian operator");
            an[j]=std::real(dt);
            r-=basis[j]*an[j];
            Orthogonalize(r,basis,j);
            bn[j+1]=Norm(r);
            auto eigen=DiagonalizeTridiagonal(an.data(),bn.data()+1,j+1);
            error=fabs(bn[j+1]*eigen.evec[j]);
            if (cIter==nMaxIter || error<tol || std::fabs(bn[j+1])<orthogonality)
            {
                an.resize(j+1);
                bn.resize(j+1);
                basis.resize(j+1);

                ener0=eigen.eval[0];
                evec=eigen.evec;
                return an.size();
            }
        }
        int j=an.size()-1;
        basis[j]=r*(1.0/bn[j]);
        r=A*basis[j]-basis[j-1]*bn[j]; cIter++;
        auto dt=Dot(basis[j],r);
        if (std::imag(dt)>tol) throw std::invalid_argument("Lanczos: non Hermitian operator");
        an[j]=std::real(dt);

        auto eigen=DiagonalizeTridiagonal(an.data(),bn.data()+1,j+1);
        ener0=eigen.eval[0];
        evec=eigen.evec;
        return an.size();
    }

    void BuildGroundState()
    {
        x0-=x0;
        for(uint i=0;i<an.size();i++)
            x0+=basis[i]*evec[i];
        x0*=1.0/Norm(x0);
    }

    void DoIt()
    {
        double error2=1;
        while(cIter < nMaxIter)
        {
            Iterate(nInnerIter);


            error2=fabs( pow(Dot(x0,A*x0),2)-Dot(x0,A*(A*x0))) ;
            if (error2<tol*tol || error<tol) break;
        }
        BuildGroundState();
        if (error2>tol*tol && error>tol)
        {
            std::cout<<"lanczos failed, error2 = "<<error2<<" ";
            std::cout<<"residual_norm = "<<error<<std::endl;
        }
    }
};

template<class scalar,class Hamiltonian, class Ket>                                      //Portal method
Lanczos<Hamiltonian,Ket,scalar> Diagonalize(const Hamiltonian& H,Ket& wf,int nIter=5000,double tol=1e-14)
{
    Lanczos<Hamiltonian,Ket,scalar> lan(H,wf,std::min(128,nIter), nIter, tol);
    lan.DoIt();
    return lan;
}


#endif // LANCZOS_H

