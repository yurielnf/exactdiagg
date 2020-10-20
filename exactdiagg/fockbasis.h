#ifndef FOCKBASIS_H
#define FOCKBASIS_H

#include<vector>
#include<map>
#include<unordered_map>
#include<bitset>
#include<iostream>
#include<limits>
#include<functional>
#include<cmath>
#include<complex>

//#include"binder.h"

template<int L, class scalar> class SymmetryGroup;

template<int L> using FockState=std::bitset<L> ;

template<int L1,int L2>
FockState<L1+L2> TensorProd(const FockState<L1>& f1,const FockState<L2>& f2)
{
    FockState<L1+L2> f1p ( f1.to_string() );
    FockState<L1+L2> f2p ( f2.to_string() );
    return (f1p<<L2) | f2p;
}

template<int L,int R>
std::array<FockState<L/R>,R> Split(FockState<L> f)
{
    static const int step=L/R;
    static const FockState<L> mask=( ulong(1) << step )-1;

    std::array<FockState<L/R>,R> a;
    for(int i=0;i<R;i++)
    {
        a[i] = (f & mask).to_ulong();
        f>>=step;
    }
    return a;
}

template<int L1,int R>
FockState<L1*R> Join(const std::array<FockState<L1>,R>& a)
{
    FockState<L1*R> f;
    for(int i=0;i<R;i++)
    {
        FockState<L1*R> mask(a[i].to_ulong());
        f|=( mask << (i*L1) );
    }
    return f;
}


template<int L>
bool BitParityLeft(FockState<L> f,int j)
{
    f>>=(j+1);
    return f.count()&1;
}


template<int L>
struct FockLess {
  bool operator() (const FockState<L>& f1,const FockState<L>& f2) const
  {return f1.to_ullong()<f2.to_ullong();}
};

template<int L>
bool operator< (const FockState<L>& f1,const FockState<L>& f2)
{return f1.to_ullong()<f2.to_ullong();}

template<int L,class scalar>
using State=std::map<FockState<L>,scalar,FockLess<L>>;


template<int L,class scalar>
double Norm(const State<L,scalar>& s)
{
    double sum=0;
    for(const auto& it:s)
        sum+=std::norm(it.second);
    return std::sqrt(sum);
}


//----------------------------------------------------------- FockBasis --------------------------

template<int L>
struct FockBasis
{
    typedef FockState<L> FState;
    std::vector<FState> vec;
    std::unordered_map<FState,int> mapa;
    int sym=0; //symmmetry sector
    std::vector<double> norma;


    int Size() const {return vec.size();}

    void Add(const FState& f,double normf)
    {
        if (mapa.find(f)!=mapa.end()) return;
        vec.push_back(f);
        mapa[f]=vec.size()-1;
        norma.push_back( FixNorm(normf) );
    }

    template<class scalar>
    void Add(FState f,const SymmetryGroup<L,scalar>& G)
    {
        if (G.nSym()==0) {Add(f,1.0); return;}
        State<L,scalar> s=G.ApplyTo(f,sym);
        double nr=Norm<L>(s);
        Add(s.begin()->first,nr);
    }

    void GenerateBaseFull()
    {
        ulong n=1<<L;
        for(ulong i=0;i<n;i++)
            Add(i,1.0);
    }

    template<class scalar>
    void GenerateBaseFull(const SymmetryGroup<L,scalar>& G)
    {
        ulong n=1<<L;
        for(ulong i=0;i<n;i++)
            Add(i,G);
    }

    template<class scalar>
    void SetSym(const SymmetryGroup<L,scalar>& G,int sym)
    {

        this->sym=sym;
        for(uint i=0;i<vec.size();i++)
        {
            State<L,scalar> s=G.ApplyTo(vec[i],sym);
            double nr=Norm<L>(s);
            norma[i]=FixNorm(nr);
        }
    }

    int PosOf(const FState& f) const {
        auto it=mapa.find(f);
        if (it!=mapa.end())
            return it->second;
        else return -1;
    }

    void Print() const
    {
        for(int i=0;i<vec.size();i++)
            std::cout<<norma[i]<<"("<<vec[i]<<") ";
        std::cout<<"\n";
    }

private:
    static double FixNorm(double nr)
    {
        static const double eps=sqrt(std::numeric_limits<double>::epsilon());
        return fabs(nr)<=eps? 0.0:(1.0/nr) ;
    }
};

template<int L1,int L2>
FockBasis<L1+L2> TensorProd(const FockBasis<L1>& b1,const FockBasis<L2>& b2)
{
    FockBasis<L1+L2> br;
    for(const FockState<L1>& f:b1.vec)
        for(const FockState<L2>& f2:b2.vec)
            br.Add(TensorProd<L1,L2>(f,f2));
    return br;
}

template<int L1,int L2,class SymmetryGroup>
FockBasis<L1+L2> TensorProd(const FockBasis<L1>& b1,const FockBasis<L2>& b2,const SymmetryGroup& G)
{
    FockBasis<L1+L2> br;
    for(const FockState<L1>& f1:b1.vec)
        for(const FockState<L2>& f2:b2.vec)
            br.Add(TensorProd<L1,L2>(f1,f2),G);
    return br;
}


//---------------------------------------------------------- FockBasisFixedCharge --------------------

inline long Combinatorial(int n, int k)
{
    if (k > n/2)  return Combinatorial(n, n-k);
    long cnk = 1;
    for(int i = 1; i <= k; ++i)
    {
        cnk *= n-i+1;
        cnk /= i;
    }
    return cnk;
}


template<int L>
struct FockBasisFixedCharge:public FockBasis<L>
{
    int nPart;
    FockBasisFixedCharge(int nPart)
        :nPart(nPart)
    {
        FockBasis<L> b;
        long dim=Combinatorial(L,nPart);
        b.vec.reserve(dim);
        b.mapa.reserve(dim);
        FockState<L> f;
        GeneraBase(f,0,nPart);
    }

protected:
    void GeneraBase(FockState<L>& f, int confPos,int nPart)
    {
        if (L<nPart+confPos) return;
        if (nPart==0) this->Add(f,1.0);
        else
        {
            GeneraBase(f,confPos+1,nPart);
            f[confPos]=1;
            GeneraBase(f,confPos+1,nPart-1);
            f[confPos]=0;
        }
    }
};


//---------------------------------------------------------- FockBasisFixedCharge with Group --------------------
#include"symmetrygroup.h"

template<int L>
struct FockBasisFixedChargeG:public FockBasis<L>
{
    int nPart;
//    const SymmetryGroup<L,scalar>& G;
    template<class scalar>
    FockBasisFixedChargeG(int nPart,const SymmetryGroup<L,scalar>& G,int sym=0)
        :nPart(nPart)
    {
        this->sym=sym;
//        FockBasis<L> b;
//        long dim=Combinatorial(L,nPart);
//        b.vec.reserve(dim);
//        b.mapa.reserve(dim);
        FockState<L> f;
        GeneraBase(f,0,nPart,G);
    }
    template<class scalar>
    FockBasisFixedChargeG(int nPart,const SymmetryGroup<L,scalar>& G,int sym,
                          std::function<bool(const FockState<L>&)> isOk)
        :nPart(nPart)
    {
        this->sym=sym;
        FockState<L> f;
        GeneraBase(f,0,nPart,G,isOk);
    }

protected:
    template<class scalar>
    void GeneraBase(FockState<L>& f, int confPos,int nPart,const SymmetryGroup<L,scalar>& G)
    {
        if (L<nPart+confPos) return;
        if (nPart==0) this->Add(f,G);
        else
        {
            GeneraBase(f,confPos+1,nPart,G);
            f[confPos]=1;
            GeneraBase(f,confPos+1,nPart-1,G);
            f[confPos]=0;
        }
    }
    template<class scalar>
    void GeneraBase(FockState<L>& f, int confPos,int nPart,const SymmetryGroup<L,scalar>& G,
                    std::function<bool(const FockState<L>&)> isOk)
    {
        if (L<nPart+confPos) return;
        if (nPart==0) {
            if (isOk(f)) this->Add(f,G);
        }
        else
        {
            GeneraBase(f,confPos+1,nPart,G,isOk);
            f[confPos]=1;
            GeneraBase(f,confPos+1,nPart-1,G,isOk);
            f[confPos]=0;
        }
    }
};




#endif // FOCKBASIS_H
