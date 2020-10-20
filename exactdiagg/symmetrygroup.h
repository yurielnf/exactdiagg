#ifndef SYMMETRYGROUP_H
#define SYMMETRYGROUP_H

//---------------------------------------------------------- ElementaryOp -------------------------------------

#include"fockbasis.h"
#include<armadillo>
#include<complex>

using cmpx=std::complex<double>;

using namespace arma;

template<int L>
using ElementaryOp = std::function<int(FockState<L>&)> ;

template<int L>
int IdentityOp(FockState<L>& f){return 1;}

template<int L>
struct TranslationOp
{
    int delta;

    TranslationOp(int delta):delta(delta){}

    int operator()(FockState<L>& f) const
    {
        if (delta==0) return 1;
        FockState<L> l=f>>delta;
        FockState<L> r=f<<(L-delta);
        f=l|r;
        bool parity=(r.count() * l.count()) & 1;//(r.count() & 1) && !(f.count() & 1);
        if (parity) return -1;
        else return 1;
    }
    TranslationOp Pow(int n) const { return TranslationOp(delta*n); }
};

template<int L>
int ReflectionOp(FockState<L>& f)
{
    std::string s=f.to_string();
    std::reverse(s.begin(), s.end());
    f=FockState<L>(s);
    if ( (f.count()/2) & 1) return -1;
    else return 1;
}

template<int L>
int ParticleHoleOp(FockState<L>& f)
{
    f=~f;
    if (f.count() & 1) return -1;
    else return 1;
}


template<int L>
struct ProductOp
{
    ElementaryOp<L> op1;
    ElementaryOp<L> op2;
    ProductOp(ElementaryOp<L> op1,ElementaryOp<L> op2): op1(op1),op2(op2) {}
    int operator()(FockState<L>& f) const
    {
        return  (op2(f)==-1) ? -op1(f) : op1(f);
    }
};

template<int L,int R, class ElementaryOp>
struct TensorPowOp
{
    ElementaryOp op;
    TensorPowOp(ElementaryOp op): op(op) {}
    int operator()(FockState<L*R>& f) const
    {
        auto a=Split<L*R,R>(f);
        int sg=1;
        for(int i=0;i<R;i++)
            if (op(a[i])==-1) sg=-sg;
        f=Join<L,R>(a);
        return sg;
    }
    TensorPowOp Pow(int n) const { return TensorPowOp(op.Pow(n)); }
};

template<int L, int R, class ElementaryOp>
TensorPowOp<L,R,ElementaryOp> TensorPow(const ElementaryOp& op) {return TensorPowOp<L,R,ElementaryOp>(op);}



template<class Oper,int L, class scalar=double>
State<L,scalar> ApplyOp(const Oper& Op,const State<L,scalar>& s)
{
    State<L,scalar> sr;
    for(const auto& it:s)
    {
        FockState<L> r=it.first;
        scalar coef=it.second;
        int sg=Op(r);
        sr[r]+=coef*scalar(sg);
    }
    return sr;
}


template<class Oper,int L,class scalar=double>
SpMat<scalar> toMatrix(const Oper& Op,const FockBasis<L> &b1, const SymmetryGroup<L,scalar>& G,const FockBasis<L> &b2)
{
    SpMat<scalar> M(b1.Size(),b2.Size());
    for(int i=0;i<b2.Size();i++)
    {
        if  (b2.norma[i]==0.0) continue;
        State<L,scalar> ssym=G.ApplyTo(b2.vec[i],b2.sym);
        State<L,scalar> res=ApplyOp<Oper,L,scalar>(Op,ssym);
        for(const auto& it:res)
        {
            FockState<L> r=it.first;
            scalar gm=G.Representant(r,b1.sym);
            int pos=b1.PosOf(r);
            if (pos!=-1 && b1.norma[pos]!=0.0)
                M(pos,i)+=gm*it.second*b1.norma[pos]*b2.norma[i];
        }
    }
    return M;
}


//---------------------------------------------------------- SymmetryGroup ----------------------------

#include<functional>

template<int L, class scalar>
class SymmetryGroup
{
public:
    std::vector< ElementaryOp<L> > T;
    std::vector< Col<scalar> > Gamma;

    void Add(ElementaryOp<L> op, const Col<scalar>& coeff)
    {
        T.push_back(op);
        Gamma.push_back(coeff);
    }

    int nSym() const {return Gamma.size()==0?0:Gamma[0].size();}

    State<L,scalar> ApplyTo(const FockState<L>& f,int sym) const
    {
        std::vector<FockState<L>> vec(T.size());
        std::vector<int> sg(T.size());
//        #pragma omp parallel for
        for(uint g=0;g<T.size();g++)
        {
            vec[g]=f;
            sg[g]=T[g](vec[g]);
        }
        State<L,scalar> state;
        for(uint g=0;g<T.size();g++)
        {
            scalar& value=state[vec[g]];
            if (sg[g]==1) value+=Gamma[g][sym];
            else          value-=Gamma[g][sym];
//            if (value==0.0) state.erase(vec[g]);
        }
        return state;
    }

    State<L,scalar> ApplyTo(const State<L,scalar>& s,int sym) const
    {
        State<L,scalar> state;
        for(auto it:s)
        {
            State<L,scalar> tmp=ApplyTo(it.first,sym);
            for(auto it2:tmp)
                state[it2.first]+=it.second*it2.second;
        }
        return state;
    }

    scalar Representant(FockState<L>& f, int sym) const
    {
        State<L,scalar> res=ApplyTo(f,sym);
        if (res.empty()) {std::cerr<<"Representant not found!"; return 0;}
        f=res.begin()->first;
        return res.begin()->second;
    }

    template<class scalar2>
    SymmetryGroup DirectProd(const SymmetryGroup<L,scalar2>& G2) const
    {
        const SymmetryGroup& G1=*this;
        if (G1.nSym()==0)
            throw std::invalid_argument("DirectProd: invalid G1");
        SymmetryGroup G;
        for(uint i=0;i<G1.Gamma.size();i++)
            for(uint j=0;j<G2.Gamma.size();j++)
            {
                auto g = [=](FockState<L>& f) { return (G2.T[j](f)==-1)?-G1.T[i](f):G1.T[i](f); };
//                auto g=ProductOp<L>(G1.T[i],G2.T[j]);
                G.Add(g,kron(G1.Gamma[i],G2.Gamma[j]));
            }
//         G.SetSym( G1.currSym*G2.nSym()+G2.currSym );
         return G;
    }
};

template<int L,class ElementaryOp>
SymmetryGroup<L,cmpx> CyclicGroupPow(const ElementaryOp& T,int nElem)
{
    SymmetryGroup<L,cmpx> G;
    for(int i=0;i<nElem;i++)
    {
        Col<cmpx> coeff(nElem);
        for(int k=0;k<nElem;k++)
        {
            cmpx phase={0.0, -2*M_PI*i*k/nElem};
            coeff[k]=std::exp(phase) * (1.0/nElem) ;
        }
        G.Add(T.Pow(i),coeff);
    }
    return G;
}

#include"fockbasis.h"

template<int L, class scalar=double,class ElementaryOp>
SymmetryGroup<L,scalar> Z2_Group(const ElementaryOp& T)
{
    SymmetryGroup<L,scalar> G;
    G.Add(IdentityOp<L>,{0.5,0.5});
    G.Add(T,{0.5,-0.5});
    return G;
}



#endif // SYMMETRYGROUP_H
