#ifndef CADENITAAA5OPEN_H
#define CADENITAAA5OPEN_H

#include"exactdiagg/all.h"
#include<string>
#include<fstream>
#include<iostream>

struct CadenitaAA5Open
{
    static const int nOrb=5;
    int L; //cantidad de celdas unidad
    double t=1.1, tp=0.1, tpp=0.1, tOO=-0.32, U=2, U3=7, J=0.838, muCu=0.5; //2.595;
    bool periodic=true;
    double mu=0,phi=0/*,angle_spin=0*/;
private:
    arma::mat delta, hop, Umat;
public:
    CadenitaAA5Open(int L)
        :L(L)
        ,delta(nOrb,nOrb,fill::zeros)
        ,hop(nOrb,nOrb,fill::zeros)
        ,Umat(2*nOrb,2*nOrb,fill::zeros)
    {}
    void Initialize()
    {
        delta(1,2)=delta(2,3)=delta(2,1)=delta(3,2) = t;
        delta(0,2)=delta(2,4)=delta(2,0)=delta(4,2) = tp;
        delta(2,2)=-muCu;
        for(int i=0;i<nOrb;i++) delta(i,i)-=mu;
        hop(0,2)=hop(4,2) = t;
        hop(1,2)=hop(3,2) = tp;
        hop(2,2)=tpp;
        hop(0,0)=hop(0,1)=hop(1,0)=hop(1,1)=tOO;
        hop(3,3)=hop(3,4)=hop(4,3)=hop(4,4)=tOO;
        for(int i=0;i<nOrb;i++)
            Umat(toInt(i,0),toInt(i,1))= (i==2) ? U3 : U;
        for(int s=0;s<2;s++)
            for(int sp=0;sp<2;sp++)
            {
                double c= (s==sp)? U-3*J : U-2*J;
                Umat(toInt(0,s),toInt(1,sp))=c;
                Umat(toInt(3,s),toInt(4,sp))=c;
            }
        // now for the first 4 orbitals

    }
//    int toInt(int i, int Ii, int spin) const { return spin+Ii*2+i*nOrb*2; }
//    int toInt(int Ii, int spin) const { return spin+Ii*2; }

//    int toInt(int i, int Ii, int spin) const { return i+Ii*L+spin*L*nOrb; }
    int toInt(int i, int Ii, int spin) const {

        if (i>=0)
            return 4+Ii+i*nOrb+spin*(4+nOrb*L);
        else
            return Ii<2 ? Ii+spin*(4+nOrb*L)
                        : Ii-1+spin*(4+nOrb*L);

    }
    int toInt(int Ii, int spin) const { return Ii+spin*nOrb; }

    QOperatorG<double> Kin() const
    {
        QOperatorG<double> h;
        for(int i=-1;i<L-1; i++)
            for(int ii=0;ii<nOrb;ii++)
                for(int jj=0;jj<nOrb;jj++)
                {
                    if (i==-1 && (ii==2 || jj==2)) continue;
                    for(int s=0;s<2;s++)
                        for(int sp=0;sp<2;sp++)
                        {
//                            double tt=hop(toInt(ii,s), toInt(jj,sp));
                            double tt=hop(ii, jj);
                            if (tt!=0)
                            {
                                int pi=toInt(i,ii,s);
                                int pj=toInt(i+1,jj,sp);
                                h.Add( FermiOp(pi,true)*FermiOp(pj,false), tt ); //tt*exp( cmpx{0,phi})
                                h.Add( FermiOp(pj,true)*FermiOp(pi,false), tt ); //tt*exp(-cmpx{0,phi})
                            }
                        }
                    double ee=delta(ii,jj);
                    if (ee!=0)
                        for(int s=0;s<2;s++)
                        {
                            int pi=toInt(i,ii,s);
                            int pj=toInt(i,jj,s);
                            h.Add( FermiOp(pi,true)*FermiOp(pj,false), ee );
                        }
                }
        return h;
    }
    QOperatorG<double> Pot() const
    {
        QOperatorG<double> h;
        for(int i=-1;i<L; i++)
            for(int ii=0;ii<nOrb;ii++)
                for(int s=0;s<2;s++)
                    for(int jj=0;jj<nOrb;jj++)
                        for(int sp=0;sp<2;sp++)
                        {
                            if (i==-1 && (ii==2 || jj==2)) continue;
                            int pi=toInt(i,ii,s);
                            int pj=toInt(i,jj,sp);
                            double coeff=Umat(toInt(ii,s),toInt(jj,sp));
                            if (coeff==0) continue;
                            h.Add( FermiOp(pi,true)*FermiOp(pi,false)*           //U n n
                                   FermiOp(pj,true)*FermiOp(pj,false), coeff);
                        }
        if (J!=0)
            for(int d=0;d<=3;d+=3)
                for(int i=-1;i<L; i++)
                    for(int ii=0;ii<2;ii++)
                    {
                        int pi=toInt(i,d+ii,0);
                        int pj=toInt(i,d+1-ii,1);
                        int pk=toInt(i,d+1-ii,0);
                        int pl=toInt(i,d+ii,1);
                        h.Add( FermiOp(pi,true)*FermiOp(pj,true)*           // -J c+c+ cc
                               FermiOp(pk,false)*FermiOp(pl,false), -J );
                        pi=toInt(i,d+ii,0);
                        pj=toInt(i,d+ii,1);
                        pk=toInt(i,d+1-ii,1);
                        pl=toInt(i,d+1-ii,0);
                        h.Add( FermiOp(pi,true)*FermiOp(pj,true)*           // J c+c+ cc
                               FermiOp(pk,false)*FermiOp(pl,false), J );
                    }
        return h;
    }
    QOperatorG<double> Ham() const { return Kin()+Pot(); }

    template<int Lt>
    std::function<bool(FockState<Lt>)> HasSz(int Sz) const
    {
        return [=](const FockState<Lt>& f) {
            int sum=0;
            for(int i=-1;i<L;i++)
                for(int ii=0;ii<nOrb;ii++)
                    if (i!=-1 || ii!=2)
                        sum+=f.test(toInt(i,ii,0))-
                                f.test(toInt(i,ii,1));
            return sum==Sz;
        };
    }

    /*template<int Lt>
    static SymmetryGroup<Lt,cmpx> SymTraslation()
    {
        const int L=Lt/(2*nOrb);
        auto T1=TranslationOp<Lt/2>(nOrb);
        auto T=TensorPow<Lt/2,2> ( T1 );
        return CyclicGroupPow<Lt>(T, L);
    }*/

    template<int L>
    static int ReflectionOpY(FockState<L>& f)
    {
        std::string s=f.to_string();
        auto s1=s.substr(0,4);
        auto s2=s.substr(4);

        auto f1=FockState<4>(s1);
        int sg1=ReflectionOp<4>(f1);

        static auto Refl=TensorPow<nOrb,(L-4)/nOrb,ElementaryOp<nOrb>> ( ReflectionOpX<nOrb> );
        auto f2=FockState<L-4>(s2);
        int sg2=Refl(f2);

        f=FockState<L>(f1.to_string()+f2.to_string());
        return sg1*sg2;
    }

    template<int Lt>
    static SymmetryGroup<Lt,double> SymReflectionY()
    {
        auto Refl=TensorPow<Lt/2,2,ElementaryOp<Lt/2>> ( ReflectionOpY<Lt/2> );
        return Z2_Group<Lt>(Refl);
    }

    template<int L>
    static int ReflectionOpX(FockState<L>& f)
    {
        std::string s=f.to_string();
        std::reverse(s.begin(), s.end());
        f=FockState<L>(s);
        if ( (f.count()/2) & 1) return -1;
        else return 1;
    }

    template<int Lt>
    static SymmetryGroup<Lt,double> SymReflectionX()
    {
        auto Refl=TensorPow<Lt/2,2,ElementaryOp<Lt/2>> ( ReflectionOpX<Lt/2> );
        return Z2_Group<Lt>(Refl);
    }

    template<int Lt>
    static SymmetryGroup<Lt,double> SymSpinFlip()
    {
        auto T1=TranslationOp<Lt>(Lt/2);
        return Z2_Group<Lt>(T1);
    }

};

#endif // CADENITAAA5OPEN_H
