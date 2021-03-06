#ifndef CADENITAAA3_H
#define CADENITAAA3_H

#include"exactdiagg/all.h"
#include<string>
#include<fstream>
#include<iostream>

struct CadenitaAA3
{
    const int nOrb=3;
    int L; //cantidad de celdas unidad
    double t=1.1, tp=0.1, tpp=0.1, tOO=-0.32, U=2, U3=7, J=0.838, muCu=0.5; //2.595;
    bool periodic=true;
    double mu=0,phi=0,angle_spin=0;

    CadenitaAA3(int L)
        :L(L)
        ,delta(nOrb,nOrb,fill::zeros)
        ,hop(nOrb,nOrb,fill::zeros)
        ,Umat(2*nOrb,2*nOrb,fill::zeros)
    {}
    void Initialize()
    {
      //delta de 0 a 2:
        delta(1,2)=delta(2,1)= t; // t_CuO
        delta(0,2)=delta(2,0)= tp; // t'
        delta(2,2)=-muCu;
        for(int i=0;i<nOrb;i++) delta(i,i)-=mu;
        hop(0,2)= t;
        hop(1,2)= tp;
        hop(2,2)=tpp;
        hop(0,0)=hop(0,1)=hop(1,0)=hop(1,1)=tOO;
        for(int i=0;i<nOrb;i++)
            Umat(toInt(i,0),toInt(i,1))= (i==2) ? U3 : U; //si i=2, elem de matriz = U3, sino =U
        for(int s=0;s<2;s++)
            for(int sp=0;sp<2;sp++)
            {
                double c= (s==sp)? U-3*J : U-2*J;
                Umat(toInt(0,s),toInt(1,sp))=c;
            }
        double q=angle_spin;
        std::vector<double> dat={cos(q/2),sin(q/2),-sin(q/2),cos(q/2)}; // col major
        mat rot(dat.data(),2,2);
        hop=kron(hop,rot);
    }

    int toInt(int i, int Ii, int spin) const { return spin+Ii*2+i*nOrb*2; }
//    int toInt(int i, int Ii, int spin) const { return i+Ii*L+spin*L*nOrb; }
//    int toInt(int i, int Ii, int spin) const { return Ii+i*nOrb+spin*nOrb*L; }
    int toInt(int Ii, int spin) const { return spin+Ii*2; }
//    int toInt(int Ii, int spin) const { return Ii+spin*nOrb; }
    QOperatorG<cmpx> Kin() const
    {
        QOperatorG<cmpx> h;
        for(int i=0;i<L-1+periodic; i++)
            for(int ii=0;ii<nOrb;ii++)
                for(int jj=0;jj<nOrb;jj++)
                {
                    for(int s=0;s<2;s++)
                        for(int sp=0;sp<2;sp++)
                        {
                            double tt=hop(toInt(ii,s), toInt(jj,sp));
                            if (tt!=0)
                            {
                                int pi=toInt(i,ii,s);
                                int pj=toInt((i+1)%L,jj,sp);
                                h.Add( FermiOp(pi,true)*FermiOp(pj,false), tt*exp( cmpx{0,phi}) );
                                h.Add( FermiOp(pj,true)*FermiOp(pi,false), tt*exp(-cmpx{0,phi}) );
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
    QOperatorG<cmpx> Pot() const
    {
        QOperatorG<cmpx> h;
        for(int i=0;i<L; i++)
            for(int ii=0;ii<nOrb;ii++)
                for(int s=0;s<2;s++)
                    for(int jj=0;jj<nOrb;jj++)
                        for(int sp=0;sp<2;sp++)
                        {
                            int pi=toInt(i,ii,s);
                            int pj=toInt(i,jj,sp);
                            double coeff=Umat(toInt(ii,s),toInt(jj,sp));
                            if (coeff==0) continue;
                            h.Add( FermiOp(pi,true)*FermiOp(pi,false)*           //U n n
                                   FermiOp(pj,true)*FermiOp(pj,false), coeff);
                        }
        if (J!=0)
            for(int i=0;i<L; i++)
                for(int ii=0;ii<2;ii++)
                {
                    int pi=toInt(i,ii,0);
                    int pj=toInt(i,1-ii,1);
                    int pk=toInt(i,1-ii,0);
                    int pl=toInt(i,ii,1);
                    h.Add( FermiOp(pi,true)*FermiOp(pj,true)*           // -J c+c+ cc
                           FermiOp(pk,false)*FermiOp(pl,false), -J );
                    pi=toInt(i,ii,0);
                    pj=toInt(i,ii,1);
                    pk=toInt(i,1-ii,1);
                    pl=toInt(i,1-ii,0);
                    h.Add( FermiOp(pi,true)*FermiOp(pj,true)*           // J c+c+ cc
                           FermiOp(pk,false)*FermiOp(pl,false), J );
                }
        return h;
    }
    QOperatorG<cmpx> Ham() const { return Kin()+Pot(); }

private:
    arma::mat delta, hop, Umat;
};



#endif // CADENITAAA3_H
