#include "all.h"
#include"exactdiagg/all.h"
#include<string>

#define NQUBIT_HUBBARD 20

using namespace std;

QOperator HubbardHam(int L,double U,bool periodic=false)
{
    QOperator H;
    for(int i=0;i<L-1+periodic;i++)
    {
        H.Add(Create(i        ) * Destroy((i+1)%L)   ,1);
        H.Add(Create((i+1)%L  ) * Destroy(i)         ,1);
        H.Add(Create(i+L      ) * Destroy((i+1)%L+L) ,1);
        H.Add(Create((i+1)%L+L) * Destroy(i+L)       ,1);
    }
    for(int i=0;i<L;i++)
        H.Add(Create(i  ) * Destroy(i  ) *
              Create(i+L) * Destroy(i+L) ,U);
    return H;
}

template<int Lt>
SymmetryGroup<Lt,cmpx> HubbardGroup()
{
    const int L=Lt/2;
    auto Tl=TranslationOp<Lt>(L);
    auto Gsf=Z2_Group<Lt,double>( Tl );
    auto Refl=TensorPow<L,2,ElementaryOp<L>> ( ReflectionOp<L> );
    auto Gr=Z2_Group<Lt,double>( Refl );
    auto Geh=Z2_Group<Lt,double> ( ParticleHoleOp<Lt> );

    auto T1=TranslationOp<L>(1);
    auto Trasl=TensorPow<L,2>( T1 );
    auto Gt=CyclicGroupPow<Lt>(Trasl, L);

    auto G=Gt.DirectProd(Geh).DirectProd(Gr).DirectProd(Gsf);
    return G;
}

template<int Lt>
struct HasSz
{
    int sz;
    HasSz(int sz):sz(sz){}
    bool operator()(const FockState<Lt>& f)
    {
        int sum=0;
        for(uint i=0;i<Lt;i++)
            if (i<Lt/2) sum+=f.test(i);
            else sum-=f.test(i);
        return sum==sz;
    }
};

void TestHubbardBasic(double U)
{
    const int Lt=NQUBIT_HUBBARD;
    const int L=Lt/2;
    int nPart=L;

    auto H=HubbardHam(L,U,true);
    auto G=HubbardGroup<Lt>();
    auto b=FockBasisFixedChargeG<Lt>(nPart,G,0,HasSz<Lt>{0});

    auto gs=FindGS<Lt>(H,b,G,nPart);
    b.SetSym(G,1);
    for(int nu=1;nu<2/*G.nSym()*/;nu++,b.SetSym(G,nu))
    {
        auto gsn=FindGS<Lt>(H,b,G,nPart);
        gs=std::min( gs, gsn );
    }
    cout<<"GS: nPart="<<gs.nPart<<" sym="<<gs.sym<<" ener="<<gs.ener;
    b.SetSym(G,gs.sym);
    QOperator nd=Create(0)*Destroy(0)*Create(L)*Destroy(L);
    auto ndv=cdot(gs.state,  nd.toMatrix(b,G,b) * gs.state);
    cout<<"\n"<<U<<" "<<real(ndv)<<" "<<gs.ener<<endl;
    gs.Save(string("gsU")+to_string(U)+".dat");
    gs.Print(b);
}

void TestHubbardQuench(double U1, double U2,double dt, int nt)
{
    EigenStateG<cmpx> gs;
    bool ok=gs.Load(string("gsU")+to_string(U1)+".dat");
    if (!ok) throw invalid_argument("TestHubbardQuench: gs not found");
    cout<<"\nGS: nPart="<<gs.nPart<<" sym="<<gs.sym<<" ener="<<gs.ener<<endl;

    const int Lt=NQUBIT_HUBBARD;
    const int L=Lt/2;
    int nPart=L;

    auto ham=HubbardHam(L,U2,true);
    auto G=HubbardGroup<Lt>();
    auto b=FockBasisFixedChargeG<Lt>(nPart,G,gs.sym,HasSz<Lt>{0});
    auto H=ham.toMatrix(b,G,b);
    auto sol=TimeEvolution(H,gs.state);

    // operadores a medir
    QOperator O=Create(0)*Destroy(0)*Create(L)*Destroy(L); //doble ocupacion
    auto nd_mat=O.toMatrix(b,G,b);
    vector<sp_cx_mat> cic0_mat(L); //operadores ci^dag c0
    for(int i=0;i<L;i++)
    {
        QOperator O=Create(i)*Destroy(0);
        cic0_mat[i]=O.toMatrix(b,G,b);
    }

    vec nd(nt);
    mat cic0(nt,L); // valores de espectacion de los operadores

    for(int i=0;i<nt;i++)
    {
        auto gt=sol.EvolveTo(i*dt);
        nd(i)=real( cdot(gt,nd_mat*gt) );
        for(int p=0;p<L;p++)
            cic0(i,p)=real( cdot(gt,cic0_mat[p]*gt) );
    }
    vec nd_w=square(abs(fft(nd)));
    mat nk_t=square(abs(fft(cic0.t()))).t();
    mat nwk=square(abs(fft2(cic0)));  //Fourier 2d
    nd.save("nd_t.dat",raw_ascii);
    nd_w.save("nd_w.dat",raw_ascii);
    cic0.save("cic0_t.dat",raw_ascii);
    nk_t.save("nk_t.dat",raw_ascii);
    nwk.save("nwk.dat",raw_ascii);
}
