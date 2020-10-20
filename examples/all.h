#ifndef TOROU_H
#define TOROU_H

void TestHamiltonianTB();
void TestHubbardBasic(double U);
void TestHubbardQuench(double U1, double U2,double dt, int nt);
void TestHubbard2Basic(double U, double tp);
void TestHubbard2V(double U, double V, double tp);
void Hernan_ED_gs(const char filename[]);

#endif // TOROU_H
