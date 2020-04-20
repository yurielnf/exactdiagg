# exactdiagg
Exact diagonalization to the quantum many-body problem using group theory and the Fock basis.
Each Fock state is represented as a bitset. To generate the sparse matrix representation of a quantum operator O, when O is applied to a bitset, the resulting Fock states are indexed in the basis using hash table. The symmetry group G splits the original problem into #G independent and smaller problems. A basis of representants is used for each irrep of G while the larger full Fock basis is never used.

Automatic Hamiltonian construction as in:

```
h.Add(Create(i)*Destroy((i+1)%L),-t);
H.Add(Create(i)*Destroy(i) * Create(i+L)*Destroy(i+L) ,U);
```

and symmetry group:

```
auto Gr=Z2_Group<L>( ReflectionOp<L> );
auto Geh=Z2_Group<L> ( ParticleHoleOp<L> );
auto T1=TranslationOp<L>(1);
auto Gt=CyclicGroupPow<L>(T1, L);

auto G=Gr.DirectProd(Geh);
```

### External dependencies

- [armadillo](https://gitlab.com/conradsnicta/armadillo-code)
- lapack
- blas
- arpack

### Getting started

This project can be edited and build using the [`qtcreator`](https://github.com/qt-creator) IDE.
